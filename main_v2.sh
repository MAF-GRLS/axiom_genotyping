mkdir -p MAF && cd MAF
work_dir=$(pwd)

## Prepare genomic resources
## get the referenece genome
genome_dir=$work_dir/refGenome
mkdir -p $genome_dir && cd $genome_dir
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/canFam3/bigZips/canFam3.fa.gz' -O canFam3.fa.gz
gunzip canFam3.fa.gz
# chromosome Y was identified by blasting the flanking sequences. The target sequence was obtained from the GenBank: https://www.ncbi.nlm.nih.gov/nuccore/KP081776.1
# Local: cp CF_y.fasta tahmed@farm.cse.ucdavis.edu:/home/tahmed/MAF/refGenome/.
(echo ">chrY" && tail -n+2 CF_y.fasta) >> canFam3.fa
canFam3_ref="$genome_dir"/canFam3.fa

cat $canFam3_ref | awk '/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}' > $genome_dir/canFam3_unwrap.fa
canFam3_ref_unwrap=$genome_dir/canFam3_unwrap.fa

conda create -n ngs
conda activate ngs
conda install -c bioconda bwa

mkdir -p $genome_dir/bwaIndex && cd $genome_dir/bwaIndex
ln -s $canFam3_ref_unwrap .
bwa index -a bwtsw canFam3_unwrap.fa

## Download genotyping CEL files from Google cloud
#You can access the bucket via this URL:
#https://console.cloud.google.com/storage/browser/entirecohortinfo.
conda update conda
conda create -n gsutil_env -c conda-forge gsutil
conda activate gsutil_env
gsutil config ## https://stackoverflow.com/questions/49302859/gsutil-serviceexception-401-anonymous-caller-does-not-have-storage-objects-list
gsutil -m cp -r "gs://entirecohortinfo/CEL_Files" "gs://entirecohortinfo/Entire Cohort Zip Files" .

## remove the spaces from the working paths
mv "Entire Cohort Zip Files" Entire_Cohort_Zip_Files
mv "Entire_Cohort_Zip_Files/Array_A/Raw data" Entire_Cohort_Zip_Files/Array_A/Raw_data
mv "Entire_Cohort_Zip_Files/Array_B/Raw data" Entire_Cohort_Zip_Files/Array_B/Raw_data
mv "CEL_Files/Array A" CEL_Files/Array_A
mv "CEL_Files/Array B" CEL_Files/Array_B

## How much data do we have?
cd $work_dir/
ls -1 CEL_Files/Array_A/ | wc -l ## 3168 samples
ls -1 CEL_Files/Array_A/ | cut -d"-" -f4 | cut -d "_" -f1 | uniq -c | wc -l ## 35 96-well plates
ls -1 CEL_Files/Array_A/ | cut -d"-" -f2,3 | uniq -c ## 4 batches (I am not sure if this prefix indicate batch)

ls -1 CEL_Files/Array_B/ | wc -l ## 3168 samples
ls -1 CEL_Files/Array_B/ | cut -d"-" -f4 | cut -d "_" -f1 | uniq -c | wc -l ## 41 96-well plates. Some plates have very few samples (e.g. plate_667 has 1 sample, plate_489 has 2 samples, plate_678 & plate_679 have 3 samples)
ls -1 CEL_Files/Array_B/ | cut -d"-" -f2,3 | uniq -c ## 7 batches including 2 sample batches < 96 samples

## unify the cel file names
ls -1 CEL_Files/Array_A/ | cut -d"-" -f4 | awk -F"[_\.]" '{if($3!="CEL")print}' ## 437_H02_2.CEL
grep 437_H02 <(ls -1 CEL_Files/Array_A/*.CEL) ## it has only one version so I will rename
mv CEL_Files/Array_A/a550771-4434582-040323-437_H02_2.CEL CEL_Files/Array_A/a550771-4434582-040323-437_H02.CEL

ls -1 CEL_Files/Array_B/ | cut -d"-" -f4 | awk -F"[_\.]" '{if($3!="CEL")print}' ## 658_E07_2.CEL && 062_C03_2.CEL
grep 658_E07 <(ls -1 CEL_Files/Array_B/*.CEL) ## it has only one version so I will rename
grep 062_C03 <(ls -1 CEL_Files/Array_B/*.CEL) ## it has only one version so I will rename
mv CEL_Files/Array_B/a550772-4435950-041523-658_E07_2.CEL CEL_Files/Array_B/a550772-4435950-041523-658_E07.CEL
mv CEL_Files/Array_B/a550772-4435952-041623-062_C03_2.CEL CEL_Files/Array_B/a550772-4435952-041623-062_C03.CEL

## Add the pilot data
gsutil -m cp -r "gs://pilot-180-genotype/original-zip-files/" .
mkdir -p original-zip-files/{Array_A,Array_B}
mv original-zip-files/{a550771-4409629-062221-343.zip,a550771-4409629-062221-344.zip} original-zip-files/Array_A/.
mv original-zip-files/{a550772-4409630-070221-050.zip,a550772-4409630-070221-051.zip} original-zip-files/Array_B/.

cd $work_dir/original-zip-files/Array_A/
unzip a550771-4409629-062221-343.zip
unzip a550771-4409629-062221-344.zip
cd $work_dir/original-zip-files/Array_B/
unzip a550772-4409630-070221-050.zip
unzip a550772-4409630-070221-051.zip

## download NEW support file (Analysis library files). We got directly from Thermo
module load rclone
mkdir -p $work_dir/lib2/{setA,setB}
rclone -v --copy-links copy remote_UCDavis_GoogleDr:MAF/New_axiom_files/Axiom_K9HDSNPA_Analysis.r2.zip /home/tahmed/MAF/lib2/setA/
rclone -v --copy-links copy remote_UCDavis_GoogleDr:MAF/New_axiom_files/Axiom_K9HDSNPA.na35.r2.a4.annot.db /home/tahmed/MAF/lib2/setA/

rclone -v --copy-links copy remote_UCDavis_GoogleDr:MAF/New_axiom_files/Axiom_K9HDSNPB_Analysis.r2.zip /home/tahmed/MAF/lib2/setB/
rclone -v --copy-links copy remote_UCDavis_GoogleDr:MAF/New_axiom_files/Axiom_K9HDSNPB.na35.r2.a4.annot.db /home/tahmed/MAF/lib2/setB/

cd $work_dir/lib2/setA
unzip -d setA_Analysis Axiom_K9HDSNPA_Analysis.r2.zip
cd $work_dir/lib2/setB
unzip -d setB_Analysis Axiom_K9HDSNPB_Analysis.r2.zip

# extract annotation csv
mamba create -n sqlite -c conda-forge sqlite
conda activate sqlite
## list the tables
sqlite3 $work_dir/lib2/setA/Axiom_K9HDSNPA.na35.r2.a4.annot.db ".tables"
#Annotations     CdfInformation  Chromosome      Information     Localization
## Extract traget tables
cd $work_dir/lib2/setA
#sqlite3 -header -csv Axiom_K9HDSNPA.na35.r2.a4.annot.db "SELECT * FROM Annotations;" > Axiom_K9HDSNPA_Annotation.r2_3.csv
sqlite3 -header -tabs Axiom_K9HDSNPA.na35.r2.a4.annot.db "SELECT * FROM Annotations;" > Axiom_K9HDSNPA_Annotation.r2_3.tab
sqlite3 -header -tabs Axiom_K9HDSNPA.na35.r2.a4.annot.db "SELECT * FROM Chromosome;" > Axiom_K9HDSNPA_Annotation.r2_3.chr.tab
sqlite3 -header -tabs Axiom_K9HDSNPA.na35.r2.a4.annot.db "SELECT * FROM Localization;" > Axiom_K9HDSNPA_Annotation.r2_3.loc.tab
sqlite3 -header -tabs Axiom_K9HDSNPA.na35.r2.a4.annot.db "SELECT * FROM Information;" > Axiom_K9HDSNPA_Annotation.r2_3.inf.tab


cd $work_dir/lib2/setB
#sqlite3 -header -csv Axiom_K9HDSNPB.na35.r2.a4.annot.db "SELECT * FROM Annotations;" > Axiom_K9HDSNPB_Annotation.r2_3.csv
sqlite3 -header -tabs Axiom_K9HDSNPB.na35.r2.a4.annot.db "SELECT * FROM Annotations;" > Axiom_K9HDSNPB_Annotation.r2_3.tab
sqlite3 -header -tabs Axiom_K9HDSNPB.na35.r2.a4.annot.db "SELECT * FROM Chromosome;" > Axiom_K9HDSNPB_Annotation.r2_3.chr.tab
sqlite3 -header -tabs Axiom_K9HDSNPB.na35.r2.a4.annot.db "SELECT * FROM Localization;" > Axiom_K9HDSNPB_Annotation.r2_3.loc.tab
sqlite3 -header -tabs Axiom_K9HDSNPB.na35.r2.a4.annot.db "SELECT * FROM Information;" > Axiom_K9HDSNPB_Annotation.r2_3.inf.tab

####
## Pedigree ped files were obtained from Slack. Therefore, I had to download locally then upload to the server
mkdir -p $work_dir/pedigree && cd $work_dir/pedigree
# local: scp sample_ped_infoA.txt tahmed@farm.cse.ucdavis.edu:/home/tahmed/MAF/pedigree/.
# local: scp sample_ped_infoB.txt tahmed@farm.cse.ucdavis.edu:/home/tahmed/MAF/pedigree/.

## Explore the files
tail -n+2 sample_ped_infoA.txt | cut -f2 | sort | uniq -c  ## 2968 GRLS & 200 Oldies
tail -n+2 sample_ped_infoB.txt | cut -f2 | sort | uniq -c  ## 2968 GRLS & 200 Oldies & GRLS Pilot

tail -n+2 sample_ped_infoA.txt | cut -f3 | sed 's/.*-//' | sort > tempA_ids
tail -n+2 sample_ped_infoB.txt | cut -f3 | sed 's/.*-//' | sort > tempB_ids
cat tempA_ids | uniq > tempA_uids
cat tempB_ids | uniq > tempB_uids
wc -l temp*ids ## B has replicates & A and B has different number of IDs
#    3168 tempA_ids
#    3168 tempA_uids
#    3360 tempB_ids
#    3358 tempB_uids
cat tempB_ids | uniq -c | sort -nr | head ## 027376 & 003698 have duplicates
comm -12 tempA_uids tempB_uids | wc -l #    3166
comm -13 tempA_uids tempB_uids | wc -l #    192
comm -23 tempA_uids tempB_uids #    027376_1 & 027376_2

grep 027376 sample_ped_infoA.txt sample_ped_infoB.txt
#sample_ped_infoA.txt:a550771-4434581-040223-137_D08.CEL    GRLS    094-027376_1    0    0    1    -9
#sample_ped_infoA.txt:a550771-4434581-040223-137_E08.CEL    GRLS    094-027376_2    0    0    1    -9
#sample_ped_infoB.txt:a550772-4435950-041523-650_D08.CEL    GRLS    094-027376    0    0    1    -9
#sample_ped_infoB.txt:a550772-4435950-041523-650_E08.CEL    GRLS    094-027376    0    0    1    -9

grep 003698 sample_ped_infoA.txt sample_ped_infoB.txt
#sample_ped_infoA.txt:a550771-4434581-040223-129_F12.CEL    GRLS    094-003698    0    0    1    -9
#sample_ped_infoB.txt:a550772-4409630-070221-050_E10.CEL    GRLS Pilot    003698    0    0    1    -9
#sample_ped_infoB.txt:a550772-4448020-102423-494_C01.CEL    GRLS    094-003698    0    0    1    -9

## Edit the info files to correct CEL file names with "_2" and to remove the "094-" prefix
sed 's/_2\.CEL/\.CEL/' sample_ped_infoA.txt > temp_infoA.txt
sed 's/_2\.CEL/\.CEL/' sample_ped_infoB.txt | sed 's/\(650_D08.*094-027376\)/\1_1/; s/\(650_E08.*094-027376\)/\1_2/;' > temp_infoB.txt
(head -n1 sample_ped_infoA.txt;ls -1 $work_dir/CEL_Files/Array_A/*.CEL | sed 's|^.*/||' | grep -Fwf - temp_infoA.txt) | awk 'BEGIN{FS=OFS="\t"}{gsub(/^094-/,"",$3);print}' > sample_ped_infoA_cohort.txt
(head -n1 sample_ped_infoB.txt;ls -1 $work_dir/CEL_Files/Array_B/*.CEL | sed 's|^.*/||' | grep -Fwf - temp_infoB.txt) | awk 'BEGIN{FS=OFS="\t"}{gsub(/^094-/,"",$3);print}' > sample_ped_infoB_cohort.txt
rm temp_infoA.txt temp_infoB.txt

## confirm the 2 files are matched
cat sample_ped_infoA_cohort.txt | cut -f3 | sort > a3.txt
cat sample_ped_infoB_cohort.txt | cut -f3 | sort > b3.txt
diff a3.txt b3.txt
rm a3.txt b3.txt

cat sample_ped_infoA_cohort.txt | cut -f2,3 | sort > a2.txt
cat sample_ped_infoB_cohort.txt | cut -f2,3 | sort > b2.txt
diff a2.txt b2.txt
rm a2.txt b2.txt
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$3]=$2;next}{if(a[$3]!=$2)print $3,a[$3],$2}' sample_ped_infoA_cohort.txt sample_ped_infoB_cohort.txt

## Family IDs are either “GRLS or Oldies”. We have 2968 GRLS & 200 Oldies in both Arrays. Each sample – as expected – has the same Family ID in both files except for 14 samples that seem to be switched. The meta-data files were updated based on the provided correct IDs
mv sample_ped_infoA_cohort.txt sample_ped_infoA_cohort.txt_wrongFamIDs
echo "040039 GRLS Oldies
040236 GRLS Oldies
040176 GRLS Oldies
040082 GRLS Oldies
040052 GRLS Oldies
040006 GRLS Oldies
040089 GRLS Oldies
003430 Oldies GRLS
000686 Oldies GRLS
015509 Oldies GRLS
030301 Oldies GRLS
033913 Oldies GRLS
030026 Oldies GRLS
004582 Oldies GRLS" | tr ' ' '\t' > famIDs_correct
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2;b[$1]=$3;next}{if(a[$3]&&a[$3]==$2)$2=b[$3];print}' famIDs_correct sample_ped_infoA_cohort.txt_wrongFamIDs > sample_ped_infoA_cohort.txt


## find sample name duplicates (We can use any one of the arrays ped files)
tail -n+2 sample_ped_infoA_cohort.txt | cut -f3 | cut -d"_" -f1 | sort | uniq -c | awk '{if($1>1)print $2}' | while read f;do echo $(grep "$f" sample_ped_infoA_cohort.txt | awk -F"\t" '{print "S"$3}' | tr '\n' '\t');done > sample_dup.txt


## compare gender across info files
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$3]=$6;next}{print $3,a[$3],$6}' sample_ped_infoA_cohort.txt sample_ped_infoB_cohort.txt | sed 's/Sex/sexA/' | sed 's/Sex/sexB/' > gender_compare.txt
tail -n+2 gender_compare.txt | awk 'BEGIN{FS=OFS="\t"}{if($2!=$3)print}'
tail -n+2 gender_compare.txt | cut -f2 | sort | uniq -c  ##    1575 1(male) // 1587 2(female) // 6 -9(unknown)

####
## Pedigree ped files for pilot data. Also, obtained from Brenna. Therefore, I had to download locally then upload to the server
mkdir -p $work_dir/ped_pilot && cd $work_dir/ped_pilot
# local: scp sample_ped_infoA.txt tahmed@farm.cse.ucdavis.edu:/home/tahmed/MAF/ped_pilot/.
# local: scp sample_ped_infoB.txt tahmed@farm.cse.ucdavis.edu:/home/tahmed/MAF/ped_pilot/.

## Explore the files
tail -n+2 sample_ped_infoA.txt | cut -f2 | sort | uniq -c  ## 192 PILOT_A
tail -n+2 sample_ped_infoB.txt | cut -f2 | sort | uniq -c  ## 192 PILOT_B

tail -n+2 sample_ped_infoA.txt | cut -f3 | sed 's/.*-//' | sort > tempA_ids
tail -n+2 sample_ped_infoB.txt | cut -f3 | sed 's/.*-//' | sort > tempB_ids
cat tempA_ids | uniq > tempA_uids
cat tempB_ids | uniq > tempB_uids
wc -l temp*ids ## No replicate IDs
#    192 tempA_ids
#    192 tempA_uids
#    192 tempB_ids
#    192 tempB_uids
comm -12 tempA_uids tempB_uids | wc -l #    192

## Check for replicate IDs with the full cohort
cat <(tail -n+2 ../pedigree/sample_ped_infoA_cohort.txt) <(tail -n+2 sample_ped_infoA_cohort.txt) | cut -f3 | sort | uniq -c | sort -k1,1nr | head


## unify the family ids & update sample “003698”
cat sample_ped_infoA.txt | sed 's/PILOT_A/PILOT/' | sed 's/003698/003698_pilot/' > sample_ped_infoA_cohort.txt
cat sample_ped_infoB.txt | sed 's/PILOT_B/PILOT/' | sed 's/003698/003698_pilot/' > sample_ped_infoB_cohort.txt

## confirm the 2 files are matched
cat sample_ped_infoA_cohort.txt | cut -f3 | sort > a3.txt
cat sample_ped_infoB_cohort.txt | cut -f3 | sort > b3.txt
diff a3.txt b3.txt
rm a3.txt b3.txt

cat sample_ped_infoA_cohort.txt | cut -f2,3 | sort > a2.txt
cat sample_ped_infoB_cohort.txt | cut -f2,3 | sort > b2.txt
diff a2.txt b2.txt
rm a2.txt b2.txt
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$3]=$2;next}{if(a[$3]!=$2)print $3,a[$3],$2}' sample_ped_infoA_cohort.txt sample_ped_infoB_cohort.txt

## find sample name duplicates (We can use any one of the arrays ped files)
tail -n+2 sample_ped_infoA_cohort.txt | cut -f3 | sed 's/A//;s/B//;' | sort | uniq -c | awk '{if($1>1)print $2}' | while read f;do echo $(grep "$f" sample_ped_infoA_cohort.txt | awk -F"\t" '{print "S"$3}' | tr '\n' '\t');done > sample_dup.txt ## 12

## compare gender across info files
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$3]=$6;next}{print $3,a[$3],$6}' sample_ped_infoA_cohort.txt sample_ped_infoB_cohort.txt | sed 's/Sex/sexA/' | sed 's/Sex/sexB/' > gender_compare.txt
tail -n+2 gender_compare.txt | awk 'BEGIN{FS=OFS="\t"}{if($2!=$3)print}'
tail -n+2 gender_compare.txt | cut -f2 | sort | uniq -c  ## 105 1(male) // 87 2(female)

## Merge full cohort and pilot datasets
mkdir -p $work_dir/All_CEL_Files/{Array_A,Array_B}
cd $work_dir/All_CEL_Files/Array_A
ln -s $work_dir/CEL_Files/Array_A/*.CEL .
ln -s $work_dir/original-zip-files/Array_A/*.CEL .
cd $work_dir/All_CEL_Files/Array_B
ln -s $work_dir/CEL_Files/Array_B/*.CEL .
ln -s $work_dir/original-zip-files/Array_B/*.CEL .

mkdir -p $work_dir/All_pedigree
# Fix the gender code (change -9 to 0) and update individual IDs "Add S before the numerical individual ID"
head -n1 $work_dir/pedigree/sample_ped_infoA_cohort.txt > $work_dir/All_pedigree/sample_ped_infoA_cohort.txt
tail -n+2 $work_dir/pedigree/sample_ped_infoA_cohort.txt | awk 'BEGIN{FS=OFS="\t"}{if($6==-9)$6=0;$3="S"$3;print}' >> $work_dir/All_pedigree/sample_ped_infoA_cohort.txt
tail -n+2 $work_dir/ped_pilot/sample_ped_infoA_cohort.txt | awk 'BEGIN{FS=OFS="\t"}{if($6==-9)$6=0;$3="S"$3;print}' >> $work_dir/All_pedigree/sample_ped_infoA_cohort.txt

head -n1 $work_dir/pedigree/sample_ped_infoB_cohort.txt > $work_dir/All_pedigree/sample_ped_infoB_cohort.txt
tail -n+2 $work_dir/pedigree/sample_ped_infoB_cohort.txt | awk 'BEGIN{FS=OFS="\t"}{if($6==-9)$6=0;$3="S"$3;print}' >> $work_dir/All_pedigree/sample_ped_infoB_cohort.txt
tail -n+2 $work_dir/ped_pilot/sample_ped_infoB_cohort.txt | awk 'BEGIN{FS=OFS="\t"}{if($6==-9)$6=0;$3="S"$3;print}' >> $work_dir/All_pedigree/sample_ped_infoB_cohort.txt

## Upon the communication with Brenna regarding the samples predicted to have a computed gender different from those in the meta-data file, she confirmed that the sample “S040136” is a female no male. The metadata file was updated accordingly
cd $work_dir/All_pedigree
for f in sample_ped_info*_cohort.txt;do
awk 'BEGIN{FS=OFS="\t"}{if($3=="S040136")$6=2;print}' $f > $f.temp;
done
mv sample_ped_infoA_cohort.txt.temp sample_ped_infoA_cohort.txt
mv sample_ped_infoB_cohort.txt.temp sample_ped_infoB_cohort.txt
########
## Command line software
## Applied BiosystemsTM Array Power Tools (APT) && SNPolisher (R package for visualization)

## Install APT software
## https://www.thermofisher.com/ca/en/home/life-science/microarray-analysis/microarray-analysis-partners-programs/affymetrix-developers-network/affymetrix-power-tools.html
## Documentation:
## a) axiom_genotyping_solution_analysis_guide.pdf
## b) The APT documentation (in HTML pages) is included in the APT downloads and available at (http://media.affymetrix.com/support/developer/powertools/changelog/index.html)

wget http://www.thermofisher.com/content/dam/LifeTech/Documents/ZIP/apt_2.11.6_linux_64_x86_binaries.zip
unzip apt_2.11.6_linux_64_x86_binaries.zip
chmod +x apt_2.11.6_linux_64_x86_binaries/bin/*
apt=$work_dir/"apt_2.11.6_linux_64_x86_binaries/bin/"

## Install SNPolisher (Locally)
## https://www.thermofisher.com/ca/en/home/life-science/microarray-analysis/microarray-analysis-partners-programs/affymetrix-developers-network/affymetrix-devnet-tools.html
## Documentation: SNPolisher documentation (in pdf) is included in the download and available at (https://assets.thermofisher.com/TFS-Assets/LSG/manuals/MAN0017790_SNPolisher_UG.pdf)

wget https://downloads.thermofisher.com/SNPolisher_3.0.zip
unzip SNPolisher_3.0.zip  ## The folder has the manual pdf
# Open R
# setwd("/Users/drtamermansour/Source/SNPolisher_3.0")
# install.packages("SNPolisher_3.0.tar.gz",repos=NULL,type="source")

##########
## Define working Paths then run the Best Practice analysis once for Array A then again for Array B
## Array_A
cel_PATH=$work_dir/"All_CEL_Files/Array_A"
set_Analysis=$work_dir/"lib2/setA/setA_Analysis/"
set_output=$work_dir/"output/setA"
qc_xml="$set_Analysis"/Axiom_K9HDSNPA.r2.apt-geno-qc.AxiomQC1.xml
genotype1_xml="$set_Analysis"/Axiom_K9HDSNPA_96orMore_Step1.r2.apt-genotype-axiom.AxiomGT1.apt2.xml
genotype2_xml="$set_Analysis"/Axiom_K9HDSNPA_96orMore_Step2.r2.apt-genotype-axiom.AxiomGT1.apt2.xml
specialSNPs="$set_Analysis"/Axiom_K9HDSNPA.r2.specialSNPs
ps2snp="$set_Analysis"/Axiom_K9HDSNPA.r2.ps2snp_map.ps
annot_db="$set_Analysis"/../Axiom_K9HDSNPA.na35.r2.a4.annot.db
annot_csv="$set_Analysis"/../Axiom_K9HDSNPA_Annotation.r2_3.tab
pedInfo=$work_dir/All_pedigree/sample_ped_infoA_cohort.txt
chr_map="$set_Analysis"/../Axiom_K9HDSNPA_Annotation.r2_3.chr.tab


## Array_B
cel_PATH=$work_dir/"All_CEL_Files/Array_B" #$work_dir/"CEL_Files/Array_A"
set_Analysis=$work_dir/"lib2/setB/setB_Analysis/" #$work_dir/"lib/setA/setA_Analysis/"
set_output=$work_dir/"output/setB" #$work_dir/"output/setA"
qc_xml="$set_Analysis"/Axiom_K9HDSNPB.r2.apt-geno-qc.AxiomQC1.xml #"$set_Analysis"/Axiom_K9HDSNPA.r2.apt-geno-qc.AxiomQC1.xml
genotype1_xml="$set_Analysis"/Axiom_K9HDSNPB_96orMore_Step1.r2.apt-genotype-axiom.AxiomGT1.apt2.xml #"$set_Analysis"/Axiom_K9HDSNPA_96orMore_Step1.r2.apt-genotype-axiom.AxiomGT1.apt2.xml
genotype2_xml="$set_Analysis"/Axiom_K9HDSNPB_96orMore_Step2.r2.apt-genotype-axiom.AxiomGT1.apt2.xml
#"$set_Analysis"/Axiom_K9HDSNPA_96orMore_Step2.r2.apt-genotype-axiom.AxiomGT1.apt2.xml
specialSNPs="$set_Analysis"/Axiom_K9HDSNPB.r2.specialSNPs #"$set_Analysis"/Axiom_K9HDSNPA.r2.specialSNPs
ps2snp="$set_Analysis"/Axiom_K9HDSNPB.r2.ps2snp_map.ps #"$set_Analysis"/Axiom_K9HDSNPA.r2.ps2snp_map.ps
annot_db="$set_Analysis"/../Axiom_K9HDSNPB.na35.r2.a4.annot.db #"$set_Analysis"/../Axiom_K9HDSNPA.na35.r2.a4.annot.db
annot_csv="$set_Analysis"/../Axiom_K9HDSNPB_Annotation.r2_3.tab
pedInfo=$work_dir/All_pedigree/sample_ped_infoB_cohort.txt #$work_dir/pedigree/sample_ped_infoA_cohort.txt
chr_map="$set_Analysis"/../Axiom_K9HDSNPB_Annotation.r2_3.chr.tab

##########
## Start the Best Practice analysis
cd $work_dir
mkdir -p "$set_output"

## 1. Group samples into batches
## Notes: axiom_genotyping_solution_analysis_guide.pdf (20)
##    1. group plates in as large a batch size as is computationally feasible, or up to 4,800 samples unless plates have known significant differences.
##    2. SNP-specific priors should be used when the total batch size is between 20 and 96 unique individuals.
##    3. The specific genotyping option for large (≥96 samples) or small (<96 samples) batch sizes must be selected in all workflows.
##    4. Each batch should contain either 15 or more distinct female samples or zero female samples
## .CEL files corresponding to each batch must be collected into a file with the full path to each .CEL file in each row and with a header line = “cel_files”.
##    5. computational resources: As a reference point, a batch size of 55 AxiomTM 96-Array Plates, each with ∼650K probe sets, requires about 16 hours to execute step 7 (Figure 7) using the apt-genotype-axiom command (“Best practices step 7: Genotype passing samples and plates using AxiomGT1.Step2” on page 97) on a LinuxTM server with the following configuration: x86_64 architecture, 16 x 3 GHz XeonTM core, and 128 GB of RAM. Note that this is without any computational parallelization.
(echo cel_files; ls -1 $cel_PATH/*.CEL) > "$set_output"/cel_list1.txt

## 2. Generate sample DQC values
## DQC is based on intensities of probe sequences for non-polymorphic genome locations
$apt/apt-geno-qc \
    --analysis-files-path "$set_Analysis" \
    --xml-file "$qc_xml" \
    --cel-files "$set_output"/cel_list1.txt \
    --out-file "$set_output"/apt-geno-qc.txt \
    --log-file "$set_output"/apt-geno-qc.log


## 3: Conduct sample QC on DQC
#DQC is based on intensities of probe sequences for non-polymorphic genome locations
#refer to the column "axiom_dishqc_DQC" in the file <OUTDIR>/apt-geno-qc.txt
#remove .CELs from the cel_list1.txt with DQC values that are < 0.82. We refer to filtered .CEL list from this step as cel_list2.txt.
(echo cel_files; awk -v p="$cel_PATH/" -F"\t" 'FNR==NR{if($18>0.82)a[p$1]=1;next}{if(a[$1])print}' <(grep -v "^#"  "$set_output"/apt-geno-qc.txt) "$set_output"/cel_list1.txt) > "$set_output"/cel_list2.txt  ## arrayA (0) // arrayB (650_E12.CEL and 658_B12.CEL which have values 0.055 and 0.092)
 

## 4: Generate sample QC call rates with APT
$apt/apt-genotype-axiom \
    --log-file "$set_output"/apt-genotype-axiom.log \
    --arg-file "$genotype1_xml" \
    --analysis-files-path "$set_Analysis" \
    --out-dir "$set_output"/step1 \
    --dual-channel-normalization true \
    --table-output false \
    --cel-files "$set_output"/cel_list2.txt

## 5: QC the samples based on QC call rate in APT
# Remove samples with a QC call rate value less than the default threshold of 97%. To execute this filter step, refer to the column "call_rate" in the file "<OUTDIR>/step1/AxiomGT1.report.txt"
# remove .CELs from the cel_list2.txt whose call rate values are less than 97%. We refer to this .CEL list as cel_list3.txt.
## Explore
grep -v "^#"  "$set_output"/step1/AxiomGT1.report.txt | tail -n+2 | cut -f1,3 | cut -d"-" -f4 | awk -F"[_\t]" '{if($3<97)print $1,$3}' ## Array_A (8) // Array_B (42)
grep -v "^#"  "$set_output"/step1/AxiomGT1.report.txt | tail -n+2 | cut -f1,3 | cut -d"-" -f4 | awk -F"[_\t]" '{if($3<95)print $1,$3}' ## Array_A (5) // Array_B (19)
grep -v "^#"  "$set_output"/step1/AxiomGT1.report.txt | tail -n+2 | cut -f1,3 | cut -d"-" -f4 | awk -F"[_\t]" '{if($3<95)print $1}' | sort | uniq -c ## Array_A (131,438x2, 137x1) // Array_B (655x6, 493x3, 494,496,650x2, 491,495,677,679x1)
grep -v "^#"  "$set_output"/step1/AxiomGT1.report.txt | tail -n+2 | cut -f1,3 | cut -d"-" -f4 | awk -F"[_\t]" '{if($3>95 && $3<97)print $1,$3}' ## Array_A (129,131,134x1) // Array_B (493x5, 495x4, 679,066,492,497x2, 658,666,064,070,151,496x1) => These are the plate that will likely fail the average plate QC
# Tamer: I will use 95% as a cutoff and accordingly I will tolerate lower average plate QC
(echo cel_files; awk -v p="$cel_PATH/" -F"\t" 'FNR==NR{if($3>95)a[p$1]=1;next}{if(a[$1])print}' <(grep -v "^#"  "$set_output"/step1/AxiomGT1.report.txt) "$set_output"/cel_list1.txt) > "$set_output"/cel_list3.txt

## 6: QC the plates
tail -n+2 "$set_output"/cel_list1.txt | sed 's|.*/||' | cut -d"-" -f4 | cut -d "_" -f1 | uniq -c | sed 's/^ *[^0-9]//g' > "$set_output"/samples_on_plate
tail -n+2 "$set_output"/cel_list3.txt | sed 's|.*/||' | cut -d"-" -f4 | cut -d "_" -f1 | uniq -c | sed 's/^ *[^0-9]//g' > "$set_output"/samples_pass_plate
awk 'BEGIN{FS=" "}FNR==NR{a[$2]=$1;next}{if(a[$2])print $2,$1/a[$2] * 100}' "$set_output"/samples_on_plate "$set_output"/samples_pass_plate  > "$set_output"/plate_pass_rate ## As expected, 3 plates on Array_A (131, 137, 438) have failing samples but all are above 97.5 // Array_B (9 plates have lost samples)
grep -v "^#"  "$set_output"/step1/AxiomGT1.report.txt | tail -n+2 | cut -f1,3 | cut -d"-" -f4 | awk -F"[_\t]" 'BEGIN{CONVFMT="%.3g"}{if($3>95){sumQC[$1]+=$3;count[$1]++}}END{for(i in sumQC)print i,sumQC[i]/count[i]}' > "$set_output"/plate_aveQC ## Array_A: 129 (ave=98.3873).  // Array_B: 679 (ave=95.8125), 667 (ave=97.675), 678 (ave=97.775)

cp "$set_output"/cel_list3.txt "$set_output"/cel_list4.txt

##  7c: Genotype passing samples and plates using AxiomGT1.Step2
$apt/apt-genotype-axiom \
    --log-file "$set_output"/step2/apt-genotype-axiom.log \
    --arg-file "$genotype2_xml" \
    --analysis-files-path "$set_Analysis" \
    --out-dir "$set_output"/step2 \
    --dual-channel-normalization true \
    --allele-summaries true \
    --genotyping-node:snp-posteriors-output true \
    --batch-folder "$set_output"/suitefiles \
    --cel-files "$set_output"/cel_list4.txt

## check the accuracy of computed gender
head -n1 $pedInfo | awk 'BEGIN{FS=OFS="\t"}{print $1,$3,$6,"computed_gender"}' > "$set_output"/step2/gender_check
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2;next}{if($6=="1")$6="male";else if($6=="2")$6="female";else $6="unknown"; if(a[$1])print $1,$3,$6,a[$1]}' <(grep -v "^#" "$set_output"/step2/AxiomGT1.report.txt | tail -n+2) <(tail -n+2 $pedInfo) >> "$set_output"/step2/gender_check
cat "$set_output"/step2/gender_check | awk 'BEGIN{FS=OFS="\t"}{if($3!="unknown" && $4=="unknown")print}' | wc -l #1A/0B
cat "$set_output"/step2/gender_check | awk 'BEGIN{FS=OFS="\t"}{if($3=="unknown" && $4=="unknown")print}' | wc -l #0A/0B
head -n1 "$set_output"/step2/gender_check > "$set_output"/step2/gender_dif
tail -n+2 "$set_output"/step2/gender_check | awk 'BEGIN{FS=OFS="\t"}{if(($3=="unknown") || ($4!="unknown" && $3!=$4))print}' | sort -k2,2 | tr '\t' ',' >> "$set_output"/step2/gender_dif

## Note: Creating smaller genotyping output files
## Ps-extract extracts data from calls, confidences, summary, references, posteriors (biallelic and multiallelic), and priors (biallelic and multiallelic) files for a supplied set of probesets and/or samples. Every input file has a matching output file.
## This could be VERY useful stratgy to prep files for visualization

## 8A: Run ps-metrics
# The ps-metrics function calculates 17 basic SNP QC metrics for each probeset
# Ps-metrics can be run on a subset of probesets or samples.
# Notes: 1) running without report-file cause warning. 2) adding report-file cause fewer recommended SNPs. 3) report-filea have samples with unknown gender
$apt/ps-metrics \
    --posterior-file "$set_output"/step2/AxiomGT1.snp-posteriors.txt \
    --call-file "$set_output"/step2/AxiomGT1.calls.txt \
    --summary-file "$set_output"/step2/AxiomGT1.summary.txt \
    --report-file "$set_output"/step2/AxiomGT1.report.txt \
    --special-snps "$specialSNPs" \
    --metrics-file "$set_output"/SNPolisher/metrics.txt \
    --output-dir "$set_output"/SNPolisher


## 8B,8C: Classify SNPs using QC metrics & Create a recommended SNP list
# ps-classification function is used to sort each SNPs/probesets into eight classes. If a genotype frequency file was not provided, the UnexpectedGenoFreq category will not appear. The category "Hemizygous" appears only when hemizygous SNPs are in the metrics file but the "special SNPs" data was not provided to ps-metrics.

# Priority order: ps-classification also identifies the best probeset per SNP using the Priority order
# Default: a) axiom_genotyping_solution_analysis_guide.pdf: PolyHighResolution, NoMinorHom, MonoHighResolution, OTV, UnexpectedGenotypeFreq, CallRateBelowThreshold, Other, and OtherMA. b) The downloaded html doc: PolyHighResolution,NoMinorHom,OTV,MonoHighResolution,CallRateBelowThreshold. c) Ps.performance.txt (output of ps-classification): PolyHighResolution, NoMinorHom, MonoHighResolution, OTV, UnexpectedGenotypeFreq, CallRateBelowThreshold, Other, OtherMA
# The priority-order argument allows the user to change the order of categories when determining which probesets are selected as the best probeset for a SNP.

# Recommended checklist: SNPs that are categorized into recommended categories are of high quality and can be used in downstream analysis. Choice of recommended categories catagories depend on the genome type (See table 5). OTV probsets may re-classified into a recommended category after the otv-caller function is used for genotyping. (re-run ps-classification on the new genotype calls; See “Adjust genotype calls for OTV SNPs” on page 40)
# Default: a) axiom_genotyping_solution_analysis_guide.pdf: for Diploid the default is PolyHighResolution, MonoHighResolution, NoMinorHom; b) The downloaded html doc: for Diploid the default is PolyHighResolution, MonoHighResolution, NoMinorHom. c) Ps.performance.txt (output of ps-classification): PolyHighResolution, NoMinorHom, MonoHighResolution, Hemizygous
# The recommended argument recives the list of categories whose SNPs will be output as recommended. If output-recommended is set to TRUE,

# The ps-classification function outputs the Ps.performance.txt file, which contains the probeset_id’s, QC metrics, hemizygous status, and an indicator if this probeset is the best for the SNP (BestProbeset), and which category the probeset belongs to (ConversionType) for each probeset. If all SNPs have one probeset, then every probeset is the best probeset by default. Any probeset that is in a recommended category and is the best probeset for a SNP will also be the BestandRecommended probeset for that SNP (see Step 8C for more details on recommended SNPs). If a ps2snp file has been provided, the snp_id column is included. Note that some versions of the performance file may have “affy_snp_id” as a column name instead of “snpid”. Column names and examples are shown in the following table (Table 4).

# for Array_A: the ps2snp map has 643641 probeset for 444805 SNP
# for Array_B: the ps2snp map has 625277 probeset for 625277 SNP
$apt/ps-classification\
    --species-type diploid \
    --metrics-file "$set_output"/SNPolisher/metrics.txt \
    --output-dir "$set_output"/SNPolisher \
    --ps2snp-file "$ps2snp"

## Array_A ## 401250 output/setA/SNPolisher/Recommended.ps
## Array_B ## 547233 output/setB/SNPolisher/Recommended.ps

## Export
## create map files for CEL file names and sample IDs
echo "Sample Filename|Alternate Sample Name" | tr '|' '\t' > "$set_output"/info_names.txt
tail -n+2 "$pedInfo" | awk 'BEGIN{FS=OFS="\t"}{print $1,$3}' >> "$set_output"/info_names.txt

$apt/apt-format-result \
      --calls-file "$set_output"/step2/AxiomGT1.calls.txt \
      --snp-list-file "$set_output"/SNPolisher/Recommended.ps \
      --annotation-file "$annot_db" \
      --snp-identifier-column Affy_SNP_ID \
      --export-chr-shortname true \
      --export-vcf-file "$set_output"/export/AxiomGT1.vcf \
      --log-file "$set_output"/export/AxiomGT1_to_vcf.log \
      --export-alternate-sample-names true \
      --sample-attributes-file "$set_output"/info_names.txt
# array A: #    WARNING  1 1824 | VCF file will be generated without the required REF and ALT columns containing necessary information. (calls file)       Count: 11328
# array B: #    WARNING  1 2423 | VCF file will be generated without the required REF and ALT columns containing necessary information. (calls file)       Count: 11188


wc -l "$set_output"/SNPolisher/Recommended.ps ## 401250 // 547233
grep -v "^#" $set_output/export/AxiomGT1.vcf | wc -l  ## 401249 // 547232
## explore available Ref alleles
grep -v "^#" $set_output/export/AxiomGT1.vcf | cut -f4 | sort | uniq -c ## REF=N 11328 // 11188
cat $set_output/export/AxiomGT1.vcf | awk 'BEGIN{FS=OFS="\t"}/^#/{print}!/^#/{if($4=="N")print}' | cut -f1-10 > $set_output/export/noRef.vcf

############## exploration
cd $work_dir/lib2/setA
paste <(head -n1 Axiom_K9HDSNPA_Annotation.r2_3.tab | tr '\t' '\n' | awk 'BEGIN{FS=OFS="\t"}{print NR,$0}') <(head -n2 Axiom_K9HDSNPA_Annotation.r2_3.tab | tail -n1 | tr '\t' '\n') ## columns 14-18: Flank, Allele_A, Allele_B, Ref_Allele, Alt_Allele && column 32: genome(i.e. assembly)
tail -n+2 Axiom_K9HDSNPA_Annotation.r2_3.tab | cut -f32 | sort | uniq -c > genome.dist ## 1597 with no defined genome
tail -n+2 Axiom_K9HDSNPA_Annotation.r2_3.tab | awk 'BEGIN{FS=OFS="\t"}{if($32=="")print $3}' | sort | uniq -c ## 1159 (unknown) and 438 (chrY) Note: the chr numbers were annotated based on Axiom_K9HDSNPA_Annotation.r2_3.chr.tab extracted from the sqlite database
tail -n+2 Axiom_K9HDSNPA_Annotation.r2_3.tab | cut -f3 | sort | uniq -c > chr.dist

tail -n+2 Axiom_K9HDSNPA_Annotation.r2_3.tab | cut -f17 | sort | uniq -c > REF.dist ## 19886 no ref
tail -n+2 Axiom_K9HDSNPA_Annotation.r2_3.tab | awk 'BEGIN{FS=OFS="\t"}{if($17=="")print $3}' | sort | uniq -c > noRef_chr.dist

tail -n+2 Axiom_K9HDSNPA_Annotation.r2_3.tab | awk 'BEGIN{FS=OFS="\t"}{if($17=="")print $15,$16}' | sort | uniq -c > noRef_Alleles.dist
tail -n+2 Axiom_K9HDSNPA_Annotation.r2_3.tab | awk 'BEGIN{FS=OFS="\t"}{if($17=="" && $15=="-")print}' | head
##############
## Update the output VCF to remove all useless/bad markers
cd "$set_output"/export/

####### Extra: start
## QC for SNPs without ref/alt
## extract all SNPs without ref/alt & select those among recommended ## column 17 = "Ref Allele"
#grep -v "^#" $annot_csv | sed 's/^"//;s/\",\"/|/g;s/"$//;' | awk -F"|" '{if($19=="---")print}' > noREF.csv
grep -Fwf "$set_output"/SNPolisher/Recommended.ps <(grep -v "^#" $annot_csv | awk -F"\t" '{if($17=="")print}') > noREF_rec.tab #11328 // 11188
rclone -v --copy-links copy noREF_rec.tab remote_UCDavis_GoogleDr:MAF/genotyping_array_v2/$(basename "$set_output") ## Thermo
cat noREF_rec.tab | cut -f3 | sort | uniq -c > noREF_rec.chr
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2;next;}{print $0,a[$2]}' "$chr_map" <(cat noREF_rec.chr | sed 's/^ *//' | tr ' ' '\t') > noREF_rec.chr.map
## Split into variants with known Chr/Pos ($3 and $4) and variants without
cat noREF_rec.tab | awk -F"\t" '{if($3!="" && $4!="")print}' > noREF_rec_known.tab ## 10809 // 11188
cat noREF_rec.tab | awk -F"\t" '{if($3=="" || $4=="")print}' > noREF_rec_unknown.tab ## 519 // 0
####### Extra: end

## Update the VCF chromosomes to match canFam3 && remove UNKNOWNCHR
# change 1,2,3,.....X,Y to chr1,chr2,...... chrY
# change MT to chrM
# change Un to chrUn_xxxxxxx
# remove UNKNOWNCHR
grep -v "^#" AxiomGT1.vcf | cut -f1 | uniq -c > chr.lst
cat chr.lst | awk 'BEGIN{OFS="\t"}{print $2,"chr"$2}' | sed 's/chrMT/chrM/' | grep -v "UNKNOWNCHR" > allChr.map
cat allChr.map | cut -f2 | while read f;do echo "##contig=<ID=$f>";done > contig_ids
cat AxiomGT1.vcf | awk -v file=contig_ids 'BEGIN{contig=1;}/^##[^contig]/{print}/^##contig/{if(contig){while((getline<file) > 0)print;contig=0;}}/^#CH/{print}!/^#/{exit}' > AxiomGT1v2.vcf
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2;next}!/#/{if(a[$1]){$1=a[$1];print}else if(a[$3]){$1=a[$3];print}}' allChr.map AxiomGT1.vcf >> AxiomGT1v2.vcf ## grep -v "^#" AxiomGT1v2.vcf | wc -l ## 400730 // 547232
## identify variants without known chromosome
grep -v "^#" AxiomGT1.vcf | awk 'BEGIN{FS=OFS="\t"}{if($1=="UNKNOWNCHR")print $3}' > noChr.lst ## 519 // 0


#conda create -n equineSNP
#conda install -c bioconda plink
#conda install -c bioconda plink2
#conda install -c bioconda bcftools
#conda install -c bioconda eigensoft
#conda install -c bioconda vcftools
conda activate equineSNP

## check the VCF integrity
bcftools sort -o AxiomGT1v2_sorted.vcf AxiomGT1v2.vcf
bcftools norm -c ws -f $canFam3_ref_unwrap AxiomGT1v2_sorted.vcf 1> AxiomGT1v2.check.vcf 2> AxiomGT1v2.check.log
# Array A
#Lines   total/split/realigned/skipped:  400730/0/10/0
#REF/ALT total/modified/added:   400730/0/10818

# Array B
#Lines   total/split/realigned/skipped:  547232/0/0/0
#REF/ALT total/modified/added:   547232/0/11188


## Explore the changed markers
paste <(grep -v "^##" AxiomGT1v2_sorted.vcf | cut -f1-5) <(grep -v "^##" AxiomGT1v2.check.vcf | cut -f1-5)  > sorted.vs.check.tab
cat sorted.vs.check.tab | awk -F"\t" '{a=$1FS$2FS$3FS$4FS$5; b=$6FS$7FS$8FS$9FS$10; if(a != b)print}' > sorted.vs.check.dif.tab ## 10818 // 11188
comm -12 <(cut -f2 noREF_rec_known.tab | sort) <(cut -f3 sorted.vs.check.dif.tab | sort) | wc -l #10805 (out of 10809 no ref) // 11188
comm -23 <(cut -f2 noREF_rec_known.tab | sort) <(cut -f3 sorted.vs.check.dif.tab | sort) | grep -Fwf - sorted.vs.check.tab ## We have 4 Y markers that keep with N allele after bcftools norm
cat sorted.vs.check.dif.tab | awk -F"\t" '{if($4 != $9)print}' > sorted.vs.check.dif.changeAlelle.tab ## 10819 // 11188
cat sorted.vs.check.dif.tab | awk -F"\t" '{if($2 != $7)print}' > sorted.vs.check.dif.changePos.tab ## 10 // 0
comm -12 <(cat sorted.vs.check.dif.changeAlelle.tab | sort) <(cat sorted.vs.check.dif.changePos.tab | sort) ## 10 (all the changePos are changeAlelle)
comm -12 <(cut -f2 noREF_rec_known.tab | sort) <(cut -f3 sorted.vs.check.dif.changePos.tab | sort) ## none of the changePos were with noREF
comm -12 <(cut -f2 noREF_rec_known.tab | sort) <(cut -f3 sorted.vs.check.dif.changeAlelle.tab | sort) | wc -l ## 10805
comm -13 <(cut -f2 noREF_rec_known.tab | sort) <(cut -f3 sorted.vs.check.dif.changeAlelle.tab | sort) | grep -v -Fwf <(cut -f3 sorted.vs.check.dif.changePos.tab ) | grep -Fwf - sorted.vs.check.dif.changeAlelle.tab ## 4 markers changed their wrong alelles
## identify remaining variants without known ref/alt alleles
grep -v "^#" AxiomGT1v2.check.vcf | awk 'BEGIN{FS=OFS="\t"}{if($4=="N")print $3}' > noRef.lst ## 4  // 0
## identify remaining variants with annotation errors
cat sorted.vs.check.dif.changeAlelle.tab | awk 'BEGIN{FS=OFS="\t"}{if($4!="N")print $3}' > changeAlelle.lst ## 14 // 0
## Conclusion:
## Among the recommended 401249 markers in array A, there are 11328 with unannotated reference allele with 519 of them having unknown chromosome.
## We exlcuded these 519 markers and tested the integrity of the remaning 400730 markers that include the remaining 10809 markers with unannotated reference allele. We found:
## a) The 10809 markers with unannotated reference allele were divided into:
##    4 Y markers that have no reference alalle in the genome and thus can't be updated
##    10805 markers that we were able to update the reference allele
## b) Additional problematic 14 merkers
##    4 markers that needed to change the reference allele
##    10 markers that needed to change their postion as well as their reference allele

## Update the alterantive allele based on array annotation
#1. update the letter format of the REF/ALT alleles in the VCF from 'acgt' to 'ACGT'
grep -v "^#" AxiomGT1v2.check.vcf | cut -f1-5 > subVCF_temp
cut -f4,5 subVCF_temp | tr 'acgt' 'ACGT' > subVCF_temp2
paste subVCF_temp subVCF_temp2 | awk 'BEGIN{FS=OFS="\t"}{if($4!=$6 || $5!=$7)print}' > diff_temp
paste subVCF_temp subVCF_temp2 | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$6,$7}' > AxiomGT1v2.check.tab
rm subVCF_temp subVCF_temp2
#2. merge REF/ALT of vcf with Flank, Allele_A, Allele_B from annotation (& adjust the Flank alleles with Allele_A/B)
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$2]=$15 FS $16 FS $14;next}{if(a[$3])print $0,a[$3]}' $annot_csv AxiomGT1v2.check.tab | tr '[/]' '\t' > AxiomGT1v2.check.ann.tab ## REF/ALT = 4/5, Allele_A/B = 6/7, Flank_A/B=9/10
cat AxiomGT1v2.check.ann.tab | awk 'BEGIN{FS=OFS="\t"}{if($6==$10 && $7==$9){a=$9;$9=$10;$10=a;}print}' > AxiomGT1v2.check.ann2.tab ## switched Flank_A/B (9/10) to match Allele_A/B (6/7)
cat AxiomGT1v2.check.ann2.tab | awk 'BEGIN{FS=OFS="\t"}{if($6!=$9 || $7!=$10)print}' | wc -l ##0 // 0
#3. update the deletion allele_A/B
cat AxiomGT1v2.check.ann2.tab | awk 'BEGIN{FS=OFS="\t"}{if($6=="-")print}' | wc -l ## 654 // 0
cat AxiomGT1v2.check.ann2.tab | awk 'BEGIN{FS=OFS="\t"}{if($7=="-")print}' | wc -l ## 0 // 0
cat AxiomGT1v2.check.ann2.tab | awk 'BEGIN{FS=OFS="\t"}{if($6=="-"){$6=substr($8,length($8),1);$7=$6$7}print}' > AxiomGT1v2.check.ann3.tab ## update Allele_A/B (6/7)
#4. update the alterantive allele based on array annotation
cat AxiomGT1v2.check.ann3.tab | awk 'BEGIN{FS=OFS="\t"}{if($5=="."){if($4==$6)$5=$7;else if($4==$7)$5=$6;}print}' > AxiomGT1v2.check.ann4.tab
cat AxiomGT1v2.check.ann4.tab | awk 'BEGIN{FS=OFS="\t"}{if(($4==$6 && $5==$7) || ($4==$7 && $5==$6))print}' > AxiomGT1v2.check.ann_updated.tab
#5. identify the variants with bad annotations
cat AxiomGT1v2.check.ann4.tab | awk 'BEGIN{FS=OFS="\t"}{if(!(($4==$6 && $5==$7) || ($4==$7 && $5==$6)))print $3}' > annError.lst ## 138 // 0
cat annError.lst | grep -Fwf - AxiomGT1v2.check.ann4.tab | awk 'BEGIN{FS=OFS="\t"}{if($5==".")print}' > annError_noRef.lst ## 124
cat annError.lst | grep -Fwf - AxiomGT1v2.check.ann4.tab | awk 'BEGIN{FS=OFS="\t"}{if($5!=".")print}' > annError_changeAllele.lst ## 14
# compare these errors to those identified earlier from exploring the output of bcftools norm
comm -12 <(cat noRef.lst changeAlelle.lst | sort) <(cat annError.lst | sort) | wc -l ## 18 = 4+14 // 0

## Exclude bad variants
cat annError.lst | grep -vFwf - AxiomGT1v2.check.vcf > AxiomGT1v2.hiQ.vcf
grep -v "^#" AxiomGT1v2.hiQ.vcf | wc -l ## 400592 // 547232

####### Extra: start
## Generate lists of annotation errors for Thermo
cat annError.lst | grep -Fwf - AxiomGT1v2.vcf | cut -f3,8 > annErrorProbes.lst
echo "chr pos probe-set ref alt CF3-ref marker" | tr ' ' '\t' > annErrorProbes.map
awk 'BEGIN{FS=OFS="\t"}FNR==NR{if($4!="N" && $4!=$10)a[$3]=$0;next}{if(a[$1])print a[$1],$2}' sorted.vs.check.dif.tab annErrorProbes.lst | sed 's/probeset_id=//' | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$11,$4,$5,$9,$3}' >> annErrorProbes.map
rclone -v --copy-links copy annErrorProbes.map remote_UCDavis_GoogleDr:MAF/genotyping_array_v2/$(basename "$set_output")

#cat check.dif.map | awk 'BEGIN{FS=OFS="\t"}{if($4==$7 && $5==$6)print $3}' > strandSwab.lst
#cat strandSwab.lst | grep -Fwf - AxiomGT1v2.vcf | cut -f3,8 > strandSwabProbes.lst
#echo "chr pos probe-set ref alt CF3-ref marker" | tr ' ' '\t' > strandSwabProbes.map
#awk 'BEGIN{FS=OFS="\t"}FNR==NR{if($4==$7 && $5==$6)a[$3]=$0;next}{if(a[$1])print a[$1],$2}' check.dif.map strandSwabProbes.lst | sed 's/probeset_id=//' | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$8,$4,$5,$6,$3}' >> strandSwabProbes.map
#rclone -v --copy-links copy strandSwabProbes.map remote_UCDavis_GoogleDr:MAF/genotyping_array_v2/$(basename "$set_output")
####### Extra: end


## There are different SNPs with the same position on the same array
## Keep one of them which has higher genotyping rate
grep -v "^#" AxiomGT1v2.hiQ.vcf | awk '{print $1"_"$2}' | sort > pos.lst
cat pos.lst | uniq -c | awk '{if($1>1)print $2}' > pos_dup.lst # 2907 // 0
cat pos.lst | uniq -c | awk '{if($1>2)print $2}' ## Not more than 2 IDs per position
awk -F"\t" 'FNR==NR{a[$1]=1;next}{if(a[$1"_"$2])print}' pos_dup.lst AxiomGT1v2.hiQ.vcf | cut -f1,2 | paste - - | awk '{if($1!=$3 || $2!=$4)print}' ## nothing (safe to paste)
echo "CHROM POS ID1 NS1 MAF1 AC_Het1 AC_Hom1 ExcHet1 ID2 NS2 MAF2 AC_Het2 AC_Hom2 ExcHet2" | tr ' '  '\t' > pos_dup.ann
awk -F"\t" 'FNR==NR{a[$1]=1;next}/^#/{print}{if(a[$1"_"$2])print}' pos_dup.lst AxiomGT1v2.hiQ.vcf | bcftools +fill-tags | bcftools query -f'%CHROM\t%POS\t%ID\t%NS\t%MAF\t%AC_Het\t%AC_Hom\t%ExcHet\n' | paste - - | cut -f1-8,11-16 >> pos_dup.ann #pos_dup.vcf
tail -n+2 pos_dup.ann | awk -F"\t" '{if($4<$10)print $3;else print $9}' > pos_dup_ToExclude.lst ## 2907 // 0
cat pos_dup_ToExclude.lst | grep -vFwf - AxiomGT1v2.hiQ.vcf > AxiomGT1v2.hiQ2.vcf
grep -v "^#" AxiomGT1v2.hiQ2.vcf | awk '{print $1"_"$2}' | sort > pos2.lst ## 397685 // 547232

####### Extra: start
tail -n+2 pos_dup.ann | cut -f3 | grep -Fwf - AxiomGT1v2.vcf | cut -f1-3,8 > ID1.lst
tail -n+2 pos_dup.ann | cut -f9 | grep -Fwf - AxiomGT1v2.vcf | cut -f1-3,8 > ID2.lst
echo "chr pos probe-set1 affy_snp_id1 probe-set2 affy_snp_id2" | tr ' ' '\t' > dupMarkersWithinArrA.map
paste ID1.lst ID2.lst | awk 'BEGIN{FS=OFS="\t"}{if($1==$5 && $2==$6)print $1,$2,$4,$3,$8,$7}' | sed 's/probeset_id=//g'  >> dupMarkersWithinArrA.map
rclone -v --copy-links copy dupMarkersWithinArrA.map remote_UCDavis_GoogleDr:MAF/genotyping_array_v2
####### Extra: end


## Export to Plink to keep the genotypes of all markers then convert the Plink files to VCF using the list in the curated VCF files
$apt/apt-format-result \
      --calls-file "$set_output"/step2/AxiomGT1.calls.txt \
      --snp-list-file "$set_output"/SNPolisher/Recommended.ps \
      --annotation-file "$annot_db" \
      --snp-identifier-column Affy_SNP_ID \
      --export-chr-shortname true \
      --export-plink-file "$set_output"/export_plink/AxiomGT1 \
      --log-file "$set_output"/export_plink/AxiomGT1_to_plink.log \
      --export-alternate-sample-names true \
      --sample-attributes-file "$set_output"/info_names.txt \
      --pedigree-file "$pedInfo"


## Update gender based on Array A predictions
# change to male(1): S000055,S016242,S020302,S020913,S028037,S040135,S040228,S040259
# change to female(2): S000054,S020301,S020391,S028073,S040085,S040195,S040234
echo "S000055,S016242,S020302,S020913,S028037,S040135,S040228,S040259" | tr ',' '\n' | awk 'BEGIN{OFS="\t"}{print $0,"1"}' > $work_dir/new_gender
echo "S000054,S020301,S020391,S028073,S040085,S040195,S040234" | tr ',' '\n' | awk 'BEGIN{OFS="\t"}{print $0,"2"}' >> $work_dir/new_gender
mv "$set_output"/export_plink/AxiomGT1.ped "$set_output"/export_plink/AxiomGT1.ped_bk
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2;next}{if(a[$2])$5=a[$2];print}' $work_dir/new_gender "$set_output"/export_plink/AxiomGT1.ped_bk > "$set_output"/export_plink/AxiomGT1.ped
mv "$set_output"/export_plink/AxiomGT1.map "$set_output"/export_plink/AxiomGT1.map_bk
awk 'BEGIN{FS=OFS="\t"}/^#/{print;next}{print "chr"$0}' "$set_output"/export_plink/AxiomGT1.map_bk > "$set_output"/export_plink/AxiomGT1.map

# select clean markers from the curated VCF files and convert Plink files to binary Plink files
grep -v "^#" "$set_output"/export/AxiomGT1v2.hiQ2.vcf | cut -f3 > "$set_output"/export/keep_markers.lst ## 397685 // 547232
plink --file "$set_output"/export_plink/AxiomGT1 --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
      --extract "$set_output"/export/keep_markers.lst \
      --make-bed --output-chr 'chrM' --out "$set_output"/export_plink/AxiomGT1.bin

# update the alleles based on AxiomGT1v2.check.ann_updated.tab
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$2]=$0;next}{if(a[$3])print a[$3],$0}' "$set_output"/export_plink/AxiomGT1.bin.bim "$set_output"/export/AxiomGT1v2.check.ann_updated.tab > "$set_output"/export_plink/AxiomGT1.bin.bim_ann
mv "$set_output"/export_plink/AxiomGT1.bin.bim "$set_output"/export_plink/AxiomGT1.bin.bim_bk
cat "$set_output"/export_plink/AxiomGT1.bin.bim_ann | awk 'BEGIN{FS=OFS="\t"}{if($6==$16)print $1,$2,$3,$8,$12,$13;else if($6==$15)print $1,$2,$3,$8,$13,$12;}' > "$set_output"/export_plink/AxiomGT1.bin.bim

# convert Plink to VCF using the list in the curated VCF files
plink --bfile "$set_output"/export_plink/AxiomGT1.bin --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
      --a2-allele "$set_output"/export/AxiomGT1v2.check.ann_updated.tab 4 3 \
      --recode vcf-iid --output-chr 'chrM' --out "$set_output"/export_plink/AxiomGT1
      
grep -v "^#" "$set_output"/export_plink/AxiomGT1.vcf | awk '{print $1"_"$2}' | sort > "$set_output"/export_plink/pos3.lst ## 397685 // 547232

#########
## Assess duplicates and genotype concordance
cd $work_dir
## shared samples
grep "^#CHROM" output/setA/export_plink/AxiomGT1.vcf | tr '\t' '\n' | tail -n+10 > A.samples ## 3355
grep "^#CHROM" output/setB/export_plink/AxiomGT1.vcf | tr '\t' '\n' | tail -n+10 > B.samples ## 3339
comm -12 <(sort A.samples) <(sort B.samples) > shared.samples ## 3338
comm -13 <(sort A.samples) <(sort B.samples) > Bonly.samples ## 1
comm -23 <(sort A.samples) <(sort B.samples) > Aonly.samples ## 17

## shared SNPs
comm -12 output/set[AB]/export_plink/pos3.lst > dup_pos # 30933
awk -F"\t" 'FNR==NR{a[$1]=1;next}/^#/{print}{if(a[$1"_"$2])print}' dup_pos output/setA/export_plink/AxiomGT1.vcf > dup_A.vcf
awk -F"\t" 'FNR==NR{a[$1]=1;next}/^#/{print}{if(a[$1"_"$2])print}' dup_pos output/setB/export_plink/AxiomGT1.vcf > dup_B.vcf
bgzip dup_A.vcf && tabix -p vcf dup_A.vcf.gz
bgzip dup_B.vcf && tabix -p vcf dup_B.vcf.gz
bcftools stats --samples-file shared.samples dup_A.vcf.gz dup_B.vcf.gz > dup_stats
records=$(cat <(echo "add line") shared.samples | wc -l)
grep -A $records "Genotype concordance by sample (SNPs)" dup_stats > dup_Geno_concord
tail -n+3 dup_Geno_concord | cut -f11 | awk '{a+=1;b+=$1}END{print b/a}' ## 0.970413

echo "CHROM POS ID1 NS1 MAF1 AC_Het1 AC_Hom1 ExcHet1" | tr ' '  '\t' > dup_A.info
zcat dup_A.vcf.gz | bcftools +fill-tags | bcftools query -f'%CHROM\t%POS\t%ID\t%NS\t%MAF\t%AC_Het\t%AC_Hom\t%ExcHet\n' >> dup_A.info
echo "CHROM POS ID2 NS2 MAF2 AC_Het2 AC_Hom2 ExcHet2" | tr ' '  '\t' > dup_B.info
zcat dup_B.vcf.gz | bcftools +fill-tags | bcftools query -f'%CHROM\t%POS\t%ID\t%NS\t%MAF\t%AC_Het\t%AC_Hom\t%ExcHet\n' >> dup_B.info

paste dup_A.info dup_B.info | awk -F"\t" '{if($2!=$10)print}' ## safe to paste
paste dup_A.info dup_B.info | cut -f1-8,11-16 > dup_AvsB.info
#head dup_AvsB.info
#grep Affx-205694835 output/setA/export/AxiomGT1.vcf | tr '\t' '\n' > 1.temp
#grep Affx-205694835 output/setB/export/AxiomGT1.vcf | tr '\t' '\n' > 2.temp
#grep "^#CHROM" output/setA/export/AxiomGT1.vcf | tr '\t' '\n' > 1n.temp
#grep "^#CHROM" output/setB/export/AxiomGT1.vcf | tr '\t' '\n' > 2n.temp
#cat 2n.temp | grep -Fwf - <(paste 1n.temp 1.temp) > 1nc.temp
#cat 1n.temp | grep -Fwf - <(paste 2n.temp 2.temp) > 2nc.temp
#paste 1nc.temp 2nc.temp  | awk '{if($1!=$3)print}'
#paste 1nc.temp 2nc.temp  | awk '{if($2!=$4)print}' | grep -v "\./\."

## Thermo
cat dup_AvsB.info | awk -F"\t" '{if($3!=$9)print}' | wc -l ## 1464
cat dup_AvsB.info | cut -f1-3,9 > dup_AvsB.map
rclone -v --copy-links copy dup_AvsB.map remote_UCDavis_GoogleDr:MAF/genotyping_array_v2
rclone -v --copy-links copy dup_A.vcf.gz remote_UCDavis_GoogleDr:MAF/genotyping_array_v2
rclone -v --copy-links copy dup_B.vcf.gz remote_UCDavis_GoogleDr:MAF/genotyping_array_v2

## define the duplicates SNPs to be removed
tail -n+2 dup_B.info | awk -F"\t" '{print $3}' > pos_dupAcrossArr_ToExclude.lst

## A) VCF files (Obsolete; if a VCF file version is needed, convert the new Plink files to VCF)
## remove duplicates SNPs from array B VCF file
cat pos_dupAcrossArr_ToExclude.lst | grep -vFwf - output/setB/export_plink/AxiomGT1.vcf > output/setB/export_plink/AxiomGT1_noDup.vcf
grep -v "^#" output/setB/export_plink/AxiomGT1_noDup.vcf | wc -l ## 516299

## prep final array VCF files to final merge
bgzip output/setA/export_plink/AxiomGT1.vcf && tabix -p vcf output/setA/export_plink/AxiomGT1.vcf.gz
bgzip output/setB/export_plink/AxiomGT1_noDup.vcf && tabix -p vcf output/setB/export_plink/AxiomGT1_noDup.vcf.gz

bcftools view --samples-file shared.samples output/setA/export_plink/AxiomGT1.vcf.gz | bgzip -c > output/setA/export_plink/AxiomGT1v2.shared.vcf.gz
tabix -p vcf output/setA/export_plink/AxiomGT1v2.shared.vcf.gz
bcftools view --samples-file shared.samples output/setB/export_plink/AxiomGT1_noDup.vcf.gz | bgzip -c > output/setB/export_plink/AxiomGT1v2.shared.vcf.gz
tabix -p vcf output/setB/export_plink/AxiomGT1v2.shared.vcf.gz

#bcftools view --samples-file Aonly.samples output/setA/export/AxiomGT1v2.hiQ2.vcf.gz | bgzip -c > output/setA/export/AxiomGT1v2.Aonly.vcf.gz
#tabix -p vcf output/setA/export/AxiomGT1v2.Aonly.vcf.gz
#bcftools view --samples-file Bonly.samples output/setB/export/AxiomGT1v2.hiQ3.vcf.gz | bgzip -c > output/setB/export/AxiomGT1v2.Bonly.vcf.gz
#tabix -p vcf output/setB/export/AxiomGT1v2.Bonly.vcf.gz

## merge VCF files
bcftools concat --allow-overlaps output/setA/export_plink/AxiomGT1v2.shared.vcf.gz output/setB/export_plink/AxiomGT1v2.shared.vcf.gz | bgzip -c > AxiomGT1v2.merge.vcf.gz
tabix -p vcf AxiomGT1v2.merge.vcf.gz

bcftools norm -c ws -f $canFam3_ref_unwrap AxiomGT1v2.merge.vcf.gz 1> AxiomGT1v2.merge.check.vcf 2> AxiomGT1v2.merge.check.log
rm AxiomGT1v2.merge.check.vcf

#Lines   total/split/realigned/skipped:  913984/0/0/0
#REF/ALT total/modified/added:   913984/0/0

#bcftools merge -m none AxiomGT1v2.merge.vcf.gz output/setA/export/AxiomGT1v2.Aonly.vcf.gz output/setB/export/AxiomGT1v2.Bonly.vcf.gz | bcftools +fill-tags | bgzip -c > AxiomGT1v2.merge2.vcf.gz
#tabix -p vcf AxiomGT1v2.merge2.vcf.gz


## B) Prep and merge the Plink files (Obsolete; See the new v2)
## remove duplicates SNPs from array B
plink --bfile output/setB/export_plink/AxiomGT1.bin --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
      --exclude pos_dupAcrossArr_ToExclude.lst \
      --make-bed --output-chr 'chrM' --out output/setB/export_plink/AxiomGT1.bin_noDup ## 516299 variants remaining

## prep final array files to final merge
cat shared.samples | grep -Fwf - output/setA/export_plink/AxiomGT1.bin.fam | cut -d" " -f1,2 > shared.samples.forPlink
plink --bfile output/setA/export_plink/AxiomGT1.bin --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
      --keep shared.samples.forPlink \
      --make-bed --output-chr 'chrM' --out output/setA/export_plink/AxiomGT1.shared ## 3338 samples remaining

plink --bfile output/setB/export_plink/AxiomGT1.bin_noDup --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
      --keep shared.samples.forPlink \
      --make-bed --output-chr 'chrM' --out output/setB/export_plink/AxiomGT1.shared ## 3338 samples remaining

## merge
plink --bfile output/setA/export_plink/AxiomGT1.shared --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
      --bmerge output/setB/export_plink/AxiomGT1.shared \
      --make-bed --output-chr 'chrM' --out AxiomGT1v2.merge
#3338 samples loaded from output/setA/export_plink/AxiomGT1.shared.fam.
#3338 samples to be merged from output/setB/export_plink/AxiomGT1.shared.fam.
#Of these, 0 are new, while 3338 are present in the base dataset.
#397685 markers loaded from output/setA/export_plink/AxiomGT1.shared.bim.
#516299 markers to be merged from output/setB/export_plink/AxiomGT1.shared.bim.
#Of these, 516299 are new, while 0 are present in the base dataset.
#Performing single-pass merge (3338 samples, 913984 variants).
#Total genotyping rate is 0.997558.



## Bv2) Prep and merge the Plink files
## --merge-mode 6 = (no merge) Report all mismatching calls.
#plink --bfile output/setB/export_plink/AxiomGT1.bin --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
      --bmerge output/setA/export_plink/AxiomGT1.bin --merge-mode 6 \
      --output-chr 'chrM' --out AxiomGT1.mismatching6
#wc -l AxiomGT1.mismatching6.diff #771445 AxiomGT1.mismatching6.diff

## --merge-mode 7 = (no merge) Report mismatching nonmissing calls.
plink --bfile output/setB/export_plink/AxiomGT1.bin --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
      --bmerge output/setA/export_plink/AxiomGT1.bin --merge-mode 7 \
      --output-chr 'chrM' --out AxiomGT1.mismatching7
      
#3339 samples loaded from output/setB/export_plink/AxiomGT1.bin.fam.
#3355 samples to be merged from output/setA/export_plink/AxiomGT1.bin.fam.
#Of these, 17 are new, while 3338 are present in the base dataset.
#547232 markers loaded from output/setB/export_plink/AxiomGT1.bin.bim.
#397685 markers to be merged from output/setA/export_plink/AxiomGT1.bin.bim.
#Of these, 368215 are new, while 29470 are present in the base dataset.
#Warning: Variants 'Affx-206088448' and 'Affx-205939940' have the same position.
#Warning: Variants 'Affx-205859745' and 'Affx-205344060' have the same position.
#Warning: Variants 'Affx-206706550' and 'Affx-205359377' have the same position.
#1460 more same-position warnings: see log file.
#Performing 1-pass diff (mode 7), writing results to AxiomGT1.mismatching7.diff
#98370860 overlapping calls, 97897061 nonmissing in both filesets.
#97593330 concordant, for a concordance rate of 0.996897.

wc -l AxiomGT1.mismatching7.diff #303732 AxiomGT1.mismatching7.diff
tail -n+2 AxiomGT1.mismatching7.diff | awk '{print $1}' | sort | uniq | wc -l ## 21849
tail -n+2 AxiomGT1.mismatching7.diff | awk '{print $2,$3}' | sort | uniq -c | sort -k1,1nr > AxiomGT1.mismatching7.diff.samples ## two samples show exceptional higher rate of mismatching (S007258 and S019740 have 11676 and 11440 mismatches respectively. The latter sample is the only sample that was predicted to be female on Array A and male on array B. This likely indicate that these samples had something wrong. Swapping them on array B did not fix the issue)

## test mismatching after excluding the likely swapped samples
echo "GRLS S007258|GRLS S019740" | tr '|' '\n' > swap_samples.lst
plink --bfile output/setA/export_plink/AxiomGT1.bin --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
      --remove swap_samples.lst \
      --make-bed --output-chr 'chrM' --out output/setA/export_plink/AxiomGT1.bin_noSwap
plink --bfile output/setB/export_plink/AxiomGT1.bin --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
      --remove swap_samples.lst \
      --make-bed --output-chr 'chrM' --out output/setB/export_plink/AxiomGT1.bin_noSwap

plink --bfile output/setB/export_plink/AxiomGT1.bin_noSwap --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
      --bmerge output/setA/export_plink/AxiomGT1.bin_noSwap --merge-mode 7 \
      --output-chr 'chrM' --out AxiomGT1_noSwap.mismatching7  ## concordance rate of 0.997132.


# investigate the same position markers (1463)
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1 FS $4]=$2;b1[$1 FS $4]=$2 FS $5 FS $6;b2[$1 FS $4]=$2 FS $6 FS $5;next}{if(a[$1 FS $4] && b1[$1 FS $4]!=$2 FS $5 FS $6 && b2[$1 FS $4]!=$2 FS $5 FS $6)print b1[$1 FS $4],$2,$5,$6}' output/setA/export_plink/AxiomGT1.bin.bim output/setB/export_plink/AxiomGT1.bin.bim > same_pos.warnings
awk 'BEGIN{FS=OFS="\t"}{if($2 FS $3 != $5 FS $6 && $2 FS $3 != $6 FS $5)print}' same_pos.warnings > same_pos.badAlleles
awk 'BEGIN{FS=OFS="\t"}{if($2 FS $3 == $5 FS $6 || $2 FS $3 == $6 FS $5)print}' same_pos.warnings > same_pos.sameAlleles

awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$2]=$14;next}{print $1,a[$1]}' $work_dir/lib2/setA/Axiom_K9HDSNPA_Annotation.r2_3.tab same_pos.sameAlleles  > same_pos.arrA.Flank
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$2]=$14;next}{print $4,a[$4]}' $work_dir/lib2/setB/Axiom_K9HDSNPB_Annotation.r2_3.tab same_pos.sameAlleles  > same_pos.arrB.Flank
paste same_pos.arrA same_pos.arrB | awk 'BEGIN{FS=OFS="\t"}{if($2!=$4)print}' > same_pos.badAlleles2


## remove same position markers SNPs from array B
grep "Warning: Variants .* have the same position" AxiomGT1.mismatching7.log | awk -F"'" 'BEGIN{OFS="\n";}{print $2,$4}' > same_pos.arrB.lst
plink --bfile output/setB/export_plink/AxiomGT1.bin --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
      --exclude same_pos.arrB.lst \
      --make-bed --output-chr 'chrM' --out output/setB/export_plink/AxiomGT1.bin_noSamePos ## 546519 variants remaining


## merge:
## --merge-mode 3 = Only overwrite calls which are nonmissing in the new file.
## i.e. Give priority to arrayA if genotypes are nonmissing
## AND exclude the likely swapped samples
plink --bfile output/setB/export_plink/AxiomGT1.bin_noSamePos --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
      --bmerge output/setA/export_plink/AxiomGT1.bin --merge-mode 3 \
      --remove swap_samples.lst \
      --make-bed --output-chr 'chrM' --out AxiomGT1v2.comp_merge

#Total genotyping rate in remaining samples is 0.994631.
#913984 variants and 3354 samples pass filters and QC.
###################################################################################
## Transform all files to map/ped
#plink --vcf AxiomGT1v2.merge.vcf.gz --double-id --chr-set 38 no-y no-xy --allow-extra-chr --allow-no-sex --make-bed --out "AxiomGT1v2.merge"
###################################################################################
## Generate map of IDs (My genotyping IDs, grls_ids (e.g. 094-000019), and public_ids (e.g. grlsH764T844))
## I got the files "dc_id.csv" & "golden_oldies_subject.csv" from Brenna
rclone -v --copy-links copy remote_UCDavis_GoogleDr:MAF/source_files/dc_id.csv .
cat dc_id.csv | sed -e "s/\r//g" | sed 's/\"//g' | tr ',' '\t' > dc_id.tab

rclone -v --copy-links copy remote_UCDavis_GoogleDr:MAF/source_files/golden_oldies_subject.csv .
cat golden_oldies_subject.csv | sed -e "s/\r//g" | sed 's/\"//g' | tr ',' '\t' > golden_oldies_subject.tab

echo grls_id public_id | tr ' ' '\t' > temp.map
awk 'BEGIN{FS=OFS="\t"}FNR==NR{print $3,$2;next}{print $0}' <(tail -n+2 dc_id.tab) <(tail -n+2 golden_oldies_subject.tab) | sort | uniq >> temp.map ## 3245 without header

cut -d" " -f2 AxiomGT1v2.comp_merge.fam | sed 's/^S/094-/;s/_.*$//;s/A$//;s/B$//' > temp.grls_ids ## 3354
cat temp.grls_ids | sort | uniq | wc -l  ## 3235

echo "Family_ID" "Individual_ID" "grls_id" | tr ' ' '\t' > exp_id.tab
paste <(cut -d" " -f1,2 AxiomGT1v2.comp_merge.fam | tr ' ' '\t') temp.grls_ids >> exp_id.tab
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2;next}{if(a[$3])print $0,a[$3]}' temp.map exp_id.tab > map_id.tab

# annotate the map by the metadata gender
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$2 FS $3]=$6;next}{print $0,a[$1 FS $2]}' All_pedigree/sample_ped_infoA_cohort.txt map_id.tab > map_id_sex.tab
##########################################################
## Make temporary versions of the genotyping files to share
mkdir -p toShare/{anonymous_ids,public_ids}

## Create a version with anonymous_ids
cat AxiomGT1v2.comp_merge.fam | awk '{$1=$2=NR;print}' > toShare/anonymous_ids/AxiomGT1v2.comp_merge.fam
cp AxiomGT1v2.comp_merge.{bim,bed} toShare/anonymous_ids/.
rclone -v --copy-links copy toShare/anonymous_ids remote_UCDavis_GoogleDr:MAF/genotypingV2/anonymous_ids

## Create a version with public_ids
tail -n+2 map_id_sex.tab | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$1,$4}' > public_ids.update.lst
plink2 --bfile AxiomGT1v2.comp_merge --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
       --update-ids public_ids.update.lst \
       --make-bed --output-chr 'chrM' --out toShare/public_ids/AxiomGT1v2.comp_merge.public_ids
rclone -v --copy-links copy toShare/public_ids remote_UCDavis_GoogleDr:MAF/genotypingV2/public_ids

##########################################################
## Identification and removal of duplicate samples
# 1. Identify the replicates
plink2 --bfile AxiomGT1v2.comp_merge --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
       --king-cutoff 0.354 \
       --out AxiomGT1v2.comp_merge.nodup

# 2. compare to the sample planned to be duplicates
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$3]+=1;next}/^Family_ID/{print}{if(a[$3]>1)print}' map_id_sex.tab map_id_sex.tab > dup_samples.tab
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1 FS $2]=1;next}{if(!a[$1 FS $2])print $0}' AxiomGT1v2.comp_merge.nodup.king.cutoff.out.id dup_samples.tab > dup_samples_remain.tab
tail -n+2 dup_samples_remain.tab | cut -f3 | sort | uniq -c | sort -k1,1nr

# 3. confirm the KING kinckship of the failed replicates
grep 094-027376 dup_samples.tab > failed_deDup.lst
plink2 --bfile AxiomGT1v2.comp_merge --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
       --keep failed_deDup.lst --make-king-table \
       --out AxiomGT1v2.comp_merge.failed_deDup

# 4. exclude replicates and failed replicates
plink2 --bfile AxiomGT1v2.comp_merge --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
       --remove <(cat AxiomGT1v2.comp_merge.nodup.king.cutoff.out.id failed_deDup.lst) \
       --make-bed --output-chr 'chrM' --out AxiomGT1v2.comp_merge.deDup
##########################################################
## Check for gender accuracy & remove samples with wrong gender identities
echo "Family_ID Individual_ID computed_sex metadata_sex" | tr ' ' '\t' > gender_conflict.lst
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1 FS $2]=$5;next}{if(a[$1 FS $2] && a[$1 FS $2]!=$5)print $1,$2,$5,a[$1 FS $2]}' map_id_sex.tab AxiomGT1v2.comp_merge.deDup.fam >> gender_conflict.lst
##########################################################
## Update sample IDs in the genotyping files to match the phenotyping files

# 1. upadate the sample IDs to match the grls_ids (e.g. 094-000019)
tail -n+2 map_id_sex.tab | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$1,$3}' > grls_id.update.lst
plink2 --bfile AxiomGT1v2.comp_merge.deDup.sexConfirm --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
       --update-ids grls_id.update.lst \
       --make-bed --output-chr 'chrM' --out AxiomGT1v2.comp_merge.deDup.sexConfirm.grls_ids

# 2. upadate the sample IDs to match the public_ids (e.g. grlsH764T844)
tail -n+2 map_id_sex.tab | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$1,$4}' > public_ids.update.lst
plink2 --bfile AxiomGT1v2.comp_merge.deDup.sexConfirm --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
       --update-ids public_ids.update.lst \
       --make-bed --output-chr 'chrM' --out AxiomGT1v2.comp_merge.deDup.sexConfirm.public_ids
##########################################################
## MAF_uploads
mkdir -p $work_dir/MAF_uploads/output/{setA,setB}/export_plink/
cp output/setA/export_plink/AxiomGT1.bin.{bed,bim,fam} MAF_uploads/output/setA/export_plink/.
cp output/setB/export_plink/AxiomGT1.bin.{bed,bim,fam} MAF_uploads/output/setB/export_plink/.
cp map_id_sex.tab MAF_uploads/.

rclone -v --copy-links copy MAF_uploads remote_UCDavis_GoogleDr:MAF/MAF_uploads
##########################################################
## to run the online tutorial
rclone -v --copy-links copy remote_UCDavis_GoogleDr:MAF/MAF_uploads/output output
rclone -v --copy-links copy remote_UCDavis_GoogleDr:MAF/Phenotypes/full_study_access.zip phenotypes
unzip full_study_access.zip -d phenotypes
