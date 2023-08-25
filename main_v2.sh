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

#conda create -n ngs
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

cd $work_dir/Entire_Cohort_Zip_Files/Array_B/ && unzip cluster_array_B.zip
cd $work_dir/Entire_Cohort_Zip_Files/Array_B/ && unzip cluster_CN_Discovery_array_B.zip
rclone -v --copy-links copy Entire_Cohort_Zip_Files/Array_B/cluster_array_B/AxiomAnalysisSuiteData/AnalysisConfiguration.analysis_settings remote_UCDavis_GoogleDr:MAF/genotyping_array/Entire_Cohort_Zip_Files/Array_B/cluster_array_B/AxiomAnalysisSuiteData
rclone -v --copy-links copy Entire_Cohort_Zip_Files/Array_B/cluster_CN_Discovery_array_B/AxiomAnalysisSuiteData/AnalysisConfiguration.analysis_settings remote_UCDavis_GoogleDr:MAF/genotyping_array/Entire_Cohort_Zip_Files/Array_B/cluster_CN_Discovery_array_B/AxiomAnalysisSuiteData
rclone -v --copy-links copy Entire_Cohort_Zip_Files/Array_B/cluster_array_B/customLib/Axiom_K9HDSNPB_96orMore_Step2.r2.apt-genotype-axiom.AxiomGT1.apt2_mod.xml remote_UCDavis_GoogleDr:MAF/genotyping_array/Entire_Cohort_Zip_Files/Array_B/cluster_array_B/customLib
rclone -v --copy-links copy Entire_Cohort_Zip_Files/Array_B/cluster_CN_Discovery_array_B/customLib/custom.cn_models remote_UCDavis_GoogleDr:MAF/genotyping_array/Entire_Cohort_Zip_Files/Array_B/cluster_CN_Discovery_array_B/customLib


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

## download support file (Analysis library files)
mkdir -p $work_dir/lib/{setA,setB}
cd $work_dir/lib/setA
## download files from: https://www.thermofisher.com/order/catalog/product/550771
## download require users to sign in and accept terms and conditions. Therefore, I had to download locally then upload to the server
## Upload to the server
# 1) Axiom Canine Genotyping Array Set A Analysis Files, r1
# local: scp TFS-Assets_LSG_Support-Files_Axiom_K9HDSNPA_Analysis.r1.zip tahmed@farm.cse.ucdavis.edu:/home/tahmed/MAF/lib/setA/.
unzip -d setA_Analysis TFS-Assets_LSG_Support-Files_Axiom_K9HDSNPA_Analysis.r1.zip
# 2) Current NetAffx Annotation Files: Axiom Canine Genotyping Array Set A Annotations, CSV format
# local: scp TFS-Assets_LSG_Support-Files_Axiom_K9HDSNPA_Annotation.r1_3.csv.zip tahmed@farm.cse.ucdavis.edu:/home/tahmed/MAF/lib/setA/.
unzip TFS-Assets_LSG_Support-Files_Axiom_K9HDSNPA_Annotation.r1_3.csv.zip
# 3) NetAffx Alignment Files: Axiom Canine Genotyping Array Set A BED File
# local: scp TFS-Assets_LSG_Support-Files_Axiom_K9HDSNPA_Annotation.r1_3.bed.zip tahmed@farm.cse.ucdavis.edu:/home/tahmed/MAF/lib/setA/.
unzip TFS-Assets_LSG_Support-Files_Axiom_K9HDSNPA_Annotation.r1_3.bed.zip
# 4) Current NetAffx Annotation Files: Axiom_K9HDSNPA Annotations, SQLite Format, r1
# local: scp Axiom_K9HDSNPA_Annotation.r1_3.db.zip tahmed@farm.cse.ucdavis.edu:/home/tahmed/MAF/lib/setA/.
unzip Axiom_K9HDSNPA_Annotation.r1_3.db.zip


cd $work_dir/lib/setB
## download files from: https://www.thermofisher.com/order/catalog/product/550772
## Upload to the server
# 1) Axiom Canine Genotyping Array Set B Analysis Files, r1
# local: scp TFS-Assets_LSG_Support-Files_Axiom_K9HDSNPB_Analysis.r1.zip tahmed@farm.cse.ucdavis.edu:/home/tahmed/MAF/lib/setB/.
unzip -d setB_Analysis TFS-Assets_LSG_Support-Files_Axiom_K9HDSNPB_Analysis.r1.zip
# 2) Current NetAffx Annotation Files: Axiom Canine Genotyping Array Set B Annotations, CSV format
# local: scp TFS-Assets_LSG_Support-Files_Axiom_K9HDSNPB_Annotation.r1_3.csv.zip tahmed@farm.cse.ucdavis.edu:/home/tahmed/MAF/lib/setB/.
unzip TFS-Assets_LSG_Support-Files_Axiom_K9HDSNPB_Annotation.r1_3.csv.zip
# 3) NetAffx Alignment Files: Axiom Canine Genotyping Array Set B BED File
# local: scp TFS-Assets_LSG_Support-Files_Axiom_K9HDSNPB_Annotation.r1_3.bed.zip tahmed@farm.cse.ucdavis.edu:/home/tahmed/MAF/lib/setB/.
unzip TFS-Assets_LSG_Support-Files_Axiom_K9HDSNPB_Annotation.r1_3.bed.zip
# 4) Current NetAffx Annotation Files: Axiom_K9HDSNPB Annotations, SQLite Format, r1
# local: scp Axiom_K9HDSNPB_Annotation.r1_3.db.zip tahmed@farm.cse.ucdavis.edu:/home/tahmed/MAF/lib/setB/.
unzip Axiom_K9HDSNPB_Annotation.r1_3.db.zip


## download NEW support file (Analysis library files). We got directly from Thermo
module load rclone
mkdir -p $work_dir/lib2/{setA,setB}
rclone -v --copy-links copy remote_UCDavis_GoogleDr:MAF/New_axiom_files/Axiom_K9HDSNPA_Analysis.r2.zip /home/tahmed/MAF/lib2/setA/
rclone -v --copy-links copy remote_UCDavis_GoogleDr:MAF/New_axiom_files/Axiom_K9HDSNPA.na35.r2.a4.annot.db /home/tahmed/MAF/lib2/setA/

rclone -v --copy-links copy remote_UCDavis_GoogleDr:MAF/New_axiom_files/Axiom_K9HDSNPB_Analysis.r2.zip /home/tahmed/MAF/lib2/setB/
rclone -v --copy-links copy remote_UCDavis_GoogleDr:MAF/New_axiom_files/Axiom_K9HDSNPB.na35.r2.a4.annot.db /home/tahmed/MAF/lib2/setB/


mamba create -n sqlite -c conda-forge sqlite
conda activate sqlite
## list the tables
sqlite3 $work_dir/lib2/setA/Axiom_K9HDSNPA.na35.r2.a4.annot.db ".tables"
#Annotations     CdfInformation  Chromosome      Information     Localization


## Extract traget tables
cd $work_dir/lib2/setA
sqlite3 -header -csv Axiom_K9HDSNPA.na35.r2.a4.annot.db "SELECT * FROM Annotations;" > Axiom_K9HDSNPA_Annotation.r2_3.csv


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
cat $work_dir/pedigree/sample_ped_infoA_cohort.txt > $work_dir/All_pedigree/sample_ped_infoA_cohort.txt
tail -n+2 $work_dir/ped_pilot/sample_ped_infoA_cohort.txt >> $work_dir/All_pedigree/sample_ped_infoA_cohort.txt

cat $work_dir/pedigree/sample_ped_infoB_cohort.txt > $work_dir/All_pedigree/sample_ped_infoB_cohort.txt
tail -n+2 $work_dir/ped_pilot/sample_ped_infoB_cohort.txt >> $work_dir/All_pedigree/sample_ped_infoB_cohort.txt
########
## Currently this section is NOT IN USE
## prep the inbred file
cd $work_dir
echo "cel_files inbred_het_penalty" | tr ' ' '\t' > inbred_test.txt
for f in cel/*/*.CEL;do
echo $(basename "$f") 4 | tr ' ' '\t';done >> inbred_test.txt


## prep the gender file
cd $work_dir
echo "cel_files gender" | tr ' ' '\t' > gender_test.txt
for f in cel/*/*.CEL;do
echo $(basename "$f");done | grep -Fwf - pedigree/sample_ped_infoB.txt | awk -F"\t" 'BEGIN{FS=OFS="\t"}{if($6==1)print $1,"male";else if($6==2)print $1,"female";else print $1,"unknown"}' >> gender_test.txt

#######
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
set_Analysis=$work_dir/"lib/setA/setA_Analysis/"
set_output=$work_dir/"output/setA"
qc_xml="$set_Analysis"/Axiom_K9HDSNPA.r1.apt-geno-qc.AxiomQC1.xml
genotype1_xml="$set_Analysis"/Axiom_K9HDSNPA_96orMore_Step1.r1.apt-genotype-axiom.AxiomGT1.apt2.xml
genotype2_xml="$set_Analysis"/Axiom_K9HDSNPA_96orMore_Step2.r1.apt-genotype-axiom.AxiomGT1.apt2.xml
specialSNPs="$set_Analysis"/Axiom_K9HDSNPA.r1.specialSNPs
ps2snp="$set_Analysis"/Axiom_K9HDSNPA.r1.ps2snp_map.ps
annot_db="$set_Analysis"/../Axiom_K9HDSNPA.r1_3.20160810.annot.db
annot_csv="$set_Analysis"/../Axiom_K9HDSNPA_Annotation.r1_3.csv
pedInfo=$work_dir/pedigree/sample_ped_infoA_cohort.txt

## Array_B
cel_PATH=$work_dir/"CEL_Files/Array_B" #$work_dir/"CEL_Files/Array_A"
set_Analysis=$work_dir/"lib/setB/setB_Analysis/" #$work_dir/"lib/setA/setA_Analysis/"
set_output=$work_dir/"output/setB" #$work_dir/"output/setA"
qc_xml="$set_Analysis"/Axiom_K9HDSNPB.r1.apt-geno-qc.AxiomQC1.xml #"$set_Analysis"/Axiom_K9HDSNPA.r1.apt-geno-qc.AxiomQC1.xml
genotype1_xml="$set_Analysis"/Axiom_K9HDSNPB_96orMore_Step1.r1.apt-genotype-axiom.AxiomGT1.apt2.xml #"$set_Analysis"/Axiom_K9HDSNPA_96orMore_Step1.r1.apt-genotype-axiom.AxiomGT1.apt2.xml
genotype2_xml="$set_Analysis"/Axiom_K9HDSNPB_96orMore_Step2.r1.apt-genotype-axiom.AxiomGT1.apt2.xml
#"$set_Analysis"/Axiom_K9HDSNPA_96orMore_Step2.r1.apt-genotype-axiom.AxiomGT1.apt2.xml
specialSNPs="$set_Analysis"/Axiom_K9HDSNPB.r1.specialSNPs #"$set_Analysis"/Axiom_K9HDSNPA.r1.specialSNPs
ps2snp="$set_Analysis"/Axiom_K9HDSNPB.r1.ps2snp_map.ps #"$set_Analysis"/Axiom_K9HDSNPA.r1.ps2snp_map.ps
annot_db="$set_Analysis"/../Axiom_K9HDSNPB.r1_3.20160810.annot.db #"$set_Analysis"/../Axiom_K9HDSNPA.r1_3.20160810.annot.db
annot_csv="$set_Analysis"/../Axiom_K9HDSNPB_Annotation.r1_3.csv
pedInfo=$work_dir/pedigree/sample_ped_infoB_cohort.txt #$work_dir/pedigree/sample_ped_infoA_cohort.txt

##########
## Start the Best Practice analysis
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
(echo cel_files; awk -v p="$cel_PATH/" -F"\t" 'FNR==NR{if($18>0.82)a[p$1]=1;next}{if(a[$1])print}' <(grep -v "^#"  "$set_output"/apt-geno-qc.txt) "$set_output"/cel_list1.txt) > "$set_output"/cel_list2.txt  ## arrayA (0) // arrayB (650_E12.CEL and 658_B12.CEL)
 

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
grep -v "^#"  "$set_output"/step1/AxiomGT1.report.txt | tail -n+2 | cut -f1,3 | cut -d"-" -f4 | awk -F"[_\t]" '{if($3<97)print $1,$3}' ## Array_A (9) // Array_B (42)
grep -v "^#"  "$set_output"/step1/AxiomGT1.report.txt | tail -n+2 | cut -f1,3 | cut -d"-" -f4 | awk -F"[_\t]" '{if($3<95)print $1,$3}' ## Array_A (5) // Array_B (19)
grep -v "^#"  "$set_output"/step1/AxiomGT1.report.txt | tail -n+2 | cut -f1,3 | cut -d"-" -f4 | awk -F"[_\t]" '{if($3<95)print $1}' | sort | uniq -c ## Array_A (131,438x2, 137x1) // Array_B (655x6, 493x3, 494,496,650x2, 491,495,677,679x1)
grep -v "^#"  "$set_output"/step1/AxiomGT1.report.txt | tail -n+2 | cut -f1,3 | cut -d"-" -f4 | awk -F"[_\t]" '{if($3>95 && $3<97)print $1,$3}' ## Array_A (129x2,131,134x1) // Array_B (493x5, 495x4, 679,066,492,497x2, 658,666,064,070,151,496x1) => These are the plate that will likely fail the average plate QC
# Tamer: I will use 95% ae a cutoff and accordingly I will tolerate lower average plate QC
(echo cel_files; awk -v p="$cel_PATH/" -F"\t" 'FNR==NR{if($3>95)a[p$1]=1;next}{if(a[$1])print}' <(grep -v "^#"  "$set_output"/step1/AxiomGT1.report.txt) "$set_output"/cel_list1.txt) > "$set_output"/cel_list3.txt

## 6: QC the plates
tail -n+2 "$set_output"/cel_list1.txt | sed 's|.*/||' | cut -d"-" -f4 | cut -d "_" -f1 | uniq -c | sed 's/^ *[^0-9]//g' > "$set_output"/samples_on_plate
tail -n+2 "$set_output"/cel_list3.txt | sed 's|.*/||' | cut -d"-" -f4 | cut -d "_" -f1 | uniq -c | sed 's/^ *[^0-9]//g' > "$set_output"/samples_pass_plate
awk 'BEGIN{FS=" "}FNR==NR{a[$2]=$1;next}{if(a[$2])print $2,$1/a[$2] * 100}' "$set_output"/samples_on_plate "$set_output"/samples_pass_plate  > "$set_output"/plate_pass_rate ## As expected, plate 129 on Array_A (131, 137, 438) have the failing samples but all are above 97.5 // Array_B (9 plates have lost samples)
grep -v "^#"  "$set_output"/step1/AxiomGT1.report.txt | tail -n+2 | cut -f1,3 | cut -d"-" -f4 | awk -F"[_\t]" 'BEGIN{CONVFMT="%.3g"}{if($3>95){sumQC[$1]+=$3;count[$1]++}}END{for(i in sumQC)print i,sumQC[i]/count[i]}' > "$set_output"/plate_aveQC ## Array_A: 129 (ave=98.3965). This is the only plate failing 2 samples // Array_B: 679 (ave=95.8125), 667 (ave=97.675), 678 (ave=97.775)

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
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2;next}{if($6=="1")$6="male";else if($6=="2")$6="female";else $6="unknown"; if(a[$1])print $1,"S"$3,$6,a[$1]}' <(grep -v "^#" "$set_output"/step2/AxiomGT1.report.txt | tail -n+2) <(tail -n+2 $pedInfo) >> "$set_output"/step2/gender_check
cat "$set_output"/step2/gender_check | awk 'BEGIN{FS=OFS="\t"}{if($3!="unknown" && $4=="unknown")print}' | wc -l #A50/B110
cat "$set_output"/step2/gender_check | awk 'BEGIN{FS=OFS="\t"}{if($3=="unknown" && $4=="unknown")print}' | wc -l #0
cat "$set_output"/step2/gender_check | awk 'BEGIN{FS=OFS="\t"}{if(($3=="unknown") || ($4!="unknown" && $3!=$4))print}' | sort -k2,2 | tr '\t' ',' > "$set_output"/step2/gender_dif

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

## Array_A ## 401504 output/setA/SNPolisher/Recommended.ps
## Array_B ## 548289 output/setB/SNPolisher/Recommended.ps

## Export
## create map files for CEL file names and sample IDs
echo "Sample Filename|Alternate Sample Name" | tr '|' '\t' > "$set_output"/info_names.txt
tail -n+2 "$pedInfo" | awk 'BEGIN{FS=OFS="\t"}{print $1,"S"$3}' >> "$set_output"/info_names.txt

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

wc -l "$set_output"/SNPolisher/Recommended.ps ## 401494 // 548289
grep -v "^#" $set_output/export/AxiomGT1.vcf | wc -l  ## 401493 // 548288


cd "$set_output"/export/
## QC for SNPs without ref/alt
## extract all SNPs without ref/alt & select those among recommended
#grep -v "^#" $annot_csv | sed 's/^"//;s/\",\"/|/g;s/"$//;' | awk -F"|" '{if($19=="---")print}' > noREF.csv
grep -Fwf "$set_output"/SNPolisher/Recommended.ps <(grep -v "^#" $annot_csv | sed 's/^"//;s/\",\"/|/g;s/"$//;' | awk -F"|" '{if($19=="---")print}') > noREF_rec.csv #1409 // 0
rclone -v --copy-links copy noREF_rec.csv remote_UCDavis_GoogleDr:MAF/genotyping_array ## Thermo
cat noREF_rec.csv | cut -d"|" -f4 | sort | uniq -c > noREF_rec.chr ## arrA (704 Unknown, 421 Unplaced contigs, 177Y, 25X, 82 all other chr)
## Split into variants with known Chr/Pos and variants without
cat noREF_rec.csv | awk -F"|" '{if($4!="---" && $5!="---")print}' > noREF_rec_known.csv
cat noREF_rec.csv | awk -F"|" '{if($4=="---" || $5=="---")print}' > noREF_rec_unknown.csv
## known: Fix chr symbol THEN prepare possible alleles & haplotypes THEN sort by chromosome
cat noREF_rec_known.csv | awk 'BEGIN{FS="|";OFS="\t"}{if($4=="Un"){z=$15;sub(/^{/,"",z);split(z,a,"_");$4=a[1]"_"a[2];}else if($4=="MT")$4="chrM";else $4="chr"$4;print $2,$4,$5,$6,$7}' | tr '[/]' '\t\t\t' | awk 'BEGIN{FS=OFS="\t"}{if($6=="-"){$6="";a1=substr($5,length($5));a2=a1$7;} else{a1=$6;a2=$7}print $1,$2,$3,$4,a1,a2,length($5),$5$6$8,$5$7$8}' | sort -k2,2 > noREF_rec_knownCor.tab ## 704
## Unknown: prepare possible alleles & haplotypes
cat noREF_rec_unknown.csv | awk 'BEGIN{FS="|";OFS="\t"}{print $2,$4,$5,$6,$7}' | tr '[/]' '\t\t\t' | awk 'BEGIN{FS=OFS="\t"}{if($6=="-"){$6="";a1=substr($5,length($5));a2=a1$7;} else{a1=$6;a2=$7}print $1,$2,$3,$4,a1,a2,length($5),$5$6$8,$5$7$8}' > noREF_rec_unknownCor.tab ## 705 (This must include the variant that already failed to export to VCF)
## Prep fasta files
cat noREF_rec_knownCor.tab noREF_rec_unknownCor.tab | awk 'BEGIN{FS="\t";OFS="\n";}{print ">"$1":a1",$8,">"$1":a2",$9}' > noREF_rec_knownCor.fa

## align by BWA
bwa mem $genome_dir/bwaIndex/canFam3_unwrap.fa noREF_rec_knownCor.fa > noREF_rec_knownCor.sam
grep -v "^@" noREF_rec_knownCor.sam | awk -F"\t" '{if($12=="NM:i:0")print $1}' > temp.align
cat temp.align | tr ':' '\t' | cut -f1 | sort | uniq -c | awk '{if($1=="1")print $2}' | grep -Fwf - temp.align | grep -Fwf - noREF_rec_knownCor.sam > noREF_rec_knownCor.uniq.sam ## 687
cat noREF_rec_knownCor.uniq.sam | awk '{if($5>=45)print}' > noREF_rec_knownCor.sig.sam ## 425 (204 Unplaced contigs, 171Y, 24X, 26 all other chr)


## Update the VCF chromosomes to match canFam3
# change 1,2,3,.....X,Y to chr1,chr2,...... chrY
# change MT to chrM
# change Un to chrUn_xxxxxxx
# ?remove UNKNOWNCHR
grep -v "^#" AxiomGT1.vcf | cut -f1 | uniq -c > chr.lst

for x in {1..38} X Y;do echo $x chr$x | tr ' ' '\t';done > allChr.map
echo MT chrM | tr ' ' '\t' >> allChr.map
for x in {1..38} X Y M;do echo "##contig=<ID=chr$x>";done > contig_ids

grep "chrUn_" $annot_csv | sed 's/^"//;s/\",\"/|/g;s/"$//;' | cut -d"|" -f2,15 | sed 's/{//' | cut -d"_" -f1,2 | tr '|' '\t' > chrUn.map
grep "^Un" AxiomGT1.vcf | cut -f3 | grep -Fwf - chrUn.map > chrUn.map.found
cat chrUn.map.found >> allChr.map
cat chrUn.map.found | cut -f2 | sort | uniq | while read f;do echo "##contig=<ID=$f>";done >> contig_ids

cat AxiomGT1.vcf | awk -v file=contig_ids 'BEGIN{contig=1;}/^##[^contig]/{print}/^##contig/{if(contig){while((getline<file) > 0)print;contig=0;}}/^#CH/{print}!/^#/{exit}' > AxiomGT1v2.vcf
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2;next}!/#/{if(a[$1]){$1=a[$1];print}else if(a[$3]){$1=a[$3];print}}' allChr.map AxiomGT1.vcf >> AxiomGT1v2.vcf ## 400789 // 548288

grep -v "^#" AxiomGT1v2.vcf | cut -f1 | uniq -c > chrv2.lst


#conda create -n equineSNP
#conda install -c bioconda plink
#conda install -c bioconda plink2
#conda install -c bioconda bcftools
#conda install -c bioconda eigensoft
#conda install -c bioconda vcftools
conda activate equineSNP
module load rclone #On Nov 2022: Module rclone/1.53.3 loaded


## check the VCF integrity
conda activate equineSNP
bcftools sort -o AxiomGT1v2_sorted.vcf AxiomGT1v2.vcf
bcftools norm -c ws -f $canFam3_ref_unwrap AxiomGT1v2_sorted.vcf 1> AxiomGT1v2.check.vcf 2> AxiomGT1v2.check.log
# Array A
#Lines   total/split/realigned/skipped:  400789/0/141/0
#REF/ALT total/modified/added:   400789/112/1191
# Array B
#Lines   total/split/realigned/skipped:  548288/0/0/0
#REF/ALT total/modified/added:   548288/0/0

diff <(grep "^#" AxiomGT1v2.vcf) <(grep "^#" AxiomGT1v2.check.vcf) > check.dif
diff <(grep -v "^#" AxiomGT1v2.vcf | cut -f1-5) <(grep -v "^#" AxiomGT1v2.check.vcf | cut -f1-5) > check.dif2
cat check.dif2 | grep "^[<>]" | sort -k2,3 -k1,1 > check.dif2.sorted
cat check.dif2.sorted | cut -f2-  | uniq | cut -f2 | uniq -c | awk '{if($1==2)print $2}' | grep -Fwf - check.dif2.sorted > check.dif2.sorted.uniq
cat check.dif2.sorted.uniq | cut -d" " -f2 | paste - - | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$9,$10}' > check.dif.map
cat check.dif.map | awk 'BEGIN{FS=OFS="\t"}{if($4==$7 && $5==$6)print}' | wc -l ## 112
cat check.dif.map | awk 'BEGIN{FS=OFS="\t"}{if($4=="N" && $5==".")print}' | wc -l ## 700
cat check.dif.map | awk 'BEGIN{FS=OFS="\t"}{if($4!="N" && $4!=$7)print}' | wc -l ## 491


## Exclude variants without known ref/alt alleles
#grep -v "^#" AxiomGT1v2.vcf | awk 'BEGIN{FS=OFS="\t"}{if($4=="N")print $1,$2}' > noPos.lst ## 705
grep -v "^#" AxiomGT1v2.vcf | awk 'BEGIN{FS=OFS="\t"}{if($4=="N")print $3}' > noPos.lst ## 705
## Exclude variants with annotation errors
#cat check.dif.map | awk 'BEGIN{FS=OFS="\t"}{if($4!="N" && $4!=$7)print $1,$2}' > annError.lst ## 491
cat check.dif.map | awk 'BEGIN{FS=OFS="\t"}{if($4!="N" && $4!=$7)print $3}' > annError.lst ## 491

#vcftools --vcf AxiomGT1v2.check.vcf --exclude-positions <(cat noPos.lst annError.lst) --recode --out AxiomGT1v2.hiQ
cat noPos.lst annError.lst | grep -vFwf - AxiomGT1v2.check.vcf > AxiomGT1v2.hiQ.vcf
grep -v "^#" AxiomGT1v2.hiQ.vcf | wc -l ## 399593 // 548288

## Generate lists of annotation errors for Thermo
cat annError.lst | grep -Fwf - AxiomGT1v2.vcf | cut -f3,8 > annErrorProbes.lst
echo "chr pos probe-set ref alt CF3-ref marker" | tr ' ' '\t' > annErrorProbes.map
awk 'BEGIN{FS=OFS="\t"}FNR==NR{if($4!="N" && $4!=$7)a[$3]=$0;next}{if(a[$1])print a[$1],$2}' check.dif.map annErrorProbes.lst | sed 's/probeset_id=//' | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$8,$6,$7,$4,$3}' >> annErrorProbes.map
rclone -v --copy-links copy annErrorProbes.map remote_UCDavis_GoogleDr:MAF/genotyping_array

cat check.dif.map | awk 'BEGIN{FS=OFS="\t"}{if($4==$7 && $5==$6)print $3}' > strandSwab.lst
cat strandSwab.lst | grep -Fwf - AxiomGT1v2.vcf | cut -f3,8 > strandSwabProbes.lst
echo "chr pos probe-set ref alt CF3-ref marker" | tr ' ' '\t' > strandSwabProbes.map
awk 'BEGIN{FS=OFS="\t"}FNR==NR{if($4==$7 && $5==$6)a[$3]=$0;next}{if(a[$1])print a[$1],$2}' check.dif.map strandSwabProbes.lst | sed 's/probeset_id=//' | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$8,$6,$7,$4,$3}' >> strandSwabProbes.map
rclone -v --copy-links copy strandSwabProbes.map remote_UCDavis_GoogleDr:MAF/genotyping_array



## There are different SNPs with the same position on the same array
## Keep one of them which has higher genotyping rate
#grep -v "^#" AxiomGT1v2.hiQ.vcf | awk '{print $3}' | sort > ids.lst
grep -v "^#" AxiomGT1v2.hiQ.vcf | awk '{print $1"_"$2}' | sort > pos.lst
cat pos.lst | uniq -c | awk '{if($1>1)print $2}' > pos_dup.lst # 2862 // 0
cat pos.lst | uniq -c | awk '{if($1>2)print $2}' ## Not more than 2 IDs per position
awk -F"\t" 'FNR==NR{a[$1]=1;next}{if(a[$1"_"$2])print}' pos_dup.lst AxiomGT1v2.hiQ.vcf | cut -f1,2 | paste - - | awk '{if($1!=$3 || $2!=$4)print}' ## nothing (safe to paste)
echo "CHROM POS ID1 NS1 MAF1 AC_Het1 AC_Hom1 ExcHet1 ID2 NS2 MAF2 AC_Het2 AC_Hom2 ExcHet2" | tr ' '  '\t' > pos_dup.ann
awk -F"\t" 'FNR==NR{a[$1]=1;next}/^#/{print}{if(a[$1"_"$2])print}' pos_dup.lst AxiomGT1v2.hiQ.vcf | bcftools +fill-tags | bcftools query -f'%CHROM\t%POS\t%ID\t%NS\t%MAF\t%AC_Het\t%AC_Hom\t%ExcHet\n' | paste - - | cut -f1-8,11-16 >> pos_dup.ann #pos_dup.vcf
tail -n+2 pos_dup.ann | awk -F"\t" '{if($4<$10)print $3;else print $9}' > pos_dup_ToExclude.lst ## 2862 // 0
cat pos_dup_ToExclude.lst | grep -vFwf - AxiomGT1v2.hiQ.vcf > AxiomGT1v2.hiQ2.vcf
grep -v "^#" AxiomGT1v2.hiQ2.vcf | awk '{print $1"_"$2}' | sort > pos2.lst

## Thermo
tail -n+2 pos_dup.ann | cut -f3 | grep -Fwf - AxiomGT1v2.vcf | cut -f1-3,8 > ID1.lst
tail -n+2 pos_dup.ann | cut -f9 | grep -Fwf - AxiomGT1v2.vcf | cut -f1-3,8 > ID2.lst
echo "chr pos probe-set1 affy_snp_id1 probe-set2 affy_snp_id2" | tr ' ' '\t' > dupMarkersWithinArrA.map
paste ID1.lst ID2.lst | awk 'BEGIN{FS=OFS="\t"}{if($1==$5 && $2==$6)print $1,$2,$4,$3,$8,$7}' | sed 's/probeset_id=//g'  >> dupMarkersWithinArrA.map
rclone -v --copy-links copy dupMarkersWithinArrA.map remote_UCDavis_GoogleDr:MAF/genotyping_array


########
cd $work_dir
grep "^#CHROM" output/setA/export/AxiomGT1.vcf | tr '\t' '\n' | tail -n+10 > A.samples
grep "^#CHROM" output/setB/export/AxiomGT1.vcf | tr '\t' '\n' | tail -n+10 > B.samples
comm -12 <(sort A.samples) <(sort B.samples) > shared.samples ## 3146
comm -13 <(sort A.samples) <(sort B.samples) > Bonly.samples ## 1
comm -23 <(sort A.samples) <(sort B.samples) > Aonly.samples ## 17

## Assess duplicates and genotype concordance
comm -12 output/set[AB]/export/pos2.lst > dup_pos # 31083
awk -F"\t" 'FNR==NR{a[$1]=1;next}/^#/{print}{if(a[$1"_"$2])print}' dup_pos output/setA/export/AxiomGT1v2.hiQ2.vcf > dup_A.vcf
awk -F"\t" 'FNR==NR{a[$1]=1;next}/^#/{print}{if(a[$1"_"$2])print}' dup_pos output/setB/export/AxiomGT1v2.hiQ2.vcf > dup_B.vcf
bgzip dup_A.vcf && tabix -p vcf dup_A.vcf.gz
bgzip dup_B.vcf && tabix -p vcf dup_B.vcf.gz
bcftools stats --samples-file shared.samples dup_A.vcf.gz dup_B.vcf.gz > dup_stats
grep -A 3147 "Genotype concordance by sample (SNPs)" dup_stats > dup_Geno_concord
tail -n+3 dup_Geno_concord | cut -f11 | awk '{a+=1;b+=$1}END{print b/a}' ## 0.954574

echo "CHROM POS ID1 NS1 MAF1 AC_Het1 AC_Hom1 ExcHet1" | tr ' '  '\t' > dup_A.info
cat dup_A.vcf | bcftools +fill-tags | bcftools query -f'%CHROM\t%POS\t%ID\t%NS\t%MAF\t%AC_Het\t%AC_Hom\t%ExcHet\n' >> dup_A.info
echo "CHROM POS ID2 NS2 MAF2 AC_Het2 AC_Hom2 ExcHet2" | tr ' '  '\t' > dup_B.info
cat dup_B.vcf | bcftools +fill-tags | bcftools query -f'%CHROM\t%POS\t%ID\t%NS\t%MAF\t%AC_Het\t%AC_Hom\t%ExcHet\n' >> dup_B.info

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
tail -n+2 dup_B.info | awk -F"\t" '{print $3}' > pos_dupAcrossArr_ToExclude.lst
cat pos_dupAcrossArr_ToExclude.lst | grep -vFwf - output/setB/export/AxiomGT1v2.hiQ2.vcf > output/setB/export/AxiomGT1v2.hiQ3.vcf

## Thermo
cat dup_AvsB.info | awk -F"\t" '{if($3!=$9)print}' | wc -l ## 1567
cat dup_AvsB.info | cut -f1-3,9 > dup_AvsB.map
rclone -v --copy-links copy dup_AvsB.map remote_UCDavis_GoogleDr:MAF/genotyping_array
rclone -v --copy-links copy dup_A.vcf.gz remote_UCDavis_GoogleDr:MAF/genotyping_array
rclone -v --copy-links copy dup_B.vcf.gz remote_UCDavis_GoogleDr:MAF/genotyping_array

## prep final array files to final merge
bgzip output/setA/export/AxiomGT1v2.hiQ2.vcf && tabix -p vcf output/setA/export/AxiomGT1v2.hiQ2.vcf.gz
bgzip output/setB/export/AxiomGT1v2.hiQ3.vcf && tabix -p vcf output/setB/export/AxiomGT1v2.hiQ3.vcf.gz

bcftools view --samples-file shared.samples output/setA/export/AxiomGT1v2.hiQ2.vcf.gz | bgzip -c > output/setA/export/AxiomGT1v2.shared.vcf.gz
tabix -p vcf output/setA/export/AxiomGT1v2.shared.vcf.gz
bcftools view --samples-file shared.samples output/setB/export/AxiomGT1v2.hiQ3.vcf.gz | bgzip -c > output/setB/export/AxiomGT1v2.shared.vcf.gz
tabix -p vcf output/setB/export/AxiomGT1v2.shared.vcf.gz

bcftools view --samples-file Aonly.samples output/setA/export/AxiomGT1v2.hiQ2.vcf.gz | bgzip -c > output/setA/export/AxiomGT1v2.Aonly.vcf.gz
tabix -p vcf output/setA/export/AxiomGT1v2.Aonly.vcf.gz
bcftools view --samples-file Bonly.samples output/setB/export/AxiomGT1v2.hiQ3.vcf.gz | bgzip -c > output/setB/export/AxiomGT1v2.Bonly.vcf.gz
tabix -p vcf output/setB/export/AxiomGT1v2.Bonly.vcf.gz


bcftools concat --allow-overlaps output/setA/export/AxiomGT1v2.shared.vcf.gz output/setB/export/AxiomGT1v2.shared.vcf.gz | bgzip -c > AxiomGT1v2.merge.vcf.gz
tabix -p vcf AxiomGT1v2.merge.vcf.gz

bcftools norm -c ws -f $canFam3_ref_unwrap AxiomGT1v2.merge.vcf.gz 1> AxiomGT1v2.merge.check.vcf 2> AxiomGT1v2.merge.check.log
rm AxiomGT1v2.merge.check.vcf
#Lines   total/split/realigned/skipped:  913936/0/0/0
#REF/ALT total/modified/added:   913936/0/0

bcftools merge -m none AxiomGT1v2.merge.vcf.gz output/setA/export/AxiomGT1v2.Aonly.vcf.gz output/setB/export/AxiomGT1v2.Bonly.vcf.gz | bcftools +fill-tags | bgzip -c > AxiomGT1v2.merge2.vcf.gz
tabix -p vcf AxiomGT1v2.merge2.vcf.gz


## Transform all files to map/ped
#zcat AxiomGT1v2.merge.vcf.gz | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf AxiomGT1v2.merge.vcf.gz --double-id --chr-set 38 no-y no-xy --allow-extra-chr --allow-no-sex --make-bed --out "AxiomGT1v2.merge"
# Total genotyping rate is 0.997573.
# 913936 variants and 3146 samples pass filters and QC.

## Calc of IBS distance for pruning of sample duplicates
## I am generating a distance matrix using --distance but we can not run this on all variants therefore I am filtrating on low --geno cutoff
plink --bfile AxiomGT1v2.merge --distance square allele-ct ibs 1-ibs --chr-set 38 no-y no-xy --allow-extra-chr --allow-no-sex --geno 0.01 --out AxiomGT1v2.merge.distance
# 42981 variants removed due to missing genotype data (--geno).
# 870955 variants and 3146 samples pass filters and QC.
# Excluding 28183 variants on non-autosomes from distance matrix calc.

# tranform the distance matrix into one colum
cat AxiomGT1v2.merge.distance.mibs | awk 'BEGIN{FS="\t";OFS="\n";}{for (i = 1; i <= NF; i++)print $i}' > temp.gen.1
# Prep the ids
awk 'BEGIN{OFS="\t"}FNR==NR{a[FNR]=$1;next}{for (i in a)print $1,a[i]}'  <(cat AxiomGT1v2.merge.distance.mdist.id | tr '\t' '.') <(cat AxiomGT1v2.merge.distance.mdist.id | tr '\t' '.') > temp.gen.2
# create a pairwise table for significantly similar samples
paste temp.gen.2 temp.gen.1 | awk 'BEGIN{FS=OFS="\t"}{if($3>0.95 && $1!=$2)print}' > ibs.gen.95_temp
cat ibs.gen.95_temp | awk 'BEGIN{FS=OFS="\t"}{if(!a[$2"."$1] && $3!="nan"){a[$1"."$2]=1;print}}' > ibs.gen.95 ## all replicats are found except S027376
# paste temp.gen.2 temp.gen.1 | grep S027376_1 | grep S027376_2 ## 0.850015

# generate a histogram
awk -v size=0.02 'BEGIN{bmin=bmax=0}{ b=int($1/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/2 }' <(paste temp.gen.2 temp.gen.1 | awk 'BEGIN{FS=OFS="\t"}{if($1!=$2 && $3!="nan")print $3}') > AxiomGT1v2.merge.distance.mibs.2per.hist
awk -v size=0.01 'BEGIN{bmin=bmax=0}{ b=int($1/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/2 }' <(paste temp.gen.2 temp.gen.1 | awk 'BEGIN{FS=OFS="\t"}{if($1!=$2 && $3!="nan")print $3}') > AxiomGT1v2.merge.distance.mibs.1per.hist


## endode the samples in the PLINK files
mv AxiomGT1v2.merge.fam AxiomGT1v2.merge.fam_original
cat AxiomGT1v2.merge.fam_original | awk '{$1=$2=NR;print}' > AxiomGT1v2.merge.fam

## upload the files to GDrive
module load rclone
rclone -v --copy-links copy AxiomGT1v2.merge.fam remote_UCDavis_GoogleDr:MAF/genotypingV1/
rclone -v --copy-links copy AxiomGT1v2.merge.bim remote_UCDavis_GoogleDr:MAF/genotypingV1/
rclone -v --copy-links copy AxiomGT1v2.merge.bed remote_UCDavis_GoogleDr:MAF/genotypingV1/


## generate a copy of the PLINK files using the public IDs
mv AxiomGT1v2.merge.fam AxiomGT1v2.merge.fam_encoded
rclone -v --copy-links copy remote_UCDavis_GoogleDr:MAF/source_files/dc_id.csv .
tail -n+2 dc_id.csv | sed -e "s/\r//g" | sed 's/\"//g' | sed 's/094-/S/' | tr ',' ' ' | cut -d" " -f2,3 > dc_id.txt ## I got the file "dc_id.csv" from Brenna
awk 'FNR==NR{a[$2]=$1;next}{if(a[$1])$1=$2=a[$1];print}' dc_id.txt <(cat AxiomGT1v2.merge.fam_original | sed 's/_1 / /g') > AxiomGT1v2.merge.fam

cat dc_id.txt | cut -d" " -f1 | sort > in.lst
cat AxiomGT1v2.merge.fam | grep -v ^S | cut -d" " -f1 | sort > out.lst
comm -32 in.lst out.lst > miss.lst
cat miss.lst | grep -Fwf - dc_id.txt > miss.txt

grep ^S AxiomGT1v2.merge.fam | cut -d" " -f1,2 > exclude.lst
conda activate equineSNP
plink --bfile AxiomGT1v2.merge --chr-set 38 no-y no-xy --allow-extra-chr --allow-no-sex --remove exclude.lst --make-bed --out AxiomGT1v2.merge.filtered

module load rclone
rclone -v --copy-links copy AxiomGT1v2.merge.filtered.fam remote_UCDavis_GoogleDr:MAF/genotypingV1_publicIDs/
rclone -v --copy-links copy AxiomGT1v2.merge.filtered.bim remote_UCDavis_GoogleDr:MAF/genotypingV1_publicIDs/
rclone -v --copy-links copy AxiomGT1v2.merge.filtered.bed remote_UCDavis_GoogleDr:MAF/genotypingV1_publicIDs/
rclone -v --copy-links copy dc_id.csv remote_UCDavis_GoogleDr:MAF/genotypingV1_publicIDs/

## A code to map back from public_ids (e.g. grlsH764T844) to grls_ids (e.g. 094-000019)
cp AxiomGT1v2.merge.filtered.fam AxiomGT1v2.merge.filtered.fam_backup
awk 'FNR==NR{a[$2]=$3;next}{if(a[$1])$1=$2=a[$1];print}' <(tail -n+2 dc_id.csv | sed -e "s/\r//g" | sed 's/\"//g' | tr ',' ' ') AxiomGT1v2.merge.filtered.fam_backup > AxiomGT1v2.merge.filtered.fam
rm AxiomGT1v2.merge.filtered.fam_backup


###############
## Fix the FAM file of the pilot genotyping data
mkdir pilot && cd pilot
module load rclone
rclone -v --copy-links copy remote_UCDavis_GoogleDr:MAF/source_files/K9_1.1M_merge.fam_fromBrenna .
cat K9_1.1M_merge.fam_fromBrenna | awk 'BEGIN{FS=OFS="\t"}{split($2,a,"A");$2=a[1];x=length($2);if(x==6)$2="S"$2;else if(x==5)$2="S0"$2;else if(x==4)$2="S00"$2;else if(x==3)$2="S000"$2; print}' | tr '\t' ' ' > K9_1.1M_merge.fam_correctIDs
awk 'FNR==NR{a[$2]=$1;next}{if(a[$2])$2=a[$2];print}' ../dc_id.txt K9_1.1M_merge.fam_correctIDs > K9_1.1M_merge.fam


rclone -v --copy-links copy K9_1.1M_merge.fam_correctIDs remote_UCDavis_GoogleDr:MAF/genotypingV1_pilot/
rclone -v --copy-links copy K9_1.1M_merge.fam remote_UCDavis_GoogleDr:MAF/genotypingV1_pilot/
