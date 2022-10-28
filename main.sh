mkdir -p MAF && cd MAF
work_dir=$(pwd)
## Download from Google cloud
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

####
## Pedigree ped files were obtained from Slack. Therefore, I had to download locally then upload to the server
mkdir -p $work_dir/pedigree && cd pedigree
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
cat sample_ped_infoA_cohort.txt | cut -f2,3 | sort > a2.txt
cat sample_ped_infoB_cohort.txt | cut -f2,3 | sort > b2.txt
diff a2.txt b2.txt
rm a2.txt b2.txt

cat sample_ped_infoA_cohort.txt | cut -f3 | sort > a3.txt
cat sample_ped_infoB_cohort.txt | cut -f3 | sort > b3.txt
diff a3.txt b3.txt
rm a3.txt b3.txt

## create map files for CEL file names and sample IDs
echo "Sample Filename|Alternate Sample Name" | tr '|' '\t' > infoA_names.txt
tail -n+2 sample_ped_infoA_cohort.txt | awk 'BEGIN{FS=OFS="\t"}{print $1,"S"$3}' >> infoA_names.txt

echo "Sample Filename|Alternate Sample Name" | tr '|' '\t' > infoB_names.txt
tail -n+2 sample_ped_infoB_cohort.txt | awk 'BEGIN{FS=OFS="\t"}{print $1,"S"$3}' >> infoB_names.txt

## find sample name duplicates
tail -n+2 infoA_names.txt | cut -f2 | cut -d"_" -f1 | sort | uniq -c | awk '{if($1>1)print $2}' | while read f;do echo $(grep "$f" infoA_names.txt | cut -f2 | tr '\n' '\t');done > sample_dup.txt

## compare gender across info files
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$3]=$6;next}{print $3,a[$3],$6}' sample_ped_infoA_cohort.txt sample_ped_infoB_cohort.txt | sed 's/Sex/sexA/' | sed 's/Sex/sexB/' > gender_compare.txt
tail -n+2 gender_compare.txt | awk 'BEGIN{FS=OFS="\t"}{if($2!=$3)print}'
tail -n+2 gender_compare.txt | cut -f2 | sort | uniq -c  ##    1575 1(male) // 1587 2(female) // 6 -9(unknown)

########
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
cel_PATH=$work_dir/"CEL_Files/Array_A"
set_Analysis=$work_dir/"lib/setA/setA_Analysis/"
set_output=$work_dir/"output/setA"
qc_xml="$set_Analysis"/Axiom_K9HDSNPA.r1.apt-geno-qc.AxiomQC1.xml
genotype1_xml="$set_Analysis"/Axiom_K9HDSNPA_96orMore_Step1.r1.apt-genotype-axiom.AxiomGT1.apt2.xml
genotype2_xml="$set_Analysis"/Axiom_K9HDSNPA_96orMore_Step2.r1.apt-genotype-axiom.AxiomGT1.apt2.xml
specialSNPs="$set_Analysis"/Axiom_K9HDSNPA.r1.specialSNPs
ps2snp="$set_Analysis"/Axiom_K9HDSNPA.r1.ps2snp_map.ps
annot_db="$set_Analysis"/../Axiom_K9HDSNPA.r1_3.20160810.annot.db
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
awk 'BEGIN{FS=" "}FNR==NR{a[$2]=$1;next}{if(a[$2])print $2,$1/a[$2] * 100}' "$set_output"/samples_on_plate "$set_output"/samples_pass_plate  > "$set_output"/plate_pass_rate ## Array_A (131, 137, 438) have the failing samples but all are above 97.5 // Array_B (9 plates have lost samples)
grep -v "^#"  "$set_output"/step1/AxiomGT1.report.txt | tail -n+2 | cut -f1,3 | cut -d"-" -f4 | awk -F"[_\t]" 'BEGIN{CONVFMT="%.3g"}{if($3>95){sumQC[$1]+=$3;count[$1]++}}END{for(i in sumQC)print i,sumQC[i]/count[i]}' > "$set_output"/plate_aveQC ## Array_A: 129 (ave=98.3965). This is the only plate failing 2 samnples // Array_B: 679 (ave=95.8125), 667 (ave=97.675), 678 (ave=97.775)

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
$apt/apt-format-result \
      --calls-file "$set_output"/step2/AxiomGT1.calls.txt \
      --snp-list-file "$set_output"/SNPolisher/Recommended.ps \
      --annotation-file lib/setA/Axiom_K9HDSNPA.r1_3.20160810.annot.db \
      --snp-identifier-column Affy_SNP_ID \
      --export-chr-shortname true \
      --export-vcf-file "$set_output"/export/AxiomGT1.vcf \
      --log-file "$set_output"/export/AxiomGT1_to_vcf.log
      
      
########
########
## Example for seting up an interactive session on Farm
srun -p bmm -t 03-00:00:00 -c 16 -n 1 -N 1 --mem=128000 --job-name "jobName" --pty bash
apt=$work_dir/"apt_2.11.6_linux_64_x86_binaries/bin/"
########
## Example for ps-extract then export
mkdir -p example/subSample
head -n10 "$set_output"/SNPolisher/Recommended.ps > example/subSample/Recommended.ps
echo "sample_list a550771-4434579-040223-600_A02.CEL" | tr ' ' '\n' > example/subSample/samples.txt
echo "Sample Filename|Alternate Sample Name" | tr '|' '\t' > example/subSample/samples_att.txt
echo "a550771-4434579-040223-600_A02.CEL|sample_1" | tr '|' '\t' >> example/subSample/samples_att.txt

$apt/ps-extract \
     --call-file "$set_output"/step2/AxiomGT1.calls.txt \
     --pid-file ./example/subSample/Recommended.ps \
     --sample-list-file ./example/subSample/samples.txt \
     --output-dir ./example/subSample_output

$apt/apt-format-result \
      --calls-file ./example/subSample_output/extract_calls.txt \
      --snp-list-file ./example/subSample/Recommended.ps \
      --annotation-file lib/setA/Axiom_K9HDSNPA.r1_3.20160810.annot.db \
      --snp-identifier-column Affy_SNP_ID \
      --export-chr-shortname true \
      --export-vcf-file ./example/subSample_output/export/AxiomGT1.vcf \
      --log-file ./example/subSample_output/export/AxiomGT1_to_vcf.log \
      --export-alternate-sample-names true \
      --sample-attributes-file ./example/subSample/samples_att.txt

