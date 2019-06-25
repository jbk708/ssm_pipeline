#!/bin/bash -l

#SBATCH -J SSM
#SBATCH -c 32
#SBATCH -t 24:00:00
#SBATCH -p high

# Name: SSM.sh
# Created by: Jonathan Blake Kirkland
# Creation date: November 30, 2018
# Updated Last: June 25, 2019
# This is the pipeline for the Sibling Subtraction Method

## Load the required software module
### MAKE SURE THAT THE BWA MODULE IS IN THE SAME LOCATION AS THE SCRIPT ###
### OR YOUR'E GONNA HAVE A BAD TIME ###

module load java/1.8
module load trimmomatic/0.36
module load samtools
module load bwa/0.7.9a

#Assign Reference Genome
#NOTE: IF YOU HAVE NOT NOT INDEXED YOUR REFERENCED GENOME BEFORE, 
#PLEASE RUN THE COMMENTED "bwa index" SECTION BELOW OR NOTHING WILL WORK!!
#Also make sure there are not special characters like brackets in your file name
#or the crick gods will frown upon you and deny your request

# Reference Geneome Location #
REF_GENOME=/home/jbk708/bwa_index_ce/WS220.64_chr.fa # Change this to your reference gneome location

#REF GENEOME INDEX

# bwa index -p $REF_GENOME -a bwtsw $REF_GENOME
# java -jar /home/mitochi/src/picard-tools-2.1.1/picard.jar CreateSequenceDictionary \
# R=$REF_GENOME \
# O=WS220.64_chr.dict
# samtools faidx $REF_GENOME

#Sample Name Folder
SAMPLE_NAME=UD557 # Make sure your folder name is the same as the sample name for ease of use
FA_RUN=FA1 # Run type that your are doing

#Assign folder for sample
SSM_FOLDER_WT=/home/jbk708/UD557/WT # Wild type strain folder
SSM_FOLDER_MUT=/home/jbk708/UD557/Mut # mutant strain folder

# WT .fastq file names (FWD and REV)
SSM_WT_FWD=1278S_FwdPaired.fastq # FWD
SSM_WT_REV=1278S_RevPaired.fastq # REV 

#Mut .fastq file names (FWD and REV)
SSM_MUT_FWD=1278P_FwdPaired.fastq # FWD
SSM_MUT_REV=1278P_RevPaired.fastq # REV

#Results Folder
SSM_RESULTS_FOLDER=/home/jbk708/${SAMPLE_NAME} # Change jbk708 to your username


### Sets the run type variables ###
###Don't mess with this unless you're sure what you're doing###
if [ "$FA_RUN" == "FA1" ]; then
	MIN_VAR_FREQ_WT=0.01
	MIN_VAR_FREQ_MUT=1.0
	MIN_HOM_FREQ=1.0;
fi

if [ "$FA_RUN" == "FA2" ]; then
	echo "FA2"
	MIN_VAR_FREQ_WT=0.1
	MIN_VAR_FREQ_MUT=0.9
	MIN_HOM_FREQ=1.0;
fi

if [ "$FA_RUN" == "FA3" ]; then
	echo "FA3"
	MIN_VAR_FREQ_WT=0.1
	MIN_VAR_FREQ_MUT=0.9
	MIN_HOM_FREQ=0.9;
fi

if [ "$FA_RUN" == "FA4" ]; then
	echo "FA4"
	MIN_VAR_FREQ_WT=0.01
	MIN_VAR_FREQ_MUT=0.9
	MIN_HOM_FREQ=0.9;
fi

if [ "$FA_RUN" == "FA5" ]; then
	echo "FA5"
	MIN_VAR_FREQ_WT=0.01
	MIN_VAR_FREQ_MUT=0.8
	MIN_HOM_FREQ=0.8;
fi

#Step 1: Trimmomatic WT 

java -jar /share/apps/Trimmomatic-0.36/trimmomatic.jar PE -threads 32 \
$SSM_FOLDER_WT/$SSM_WT_FWD \
$SSM_FOLDER_WT/$SSM_WT_REV \
$SSM_FOLDER_WT/${SAMPLE_NAME}_FwdPaired.fastq \
$SSM_FOLDER_WT/${SAMPLE_NAME}_FwdSingle.fastq \
$SSM_FOLDER_WT/${SAMPLE_NAME}_RevPaired.fastq \
$SSM_FOLDER_WT/${SAMPLE_NAME}_RevSingle.fastq \
ILLUMINACLIP:adapter.txt:2:30:10 HEADCROP:8 SLIDINGWINDOW:4:25 MINLEN:75

# #Step 1.1: Trimmomatic Mut

java -jar /share/apps/Trimmomatic-0.36/trimmomatic.jar PE -threads 32 \
$SSM_FOLDER_MUT/$SSM_MUT_FWD \
$SSM_FOLDER_MUT/$SSM_MUT_REV \
$SSM_FOLDER_MUT/${SAMPLE_NAME}_Mut_FwdPaired.fastq \
$SSM_FOLDER_MUT/${SAMPLE_NAME}_Mut_FwdSingle.fastq \
$SSM_FOLDER_MUT/${SAMPLE_NAME}_Mut_RevPaired.fastq \
$SSM_FOLDER_MUT/${SAMPLE_NAME}_Mut_RevSingle.fastq \
ILLUMINACLIP:adapter.txt:2:30:10 HEADCROP:8 SLIDINGWINDOW:4:25 MINLEN:75

#Step 2 Map BWA WT

bwa mem -A 5 -B 10 -O 60 -t 32 \
$REF_GENOME \
$SSM_FOLDER_WT/${SAMPLE_NAME}_FwdPaired.fastq \
$SSM_FOLDER_WT/${SAMPLE_NAME}_RevPaired.fastq \
> $SSM_FOLDER_WT/${SAMPLE_NAME}_mapped_reads.sam

#Step 2.1 Map BWA Mut

bwa mem -A 5 -B 10 -O 60 -t 32 \
$REF_GENOME \
$SSM_FOLDER_MUT/${SAMPLE_NAME}_Mut_FwdPaired.fastq \
$SSM_FOLDER_MUT/${SAMPLE_NAME}_Mut_RevPaired.fastq \
> $SSM_FOLDER_MUT/${SAMPLE_NAME}_Mut_mapped_reads.sam

#Step 3 Samtools View WT

samtools view -o $SSM_FOLDER_WT/${SAMPLE_NAME}_filterMapReads.bam \
-h -b -F 0x30c \
$SSM_FOLDER_WT/${SAMPLE_NAME}_mapped_reads.sam

#Step 3.1 Samtools View Mut

samtools view -o $SSM_FOLDER_MUT/${SAMPLE_NAME}_Mut_filterMapReads.bam \
-h -b -F 0x30c \
$SSM_FOLDER_MUT/${SAMPLE_NAME}_Mut_mapped_reads.sam

#Step 4 Picard Tools: Add or Replace Read Groups for WT

java -Xmx32G -jar /home/mitochi/src/picard-tools-2.1.1/picard.jar \
AddOrReplaceReadGroups \
I=$SSM_FOLDER_WT/${SAMPLE_NAME}_filterMapReads.bam \
O=$SSM_FOLDER_WT/${SAMPLE_NAME}_filterMapReads_rg.bam \
RGSM=rgSM RGLB=rgLB RGPL=illumina RGPU=rgPU \

#Step 4.1 Picard Tools: Add or Replace Read Groups for Mut

java -Xmx32G -jar /home/mitochi/src/picard-tools-2.1.1/picard.jar \
AddOrReplaceReadGroups \
I=$SSM_FOLDER_MUT/${SAMPLE_NAME}_Mut_filterMapReads.bam \
O=$SSM_FOLDER_MUT/${SAMPLE_NAME}_Mut_filterMapReads_rg.bam \
RGSM=rgSM RGLB=rgLB RGPL=illumina RGPU=rgPU \

#Step 5 samtools sort and index for WT

samtools sort -l 6 -@ 32 \
-o $SSM_FOLDER_WT/${SAMPLE_NAME}_SortFilRG.bam \
-T $SSM_FOLDER_WT/${SAMPLE_NAME}_sorted \
$SSM_FOLDER_WT/${SAMPLE_NAME}_filterMapReads_rg.bam
samtools index $SSM_FOLDER_WT/${SAMPLE_NAME}_SortFilRG.bam

#Step 5.1 samtools sort and index for Mut

samtools sort -l 6 -@ 32 \
-o $SSM_FOLDER_MUT/${SAMPLE_NAME}_Mut_SortFilRG.bam \
-T $SSM_FOLDER_MUT/${SAMPLE_NAME}_Mut_sorted \
$SSM_FOLDER_MUT/${SAMPLE_NAME}_Mut_filterMapReads_rg.bam
samtools index $SSM_FOLDER_MUT/${SAMPLE_NAME}_Mut_SortFilRG.bam

# #Step 6 GATK Realigner Target Creator for WT

java -Xmx32G -jar /home/mitochi/src/GenomeAnalysisTK/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R $REF_GENOME \
-I $SSM_FOLDER_WT/${SAMPLE_NAME}_SortFilRG.bam \
-o IndelRealigner.intervals


# #Step 6.1 GATK Indel Realigner for WT

java -Xmx32G -jar /home/mitochi/src/GenomeAnalysisTK/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R $REF_GENOME \
-I $SSM_FOLDER_WT/${SAMPLE_NAME}_SortFilRG.bam \
-targetIntervals IndelRealigner.intervals \
-o $SSM_FOLDER_WT/${SAMPLE_NAME}_RealignedReads.bam


# #Step 7 GATK Realigner Target Creator for Mut

java -Xmx32G -jar /home/mitochi/src/GenomeAnalysisTK/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R $REF_GENOME \
-I $SSM_FOLDER_MUT/${SAMPLE_NAME}_Mut_SortFilRG.bam \
-o IndelRealigner.intervals

# #Step 7.1 GATK Indel Realigner for Mut

java -Xmx32G -jar /home/mitochi/src/GenomeAnalysisTK/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R $REF_GENOME \
-I $SSM_FOLDER_MUT/${SAMPLE_NAME}_Mut_SortFilRG.bam \
-targetIntervals IndelRealigner.intervals \
-o $SSM_FOLDER_MUT/${SAMPLE_NAME}_Mut_RealignedReads.bam

#Step 8 Samtools pileup for WT (TAKES FORRREVER)

samtools mpileup -C 50 \
-f $REF_GENOME \
-o $SSM_FOLDER_WT/${SAMPLE_NAME}_pileup \
$SSM_FOLDER_WT/${SAMPLE_NAME}_RealignedReads.bam

#Step 8.1 Samtools pileup for Mut (TAKES FORRREVER)

samtools mpileup -C 50 \
-f $REF_GENOME \
-o $SSM_FOLDER_MUT/${SAMPLE_NAME}_Mut_pileup_noPicard \
$SSM_FOLDER_MUT/${SAMPLE_NAME}_Mut_RealignedReads.bam

# Step 9 VarScan SNP Pileup for WT

java -jar VarScan.jar mpileup2snp \
$SSM_FOLDER_WT/${SAMPLE_NAME}_pileup \
--output-vcf 1 \
--min-coverage 1 \
--min-reads2 1 \
--min-avg-qual 24 \
--min-var-freq MIN_VAR_FREQ_WT \
--p-value 0.99 \
--variants 1 \
--min-freq-for-hom $MIN_HOM_FREQ \
> $SSM_FOLDER_WT/${SAMPLE_NAME}_snp_pileup_${FA_RUN}.vcf 

# Step 9.1 VarScan Indel Pileup for WT

java -jar VarScan.jar mpileup2indel \
$SSM_FOLDER_WT/${SAMPLE_NAME}_pileup \
--min-coverage 1 \
--min-reads2 1 \
--min-avg-qual 24 \
--min-var-freq MIN_VAR_FREQ_WT \
--p-value 0.99 \
--output-vcf 1 \
--variants 1 \
--min-freq-for-hom $MIN_HOM_FREQ \
> $SSM_FOLDER_WT/${SAMPLE_NAME}_indel_pileup_${FA_RUN}.vcf

#Step 10 VarScan SNP Pileup for Mut

java -jar VarScan.jar mpileup2snp \
$SSM_FOLDER_MUT/${SAMPLE_NAME}_Mut_pileup_noPicard \
--min-coverage 8 \
--min-reads2 1 \
--min-avg-qual 24 \
--min-var-freq $MIN_VAR_FREQ_MUT \
--p-value 0.99 \
--output-vcf 1 \
--variants 1 \
--min-freq-for-hom $MIN_HOM_FREQ \
> $SSM_FOLDER_MUT/${SAMPLE_NAME}_Mut_snp_pileup_${FA_RUN}.vcf

# Step 10.1 VarScan Indel Pileup for Mut

java -jar VarScan.jar mpileup2indel \
$SSM_FOLDER_MUT/${SAMPLE_NAME}_Mut_pileup_noPicard \
--min-coverage 8 \
--min-reads2 1 \
--min-avg-qual 24 \
--min-var-freq $MIN_VAR_FREQ_MUT \
--p-value 0.99 \
--output-vcf 1 \
--variants 1 \
--min-freq-for-hom $MIN_HOM_FREQ \
> $SSM_FOLDER_MUT/${SAMPLE_NAME}_Mut_indel_pileup_${FA_RUN}.vcf

# Step 11 VCF Merge WT (Combines INDEL and SNP)

/home/jbk708/apps/vcflib/bin/vcfcombine \
$SSM_FOLDER_WT/${SAMPLE_NAME}_snp_pileup_${FA_RUN}.vcf  \
$SSM_FOLDER_WT/${SAMPLE_NAME}_indel_pileup_${FA_RUN}.vcf \
>$SSM_FOLDER_WT/${SAMPLE_NAME}_combo_pileup_${FA_RUN}.vcf

#Step 11.1 VCF Merge Mut (Combines INDEL and SNP)

/home/jbk708/apps/vcflib/bin/vcfcombine \
$SSM_FOLDER_MUT/${SAMPLE_NAME}_Mut_snp_pileup_${FA_RUN}.vcf  \
$SSM_FOLDER_MUT/${SAMPLE_NAME}_Mut_indel_pileup_${FA_RUN}.vcf \
>$SSM_FOLDER_MUT/${SAMPLE_NAME}_Mut_combo_pileup_${FA_RUN}.vcf

# # Step 12 VCF-VCF Intersect (Performs the Subraction)

/home/jbk708/apps/vcflib/bin/vcfintersect \
--intersect-vcf $SSM_FOLDER_WT/${SAMPLE_NAME}_combo_pileup_${FA_RUN}.vcf \
--invert \
-r /home/jbk708/bwa_index_ce/WS220.64_chr.fa \
--window-size 30 \
$SSM_FOLDER_MUT/${SAMPLE_NAME}_Mut_combo_pileup_${FA_RUN}.vcf \
>$SSM_RESULTS_FOLDER/${SAMPLE_NAME}_subtracted_${FA_RUN}.vcf

### RUN snpEff on Galaxy!!! ###



