#!/bin/bash -l

#SBATCH -J SNP_Marker
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
SAMPLE_NAME=UD278 # Make sure your folder name is the same as the sample name for ease of use
FA_RUN=FA1 # Run type that your are doing

#Assign folder for sample
SSM_FOLDER=/home/jbk708/emu_strains/UD278/ 

# WT .fastq file names (FWD and REV)
SSM_FWD=UD278_FWD.fastq # FWD
SSM_REV=UD278_REV.fastq # REV 

### Sets the run type variables ###
###Don't mess with this unless you're sure what you're doing###
if [ "$FA_RUN" == "FA1" ]; then
	MIN_VAR_FREQ=0.01
	MIN_HOM_FREQ=1.0;
fi

if [ "$FA_RUN" == "FA2" ]; then
	echo "FA2"
	MIN_VAR_FREQ=0.1
	MIN_HOM_FREQ=1.0;
fi

if [ "$FA_RUN" == "FA3" ]; then
	echo "FA3"
	MIN_VAR_FREQ_=0.1
	MIN_HOM_FREQ=0.9;
fi

if [ "$FA_RUN" == "FA4" ]; then
	echo "FA4"
	MIN_VAR_FREQ=0.01
	MIN_HOM_FREQ=0.9;
fi

if [ "$FA_RUN" == "FA5" ]; then
	echo "FA5"
	MIN_VAR_FREQ=0.01
	MIN_HOM_FREQ=0.8;
fi

#Step 1: Trimmomatic

java -jar /share/apps/Trimmomatic-0.36/trimmomatic.jar PE -threads 32 \
$SSM_FOLDER/$SSM_FWD \
$SSM_FOLDER/$SSM_REV \
$SSM_FOLDER/${SAMPLE_NAME}_FwdPaired.fastq \
$SSM_FOLDER/${SAMPLE_NAME}_FwdSingle.fastq \
$SSM_FOLDER/${SAMPLE_NAME}_RevPaired.fastq \
$SSM_FOLDER/${SAMPLE_NAME}_RevSingle.fastq \
ILLUMINACLIP:adapter.txt:2:30:10 HEADCROP:8 SLIDINGWINDOW:4:25 MINLEN:75

#Step 2 Map BWA

bwa mem -A 5 -B 10 -O 60 -t 32 \
$REF_GENOME \
$SSM_FOLDER/${SAMPLE_NAME}_FwdPaired.fastq \
$SSM_FOLDER/${SAMPLE_NAME}_RevPaired.fastq \
> $SSM_FOLDER/${SAMPLE_NAME}_mapped_reads.sam

#Step 3 Samtools View 

samtools view -o $SSM_FOLDER/${SAMPLE_NAME}_filterMapReads.bam \
-h -b -F 0x30c \
$SSM_FOLDER/${SAMPLE_NAME}_mapped_reads.sam

java -Xmx32G -jar /home/mitochi/src/picard-tools-2.1.1/picard.jar \
AddOrReplaceReadGroups \
I=$SSM_FOLDER/${SAMPLE_NAME}_filterMapReads.bam \
O=$SSM_FOLDER/${SAMPLE_NAME}_filterMapReads_rg.bam \
RGSM=rgSM RGLB=rgLB RGPL=illumina RGPU=rgPU \

#Step 4 samtools sort and index

samtools sort -l 6 -@ 32 \
-o $SSM_FOLDER/${SAMPLE_NAME}_SortFilRG.bam \
-T $SSM_FOLDER/${SAMPLE_NAME}_sorted \
$SSM_FOLDER/${SAMPLE_NAME}_filterMapReads_rg.bam
samtools index $SSM_FOLDER/${SAMPLE_NAME}_SortFilRG.bam

# #Step 5 GATK Realigner Target Creator 

java -Xmx32G -jar /home/mitochi/src/GenomeAnalysisTK/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R $REF_GENOME \
-I $SSM_FOLDER/${SAMPLE_NAME}_SortFilRG.bam \
-o IndelRealigner.intervals

# #Step 6 GATK Indel Realigner 

java -Xmx32G -jar /home/mitochi/src/GenomeAnalysisTK/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R $REF_GENOME \
-I $SSM_FOLDER/${SAMPLE_NAME}_SortFilRG.bam \
-targetIntervals IndelRealigner.intervals \
-o $SSM_FOLDER/${SAMPLE_NAME}_RealignedReads.bam

#Step 7 Samtools pileup (TAKES FORRREVER)

samtools mpileup -C 50 \
-f $REF_GENOME \
-o $SSM_FOLDER/${SAMPLE_NAME}_pileup \
$SSM_FOLDER/${SAMPLE_NAME}_RealignedReads.bam

# Step 8 VarScan SNP Pileup

java -jar VarScan.jar mpileup2snp \
$SSM_FOLDER/${SAMPLE_NAME}_pileup \
--output-vcf 1 \
--min-coverage 1 \
--min-reads2 1 \
--min-avg-qual 24 \
--min-var-freq $MIN_VAR_FREQ \
--p-value 0.99 \
--variants 1 \
--min-freq-for-hom $MIN_HOM_FREQ \
> $SSM_FOLDER/${SAMPLE_NAME}_snp_pileup_${FA_RUN}.vcf 

# Step 9 VarScan Indel Pileup

java -jar VarScan.jar mpileup2indel \
$SSM_FOLDER/${SAMPLE_NAME}_pileup \
--min-coverage 1 \
--min-reads2 1 \
--min-avg-qual 24 \
--min-var-freq $MIN_VAR_FREQ \
--p-value 0.99 \
--output-vcf 1 \
--variants 1 \
--min-freq-for-hom $MIN_HOM_FREQ \
> $SSM_FOLDER/${SAMPLE_NAME}_indel_pileup_${FA_RUN}.vcf

# Step 10 VCF Merge (Combines INDEL and SNP)

/home/jbk708/apps/vcflib/bin/vcfcombine \
$SSM_FOLDER/${SAMPLE_NAME}_snp_pileup_${FA_RUN}.vcf  \
$SSM_FOLDER/${SAMPLE_NAME}_indel_pileup_${FA_RUN}.vcf \
>$SSM_FOLDER/${SAMPLE_NAME}_combo_pileup_${FA_RUN}.vcf

### RUN snpEff on Galaxy!!! ###














