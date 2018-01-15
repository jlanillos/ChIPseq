#!/bin/bash

# This is the template sbatch script for ChIP-seq analysis using Bowtie2. B2491 jlanillos 2017
# The script is done to analyze three differen sample: 2 replicates and a control sample

# How to run in batch:  sbatch –A g2017020 –t 2:00:00 –p core –n 8 –o my_job_name.out sbatch_chipseq.sh

# Don't forget to add the appropriate modules! See THE MODULE SYSTEM section in the instructions.
module add bioinfo-tools
module add bowtie2
module add samtools
module add BEDTools
module add bcftools
module add MACS
# Location of the index ref genome: /proj/sllstore2017067/project6/Javi/RefGenome/GRCh38/all_chr_ref_GRCh38.fa

# cd to the directory where the data is:

cd /proj/sllstore2017067/project6/Javi/data


# 1. Bowtie2 indexing (must be done only once):

bowtie2-build /proj/sllstore2017067/project6/Javi/RefGenome/GRCh38/all_chr_ref_GRCh38.fa /proj/sllstore2017067/project6/Javi/data/bowtie2_indx/bowtie2_index_refgen_GRCh38

# The aligments (BAM and SAM) are stored in a new folder Javi/data/alignments
#Change directory to bowtie2 index files

cd /proj/sllstore2017067/project6/Javi/data/bowtie2_indx

# Bowtie2 alignment:
#(control, sample1 and sample2)

bowtie2 -x /proj/sllstore2017067/project6/Javi/data/bowtie2_indx/bowtie2_index_refgen_GRCh38 -U control.fq -S /proj/sllstore2017067/project6/Javi/data/alignments/control_bowtie2.sam
bowtie2 -x /proj/sllstore2017067/project6/Javi/data/bowtie2_indx/bowtie2_index_refgen_GRCh38 -U sample1.fq -S /proj/sllstore2017067/project6/Javi/data/alignments/sample1_bowtie2.sam
bowtie2 -x /proj/sllstore2017067/project6/Javi/data/bowtie2_indx/bowtie2_index_refgen_GRCh38 -U sample2.fq -S /proj/sllstore2017067/project6/Javi/data/alignments/sample2_bowtie2.sam

# Change dir to alignments (SAM)
cd /proj/sllstore2017067/project6/Javi/data/alignments

# SamTools

# Indexing SAM files. NECESSARY? Uncomment only if needed
#samtools faidx /proj/sllstore2017067/project6/Javi/RefGenome/GRCh38/all_chr_ref_GRCh38.fa

# Convert SAM to BAM:
samtools view -S -b control_bowtie2.sam > control_bowtie2.bam
samtools view -S -b sample1_bowtie2.sam > sample1_bowtie2.bam
samtools view -S -b sample2_bowtie2.sam > sample2_bowtie2.bam

# Remove unmapped reads and save the mapped reads only for the peak finders
samtools view -b -F 4 control_bowtie2.bam > control_mapped_bowtie2.bam
samtools view -b -F 4 sample1_bowtie2.bam > sample1_mapped_bowtie2.bam
samtools view -b -F 4 sample2_bowtie2.bam > sample2_mapped_bowtie2.bam

# Sort the BAM file in order to browse the alignment
samtools sort control_mapped_bowtie2.bam -o control_mapped_bowtie2.sorted.bam
samtools sort sample1_mapped_bowtie2.bam -o sample1_mapped_bowtie2.sorted.bam
samtools sort sample2_mapped_bowtie2.bam -o sample2_mapped_bowtie2.sorted.bam

# Remove potential duplicates (rmdup or PicardTools?)
samtools rmdup -s control_mapped_bowtie2.sorted.bam control_mapped_bowtie2.sorted.duplicates_filtered.bam
samtools rmdup -s sample1_mapped_bowtie2.sorted.bam sample1_mapped_bowtie2.sorted.duplicates_filtered.bam
samtools rmdup -s sample2_mapped_bowtie2.sorted.bam sample2_mapped_bowtie2.sorted.duplicates_filtered.bam

# Index BAM file to produce a BAI (uncomment if necessary)
#samtools index control_mapped_bowtie2.sorted.duplicates_filtered.bam
#samtools index sample1_mapped_bowtie2.sorted.duplicates_filtered.bam
#samtools index sample2_mapped_bowtie2.sorted.duplicates_filtered.bam

# Convert BAM to BED for peak finders (comment only if your MACS2 input is BAM)
# Make a directory to store the final BED files for the peak finder
mkdir BED
#bedtools bamtobed -i control_mapped_bowtie2.sorted.duplicates_filtered.bam > /proj/sllstore2017067/project6/Javi/data/alignments/BED/control.bed
#bedtools bamtobed -i sample1_mapped_bowtie2.sorted.duplicates_filtered.bam > /proj/sllstore2017067/project6/Javi/data/alignments/BED/sample1.bed
#bedtools bamtobed -i sample2_mapped_bowtie2.sorted.duplicates_filtered.bam > /proj/sllstore2017067/project6/Javi/data/alignments/BED/sample2.bed


### MACS2 ###
# cd to the BED/BAM files location
cd /proj/sllstore2017067/project6/Javi/data/alignments/

macs2 callpeak -t merged_mapped_bowtie2.sorted.duplicates_filt.bam -c control_mapped_bowtie2.sorted.duplicates_filtered.bam -f BAM -g hs -q 0.01 --outdir /proj/sllstore2017067/project6/Javi/data/MACS2 -n sample1_MACS2 -B --trackline --SPMR
macs2 callpeak -t sample2_mapped_bowtie2.sorted.duplicates_filtered.bam -c control_mapped_bowtie2.sorted.duplicates_filtered.bam -f BAM -g hs -q 0.01 --outdir /proj/sllstore2017067/project6/Javi/data/MACS2 -n sample2_MACS2 -B --trackline --SPMR

#DEFAULT: macs2 callpeak -t sample2.bed -c control.bed -f BED -g hs --outdir /proj/sllstore2017067/project6/Javi/data/MACS2 -n sample2_MACS2
## default: #macs2 callpeak -t sample2.bed -c control.bed -f BED -g hs --outdir /proj/sllstore2017067/project6/Javi/data/MACS2 -n sample2_MACS2


### HOMER motif finding ###
# Install and save HOMER in a location of interest.
# To run this step, please, prepare the peak files with the right format (HOMERpeaks.txt or .bed file as speficified in the HOMER website: http://homer.ucsd.edu/homer/ngs/peakMotifs.html

export PATH=$PATH:/path/to/HOMER/bin
cd /path/to/HOMER/bin
# Create output directories
mkdir /proj/sllstore2017067/project6/Javi/data/HOMER/output_sample1
mkdir /proj/sllstore2017067/project6/Javi/data/HOMER/output_sample2

# Run Motif discovery "findMotifsGenome.pl"
# The input files called HOMER_inputBED.sample1 are specifically modified from MACS2 to be compatible with HOMER
perl /path/to/HOMER/bin/findMotifsGenome.pl /proj/sllstore2017067/project6/Javi/data/HOMER/HOMER_inputBED.sample1 hg38 /proj/sllstore2017067/project6/Javi/data/HOMER/output_sample1
perl /path/to/HOMER/bin/findMotifsGenome.pl /proj/sllstore2017067/project6/Javi/data/HOMER/HOMER_inputBED.sample2 hg38 /proj/sllstore2017067/project6/Javi/data/HOMER/output_sample2
#############
