#!/bin/bash

# This is the template sbatch script for ChIP-seq analysis. B2491 jlanillos 2017

# sbatch –A g2017020 –t 2:00:00 –p core –n 8 –o my_job_name.out sbatch_template.sh
# Location of the Reference Genome: /proj/sllstore2017067/project6/Andre/project/ref_genome/GRCh38_r90.all.fa

# Don't forget to add the appropriate modules! See THE MODULE SYSTEM section in the instructions.
module add bioinfo-tools
module add bwa/0.7.15
module add samtools
module add BEDTools


# cd to the directory where the data is:
cd /proj/sllstore2017067/project6/Javi/data

# Then list the commands that you would like to run. 
###cd /proj/sllstore2017067/project6/Andre/project/ref_genome/
#bwa index -a bwtsw /proj/g2017020/BB2491/LAB2/all_chr.hg19.fa
#bwa aln /proj/g2017020/BB2491/LAB2/all_chr.hg19.fa control.fq > control.filt.bwa_aln.default.sai
bwa samse /proj/g2017020/BB2491/LAB2/all_chr.hg19.fa control.filt.bwa_aln.default.sai control.filt.fastq > control.filt.bwa_aln.default.sam


# SamTools Indexing

# Creates a .fai. Available at: /proj/g2017020/BB2491/LAB2/all_chr.hg19.fa.fai
samtools faidx /proj/g2017020/BB2491/LAB2/all_chr.hg19.fa > all_chr.hg19.fa.fai

# Convert SAM to BAM:
samtools import /proj/g2017020/BB2491/LAB2//all_chr.hg19.fa.fai control.filt.bwa_aln.default.sam control.filt.bwa_aln.default.bam

# Sort the BAM file in order to browse the alignment
samtools sort -@ 8 -m 6G control.filt.bwa_aln.default.bam -o control.filt.bwa_aln.default.sorted.bam
# Index the BAM file for fast access
samtools index control.filt.bwa_aln.default.sorted.bam
# View using samtools view
samtools view control.filt.bwa_aln.default.sorted.bam

# In case you want to direct your output towards certain regions, modify and uncomment the next line:
samtools view ERR001014.filt.bwa_aln.default.sorted.bam "chr22_hg19:18000000-18080000" > chr22_part.ERR001014.filt.bwa_aln.default.sam


# BEDtools to convert the BAM aligned file to BED file:
# To restrict by Chromosome
#bamToBed -i ERR001014.filt.bwa_aln.default.sorted.bam | awk '/^chr22_hg19/ {if ($2>=18000000 && $2<=18080000) print "chr22", substr($0,11)}' > chr22_part.ERR001014.filt.bwa_aln.default.sorted.bed
bamToBed -i control.filt.bwa_aln.default.sorted.bam > control.bed

 
