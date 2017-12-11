#!/bin/bash

# This is the template sbatch script for ChIP-seq analysis using Bowtie2. B2491 jlanillos 2017

# sbatch –A g2017020 –t 2:00:00 –p core –n 8 –o my_job_name.out sbatch_template.sh

# Don't forget to add the appropriate modules! See THE MODULE SYSTEM section in the instructions.
module add bioinfo-tools
module add bowtie2
module add samtools
module add BEDTools
module add bcftools
# Location of the index ref genome: /proj/sllstore2017067/project6/Javi/RefGenome/GRCh38/all_chr_ref_GRCh38.fa

# cd to the directory where the data is:
cd /proj/sllstore2017067/project6/Javi/data


# Bowtie2 indexing:

###bowtie2-build /proj/sllstore2017067/project6/Javi/RefGenome/GRCh38/all_chr_ref_GRCh38.fa /proj/sllstore2017067/project6/Javi/data/bowtie2_indx/bowtie2_index_refgen_GRCh38
 # The aligments (BAM and SAM)  are stored in a new folder Javi/data/alignments
#Change directory to bowtie2 index files
cd /proj/sllstore2017067/project6/Javi/data/bowtie2_indx
# Bowtie2 alignment:
#(control, sample1 and sample2)
#bowtie2 -x /proj/sllstore2017067/project6/Javi/data/bowtie2_indx/bowtie2_index_refgen_GRCh38 -U control.fq -S /proj/sllstore2017067/project6/Javi/data/alignments/control_bowtie2.sam
#bowtie2 -x /proj/sllstore2017067/project6/Javi/data/bowtie2_indx/bowtie2_index_refgen_GRCh38 -U sample1.fq -S /proj/sllstore2017067/project6/Javi/data/alignments/sample1_bowtie2.sam
#bowtie2 -x /proj/sllstore2017067/project6/Javi/data/bowtie2_indx/bowtie2_index_refgen_GRCh38 -U sample2.fq -S /proj/sllstore2017067/project6/Javi/data/alignments/sample2_bowtie2.sam

# Change dir to alignments (SAM)
cd /proj/sllstore2017067/project6/Javi/data/alignments

# SamTools

# Indexing SAM files. NECESSARY????
##samtools faidx /proj/sllstore2017067/project6/Javi/RefGenome/GRCh38/all_chr_ref_GRCh38.fa

# Convert SAM to BAM:
#samtools view -S -b control_bowtie2.sam > control_bowtie2.bam
#samtools view -S -b sample1_bowtie2.sam > sample1_bowtie2.bam
#samtools view -S -b sample2_bowtie2.sam > sample2_bowtie2.bam


# Sort the BAM file in order to browse the alignment
##samtools sort control_bowtie2.bam -o control_bowtie2.sorted.bam
#samtools sort sample1_bowtie2.bam -o sample1_bowtie2.sorted.bam
#samtools sort sample2_bowtie2.bam -o sample2_bowtie2.sorted.bam
# Remove potential duplicates (rmdup or PicardTools?)
#samtools rmdup -s control_bowtie2.sorted.bam control_bowtie2.duplicates_filtered.bam
#samtools rmdup -s sample1_bowtie2.sorted.bam sample1_bowtie2.duplicates_filtered.bam
#samtools rmdup -s sample2_bowtie2.sorted.bam sample2_bowtie2.duplicates_filtered.bam

# Index BAM file to produce a BAI
#samtools index control_bowtie2.duplicates_filtered.bam
#samtools index sample1_bowtie2.duplicates_filtered.bam
#samtools index sample2_bowtie2.duplicates_filtered.bam

# Convert BAM to BED for peak finders
#bedtools bamtobed -i control_bowtie2.duplicates_filtered.bam > control.bed
#bedtools bamtobed -i sample1_bowtie2.duplicates_filtered.bam > sample1.bed
#bedtools bamtobed -i sample2_bowtie2.duplicates_filtered.bam > sample2.bed

# Generate VCF variants:
samtools mpileup -uf control_bowtie2.duplicates_filtered.bam | bcftools view -Ov - > control_bowtie2.bcf
samtools mpileup -uf sample1_bowtie2.duplicates_filtered.bam | bcftools view -Ov - > sample1_bowtie2.bcf
samtools mpileup -uf sample2_bowtie2.duplicates_filtered.bam | bcftools view -Ov - > sample2_bowtie2.bcf

# To get the unmapped reads from a bam file use : (WHICH INPUT SHOULD I USE??)
samtools view -f 4 control_bowtie2.sorted.bam > control_bowtie2_unmapped.sam
samtools view -f 4 sample1_bowtie2.sorted.bam > sample1_bowtie2_unmapped.sam
samtools view -f 4 sample2_bowtie2.sorted.bam > sample2_bowtie2_unmapped.sam




# To get only the mapped reads use the parameter 'F', which works like -v of grep and skips the alignments for a specific flag.
##samtools view -b -F 4 control_bowtie2.sorted.bam > control_bowtie2_mapped.bam


#??? Index the BAM file for fast access
##samtools index control.filt.bwa_aln.default.sorted.bam

#??? View using samtools view
##samtools view control.filt.bwa_aln.default.sorted.bam

#??? In case you want to direct your output towards certain regions, modify and uncomment the next line:
#samtools view ERR001014.filt.bwa_aln.default.sorted.bam "chr22_hg19:18000000-18080000" > chr22_part.ERR001014.filt.bwa_aln.default.sam


# BEDtools to convert the BAM aligned file to BED file:
# To restrict by Chromosome
#bamToBed -i ERR001014.filt.bwa_aln.default.sorted.bam | awk '/^chr22_hg19/ {if ($2>=18000000 && $2<=18080000) print "chr22", substr($0,11)}' > chr22_part.ERR001014.filt.bwa_aln.default.sorted.bed
