# ChIPseq
Project on HT data analysis of ChIP-seq data (ERC-alpha and beta)

Content of this repo:

- SampleID.txt: samples were renamed to easy handling. This doc contains the original names
- sbatch_chipseq.sh: script to submit to the CLuster containing the pipeline (tools used) of this project:
    Bowtie2 (indexing and alignment) --> (samtools) --> MACS2 peak finder --> HOMER motifs discovery (if you want, you can add AnnotatePeaks.pl script from HOMER)
- LOGBOOK folder contains all the submitted jobs (.log .out files) during the analysis
- FastQC: fastQC results of my specific samples before trimming
    
Data and results is not public yet in this project.
For more info about the tools used and the source of code, please READ THE WIKI of this project.
