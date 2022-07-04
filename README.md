# GATK-best-practice vs Optimized Pipeline for Improved Detection of Germline Exomes
# Introduction
The next-generation sequencing (NGS) technology represents a significant advance in genomics
and medical diagnosis. Nevertheless, the time it takes to perform sequencing, data analysis,
and variant interpretation is a bottleneck in using the next-generation sequencing in precision
medicine. For accurate and efficient performance in a clinical diagnostic lab practice, a consistent data analysis pipeline is necessary to avoid false variant calls and achieve optimum
accuracy. The purpose of this study is to compare the performance of two NGS data analysis
pipeline compartments, including short-read mapping (BWA-MEM, and BWA-MEM2), and
variant calling (GATK-HaplotypeCaller and DRAGEN-GATK) on Whole Exome Sequencing
(WES) data. Here we provide two WES analysis piplines start from fatsq till vcf files. Then , you can  compare the performance of two WES data analysis
pipeline compartments using VCAT from illumina or any vcf comparisioin tool.
# Requirements
* GATK-4.2.2.0
* Trimmomatic-0.36
* fastqc
* bwa-mem-2.2.1
* samtools
* bcftools
* picard-tools-2.2.1
* Annovar
# How to use
* Create Sample List: 
a simple text file list all fastq you need to analyze.
* How to run: as any shell script use ./name of the .sh file
