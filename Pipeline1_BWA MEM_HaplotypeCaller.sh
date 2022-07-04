#!/bin/bash

#==================================================================================
 # ======= Loading Modules======================================================
module use /app/Genome/modules
module load samtools-1.2

module use /app/utils/modules
module load jdk-1.8

S


timer1=$(date +'%s')

#==================================================================================
# ======= Read Fastq from Sample List=======================================

parameters=`sed -n "${PBS_ARRAY_INDEX} p" SampleList.txt`
 i=$parameters


mkdir Alignment_Germline_Aziz
mkdir Output_Germline_Aziz
mkdir Alignment_Germline_Aziz/FastQC

#==================================================================================
# ======= SAM and BAM converting =======================================
  for j in 1 2 
      do
java -jar /app/Genome/Trimmomatic-0.36/trimmomatic-0.36.jar PE ${i}_L00${j}_R1_001.fastq.gz ${i}_L00${j}_R2_001.fastq.gz Output_Germline_Aziz/${i}_L00${j}_R1_001.trim.fastq.gz Output_Germline_Aziz/${i}_L00${j}_R1_001_1un.trim.fastq.gz Output_Germline_Aziz/${i}_L00${j}_R2_001.trim.fastq.gz Output_Germline_Aziz/${i}_L00${j}_R2_001_2un.trim.fastq.gz SLIDINGWINDOW:4:20 MINLEN:20 
 
/app/Genome/FastQC/fastqc Output_Germline_Aziz/${i}_L00${j}_R1_001.trim.fastq.gz Output_Germline_Aziz/${i}_L00${j}_R2_001.trim.fastq.gz  -o Alignment_Germline_Aziz/FastQC 
 
/app/Genome/bwa-0.7.12/bwa mem -M -t 16 -R  "@RG\tID:$i\tSM:$i\tPL:ILLUMINA\tPI:330" ../../genome/genome.fa Output_Germline_Aziz/${i}_L00${j}_R1_001.trim.fastq.gz Output_Germline_Aziz/${i}_L00${j}_R2_001.trim.fastq.gz > Output_Germline_Aziz/${i}_L00${j}.sam
      done 


      for f in 1 2 
      do

     samtools view -b -S -o  Output_Germline_Aziz/${i}_L00$f.bam Output_Germline_Aziz/${i}_L00$f.sam
		 samtools sort Output_Germline_Aziz/${i}_L00$f.bam Output_Germline_Aziz/${i}_L00$f.sorted
      done
#========================================================================================
#========= merging bams
  
samtools merge  Output_Germline_Aziz/$i.bam Output_Germline_Aziz/${i}_L001.sorted.bam Output_Germline_Aziz/${i}_L002.sorted.bam 
 

#============================================================================================
#======== Mark duplication
java -Xmx20g -jar /app/Genome/picard-tools-2.2.1/picard.jar MarkDuplicates INPUT=Output_Germline_Aziz/$i.bam OUTPUT=Output_Germline_Aziz/$i.rmdup.bam REMOVE_DUPLICATES=FALSE METRICS_FILE=Output_Germline_Aziz/$i.rmdup.metrics.txt

# BAM QC
#/app/Genome/qualimap_v2.2.1/qualimap bamqc -bam Alignment_Germline_Aziz/$i.bam  -outfile ${i}_BAMQC  -outformat PDF --java-mem-size=4G 


# indexing       
samtools index Output_Germline_Aziz/$i.rmdup.bam

#===============================================

/app/Genome/gatk-4.1.8.1/gatk BaseRecalibrator -R ../../genome/genome.fa -I Output_Germline_Aziz/$i.rmdup.bam  -L ../../recalset/nextera_dna_exome_targeted_regions_manifest_v1_2.bed --known-sites ../../recalset/dbsnp_132.hg19.vcf -known-sites ../../recalset/Mills_Devine_2hit.indels.hg19.vcf --known-sites ../../recalset/1000G_omni2.5.hg19.sites.vcf --known-sites ../../recalset/hapmap_3.3.hg19.sites.vcf -O Output_Germline_Aziz/${i}.recal_data.table 

 /app/Genome/gatk-4.1.8.1/gatk  ApplyBQSR  -R ../../genome/genome.fa  -I Output_Germline_Aziz/$i.rmdup.bam  --bqsr-recal-file Output_Germline_Aziz/${i}.recal_data.table  -O Alignment_Germline_Aziz/$i.bam

#==============================================

#===============================================================================

#===============================================================================

#Variant Calling
#Germaline
/app/Genome/GATK-4.1.4/gatk  HaplotypeCaller -R ../../genome/genome.fa -I Alignment_Germline_Aziz/$i.bam  --dbsnp ../../recalset/dbsnp_132.hg19.vcf  -L ../../recalset/nextera_dna_exome_targeted_regions_manifest_v1_2.bed -O Output_Germline_Aziz/${i}_HaplotypeCall.vcf -A StrandBiasBySample



#===================================================================================================
# FILTERING 
# Filtering by variants features
/app/Genome/GATK-4.1.4/gatk  VariantFiltration -R ../../genome/genome.fa -V  Output_Germline_Aziz/${i}_HaplotypeCall.vcf -filter "QUAL < 10.4139" --filter-name "lowQUAL" --filter-expression "DP < 2"   --filter-name "LowDP" --filter-expression "MQ < 40.0"  --filter-name "LowMQ" --filter-expression "QD < 2.0"  --filter-name "LowQD" -O Alignment_Germline_Aziz/${i}_HaplotypeCaller.vcf



/app/Genome/bcftools-1.3/bin/bcftools view -f PASS  Alignment_Germline_Aziz/${i}_HaplotypeCaller.vcf > Output_Germline_Aziz/${i}_HaplotypeCaller_Filtered2.vcf

 perl /fefs1/BDC/annovar/table_annovar.pl Output_Germline_Aziz/${i}_HaplotypeCaller_Filtered2.vcf  /fefs1/BDC/annovar/humandb/ -buildver hg19 -out Output_Germline_Aziz/${i}_HaplotypeCaller_Filtered3.vcf -protocol exac03,esp6500siv2_all,gnomad_exome -operation f,f,f -nastring . -vcfinput

/app/Genome/bcftools-1.3/bin/bcftools filter -e 'INFO/ExAC_ALL>0.02 | INFO/gnomAD_exome_ALL>0.02 | INFO/esp6500siv2_all>0.02' Output_Germline_Aziz/${i}_HaplotypeCaller_Filtered3.vcf.hg19_multianno.vcf > Alignment_Germline_Aziz/${i}_HaplotypeCaller_Filtered.vcf



#===============================================================================
timer2=$(date +'%s')
diff=$(($timer2-$timer1))
echo "Total excution time: $(($diff / 3600 )) hours $((($diff % 3600) / 60)) minutes $(($diff % 60)) seconds" 