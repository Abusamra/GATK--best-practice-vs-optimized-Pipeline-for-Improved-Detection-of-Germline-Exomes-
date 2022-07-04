#!/bin/bash

#==================================================================================
 # ======= Loading Modules======================================================
module use /app/Genome/modules
module load samtools-1.2

module use /app/utils/modules
module load jdk-1.8

module use /app/utils/modules/ 
module load gcc-5.3.0 boost/gcc/1_69_0



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
java -jar /app/Genome/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 20 ${i}_L00${j}_R1_001.fastq.gz ${i}_L00${j}_R2_001.fastq.gz Output_Germline_Aziz/${i}_L00${j}_R1_001.trim.fastq.gz Output_Germline_Aziz/${i}_L00${j}_R1_001_1un.trim.fastq.gz Output_Germline_Aziz/${i}_L00${j}_R2_001.trim.fastq.gz Output_Germline_Aziz/${i}_L00${j}_R2_001_2un.trim.fastq.gz SLIDINGWINDOW:4:20 MINLEN:20 
 
/app/Genome/FastQC/fastqc Output_Germline_Aziz/${i}_L00${j}_R1_001.trim.fastq.gz Output_Germline_Aziz/${i}_L00${j}_R2_001.trim.fastq.gz  -o Alignment_Germline_Aziz/FastQC 
 
#/app/Genome/dragmap-1.2.1/dragen-os -r /fefs1/BDC/Analysis/NGS_Analysis/recalset/ -1 Output_Germline_Aziz/${i}_L00${j}_R1_001.trim.fastq.gz -2 Output_Germline_Aziz/${i}_L00${j}_R2_001.trim.fastq.gz  --output-directory  Output_Germline_Aziz/  --output-file-prefix ${i}_L00${j}
/app/Genome/bwa-mem-2.2.1/bwa-mem2 mem -t 20 -R  "@RG\tID:$i\tSM:$i\tPL:ILLUMINA\tPI:330" ../../recalset/genome/genome.fa Output_Germline_Aziz/${i}_L00${j}_R1_001.trim.fastq.gz Output_Germline_Aziz/${i}_L00${j}_R2_001.trim.fastq.gz > Output_Germline_Aziz/${i}_L00${j}.sam
      done 


      for f in 1 2 
      do

     samtools view -b -S -o  Output_Germline_Aziz/${i}_L00$f.bam Output_Germline_Aziz/${i}_L00$f.sam
		 samtools sort  Output_Germline_Aziz/${i}_L00$f.bam Output_Germline_Aziz/${i}_L00$f.sorted
      done
#========================================================================================
#========= merging bams
  
samtools merge  Output_Germline_Aziz/$i.bam Output_Germline_Aziz/${i}_L001.sorted.bam Output_Germline_Aziz/${i}_L002.sorted.bam 
 

#============================================================================================
#======== Mark duplication
samtools index  Output_Germline_Aziz/$i.bam

/app/Genome/gatk-4.2.2.0/gatk MarkDuplicatesSpark  -I Output_Germline_Aziz/$i.bam -O Output_Germline_Aziz/$i.rmdup.bam --remove-all-duplicates false -M Output_Germline_Aziz/$i.rmdup.metrics.txt -L ../../recalset/nextera_dna_exome_targeted_regions_manifest_v1_2.bed --spark-master local[10]


# indexing       
samtools index  Output_Germline_Aziz/$i.rmdup.bam

timer3=$(date +'%s')
diff=$(($timer3-$timer1))
echo "Aligner excution time: $(($diff / 3600 )) hours $((($diff % 3600) / 60)) minutes $(($diff % 60)) seconds"
#===============================================
### BQSR=======================================

/app/Genome/gatk-4.2.2.0/gatk BQSRPipelineSpark -R ../../genome/genome.fa -I Output_Germline_Aziz/$i.rmdup.bam  -L ../../recalset/nextera_dna_exome_targeted_regions_manifest_v1_2.bed --known-sites ../../recalset/dbsnp_132.hg19.vcf -known-sites ../../recalset/Mills_Devine_2hit.indels.hg19.vcf --known-sites ../../recalset/1000G_omni2.5.hg19.sites.vcf --known-sites ../../recalset/hapmap_3.3.hg19.sites.vcf -O Alignment_Germline_Aziz/$i.bam --spark-master local[10]
#==============================================

#===============================================================================
#===============================================================================
timer4=$(date +'%s')
#Variant Calling
#Germaline

/app/Genome/gatk-4.2.2.0/gatk CalibrateDragstrModel -I Alignment_Germline_Aziz/$i.bam  -O Output_Germline_Aziz/$i-DRgstr.txt -R ../../genome/genome.fa -str  ../../recalset/hg19.str.zip

/app/Genome/gatk-4.2.2.0/gatk  HaplotypeCallerSpark -R ../../genome/genome.fa -I Alignment_Germline_Aziz/$i.bam   -O Output_Germline_Aziz/${i}_HaplotypeCall.vcf -L ../../recalset/nextera_dna_exome_targeted_regions_manifest_v1_2.bed --dragstr-params-path Output_Germline_Aziz/$i-DRgstr.txt  --dragen-mode true  --spark-master local[10]


#===================================================================================================
# FILTERING 
# Filtering by variants features
/app/Genome/GATK-4.1.4/gatk  VariantFiltration -R ../../genome/genome.fa -V  Output_Germline_Aziz/${i}_HaplotypeCall.vcf -filter "QUAL < 10.4139" --filter-name "lowQUAL" --filter-expression "DP < 2"   --filter-name "LowDP" --filter-expression "MQ < 40.0"  --filter-name "LowMQ" --filter-expression "QD < 2.0"  --filter-name "LowQD" -O Alignment_Germline_Aziz/${i}_HaplotypeCaller.vcf



/app/Genome/bcftools-1.3/bin/bcftools view -f PASS  Alignment_Germline_Aziz/${i}_HaplotypeCaller.vcf > Output_Germline_Aziz/${i}_HaplotypeCaller_Filtered2.vcf

 perl /fefs1/BDC/annovar/table_annovar.pl Output_Germline_Aziz/${i}_HaplotypeCaller_Filtered2.vcf  /fefs1/BDC/annovar/humandb/ -buildver hg19 -out Output_Germline_Aziz/${i}_HaplotypeCaller_Filtered3.vcf -protocol exac03,esp6500siv2_all,gnomad_exome -operation f,f,f -nastring . -vcfinput

/app/Genome/bcftools-1.3/bin/bcftools filter -e 'INFO/ExAC_ALL>0.02 | INFO/gnomAD_exome_ALL>0.02 | INFO/esp6500siv2_all>0.02' Output_Germline_Aziz/${i}_HaplotypeCaller_Filtered3.vcf.hg19_multianno.vcf > Alignment_Germline_Aziz/${i}_HaplotypeCaller_Filtered.vcf

timer5=$(date +'%s')
diff=$(($timer5-$timer4))
echo "Caller excution time: $(($diff / 3600 )) hours $((($diff % 3600) / 60)) minutes $(($diff % 60)) seconds"
#===============================================================================
#Depth of Coverage

/app/Genome/gatk-4.2.2.0/gatk DepthOfCoverage \
   -R ../../genome/genome.fa \
   -O Output_Germline_Aziz/${i}.coverage \
   -I Alignment_Germline_Aziz/$i.bam \
   -L ../../recalset/nextera_dna_exome_targeted_regions_manifest_v1_2.bed
   
   
   mv Output_Germline_Aziz/${i}.coverage.sample_interval_summary  Alignment_Germline_Aziz/${i}.coverage.csv

#===============================================================================


timer2=$(date +'%s')
diff=$(($timer2-$timer1))
echo "Total excution time: $(($diff / 3600 )) hours $((($diff % 3600) / 60)) minutes $(($diff % 60)) seconds"
#==========================================================================
#==========================================================================

 