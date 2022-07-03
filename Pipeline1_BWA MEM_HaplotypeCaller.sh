#!/bin/bash

#==================================================================================
 # ======= Loading Modules======================================================
module use /app/Genome/modules
module load samtools-1.2

module use /app/utils/modules
module load jdk-1.8




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


#Depth of Coverage

/app/Genome/gatk-4.2.2.0/gatk DepthOfCoverage \
   -R ../../genome/genome.fa \
   -O Output_Germline_Aziz/${i}.coverage \
   -I Alignment_Germline_Aziz/$i.bam \
   -L ../../recalset/nextera_dna_exome_targeted_regions_manifest_v1_2.bed
   
   
   mv Output_Germline_Aziz/${i}.coverage.sample_interval_summary  Alignment_Germline_Aziz/${i}.coverage.csv

#==========================================================================
#==========================================================================
# Generating Ageregate Report
java -Xmx20g -jar /app/Genome/picard-tools-2.2.1/picard.jar CollectAlignmentSummaryMetrics \
          R=../../genome/genome.fa  \
          I=Alignment_Germline_Aziz/$i.bam \
          O=Output_Germline_Aziz/${i}_output1.txt
          
java -Xmx20g -jar /app/Genome/picard-tools-2.2.1/picard.jar CollectQualityYieldMetrics \
       I=Alignment_Germline_Aziz/$i.bam\
       O=Output_Germline_Aziz/${i}_output2.txt 
       
       java -Xmx20g -jar /app/Genome/picard-tools-2.2.1/picard.jar CollectVariantCallingMetrics  \
       DBSNP=../../recalset/dbsnp_132.hg19.vcf \
       I=Alignment_Germline_Aziz/${i}_HaplotypeCaller.vcf\
       O=Output_Germline_Aziz/${i}_output3.txt 
       
      java -Xmx20g -jar /app/Genome/picard-tools-2.2.1/picard.jar CollectHsMetrics   I=Alignment_Germline_Aziz/$i.bam  O=Output_Germline_Aziz/${i}_collect_hs_metrics.txt R=../../genome/genome.fa  BAIT_INTERVALS=../../recalset/Probe.interval_list  TARGET_INTERVALS=../../recalset/Traget.interval_list PER_BASE_COVERAGE=Output_Germline_Aziz/${i}_base-coverage PER_TARGET_COVERAGE=Output_Germline_Aziz/${i}_targetcoverage
     
     a="$(awk '{print $23}' Output_Germline_Aziz/${i}_collect_hs_metrics.txt)"
awk '$4 > $a * 2/10' Output_Germline_Aziz/${i}_base-coverage >Output_Germline_Aziz/${i}_Uniformity_RD.txt
     b=$(wc -l < Output_Germline_Aziz/${i}_Uniformity_RD.txt)
     c=$(wc -l < Output_Germline_Aziz/${i}_base-coverage)
     div=`echo $b / $c *100 |bc -l`

     
     
java -Xmx20g -jar /app/Genome/picard-tools-2.2.1/picard.jar CollectInsertSizeMetrics   I=Alignment_Germline_Aziz/$i.bam O=Output_Germline_Aziz/${i}_insert_size_metrics.txt   H=insert_size_histogram.pdf  M=0.5

/app/Genome/bedtools2/bin/bedtools intersect -abam Alignment_Germline_Aziz/$i.bam -b ../../recalset/nextera_dna_exome_targeted_regions_manifest_v1_2.bed -bed | wc > Output_Germline_Aziz/${i}_On_Target_Read.txt
      
  sed -sn 8p Output_Germline_Aziz/${i}_insert_size_metrics.txt>Output_Germline_Aziz/${i}_insert_size_metrics2.txt 

sed -sn 8,10p Output_Germline_Aziz/${i}_output1.txt> Output_Germline_Aziz/${i}_output11.txt
sed -sn 8p Output_Germline_Aziz/${i}_output2.txt> Output_Germline_Aziz/${i}_output22.txt
sed -sn 8p Output_Germline_Aziz/${i}_output3.txt.variant_calling_summary_metrics> Output_Germline_Aziz/${i}_output33.txt
sed -sn 8p Output_Germline_Aziz/${i}_collect_hs_metrics.txt> Output_Germline_Aziz/${i}_collect_hs_metrics4.txt 
sed -sn 9p Output_Germline_Aziz/${i}.rmdup.metrics.txt>Output_Germline_Aziz/${i}.rmdup.metrics2.txt


 
 
  echo '<!DOCTYPE html> <html> <head>  <style> tbody {display: flex;}table, th, td {display: flex;  border-collapse: collapse;} th, td {  border: 0.5px solid black;padding: 5px;  text-align: left;    }  img {display: block; margin-left: auto;margin-right: auto;}</style> </head>  <body> <img src="/fefs1/BDC/Analysis/NGS_Analysis/recalset/Logo.png" width="990"  height="130" > <br><br><br><br><br><br> <h1 align="center">Bioinformatics Unit</h1> <h2 align="center">GenaTi</h2> <h1 align="center" style="color:orange;"> Enrichment Sequencing Report</h1> <br> <h2 align="center" >Sample:'>Alignment_Germline_Aziz/${i}_Report.html
  

 echo $i >>Alignment_Germline_Aziz/${i}_Report.html 
 echo '</h2> <br> <p align="center">Workflow : WES-Germline Pipeline v2.0</p> <p align="center">Report Date:'>>Alignment_Germline_Aziz/${i}_Report.html  
 cat date1.txt >> Alignment_Germline_Aziz/${i}_Report.html  
 #<div  id="current_date" <p align="center"><script> date = new Date();year = date.getFullYear(); month = date.getMonth() + 1; day = date.getDate(); document.getElementById("current_date").innerHTML = month + "/" + day + "/" + year; </script> </div></p>'

#-Sample Information--------------------------------------------------------------- 
 echo ' <hr> <h2 align="left" style="color:blue;"> Sample Information </h2> <table > <tr > <th>Sample ID</th> <th>Run Folder </th> <th>Total PF Reads </th><th>Total PF UNIQUE READS* </th><th>Total PF Bases </th> <th>Percent Q30</th> <th>Mean Read Length</th> <th>Adapter trimmed</th></tr>'>> Alignment_Germline_Aziz/${i}_Report.html 
 echo ' <tr> <td>' >>Alignment_Germline_Aziz/${i}_Report.html 
 echo $i >>Alignment_Germline_Aziz/${i}_Report.html 
 echo ' </td> <td>'>> Alignment_Germline_Aziz/${i}_Report.html 
 echo "${PWD##*/}">> Alignment_Germline_Aziz/${i}_Report.html 
 echo ' </td> <td> '>> Alignment_Germline_Aziz/${i}_Report.html 
 PR="$(awk '{print $2}' Output_Germline_Aziz/${i}_output22.txt)"
 PB="$(awk '{print $4}' Output_Germline_Aziz/${i}_output22.txt)"
  RL="$(awk '{print $3}' Output_Germline_Aziz/${i}_output22.txt)"
   Q30="$(awk '{print $8}' Output_Germline_Aziz/${i}_output22.txt)"
   UR="$(awk '{print $8}'  Output_Germline_Aziz/${i}_collect_hs_metrics4.txt)"
  echo $PR>>Alignment_Germline_Aziz/${i}_Report.html 
   echo ' </td> <td>'>> Alignment_Germline_Aziz/${i}_Report.html 
     echo $UR>>Alignment_Germline_Aziz/${i}_Report.html 
   echo ' </td> <td>'>> Alignment_Germline_Aziz/${i}_Report.html 
   echo $PB>>Alignment_Germline_Aziz/${i}_Report.html 
    echo ' </td> <td>'>> Alignment_Germline_Aziz/${i}_Report.html 
 echo  " $Q30 / $PB *100" |bc -l  >> Alignment_Germline_Aziz/${i}_Report.html  
      echo ' </td> <td>'>> Alignment_Germline_Aziz/${i}_Report.html 
      echo $RL>>Alignment_Germline_Aziz/${i}_Report.html 
       echo ' </td><td> "Yes"</td> </table>'>> Alignment_Germline_Aziz/${i}_Report.html 
echo '* total PF_UNIQUE_READS: The number of PF reads that are not marked as duplicates.'>>    Alignment_Germline_Aziz/${i}_Report.html 

#------------------------------------------------------    
#---------Gender check-----------------------------
samtools idxstats Alignment_Germline_Aziz/$i.bam 2>/dev/null | awk '{ tmp=($3)/($2) ; printf"%s %0.5f\n", $1, tmp }' 2>/dev/null | head -23 > Output_Germline_Aziz/$i-autosomalcov.txt



min=$(echo $(sort -n -k 2 Output_Germline_Aziz/$i-autosomalcov.txt | head -n 1)| awk '{print $2}')

xcov=$(echo "scale=4; $(samtools idxstats Alignment_Germline_Aziz/$i.bam| grep "X" | cut -f 3)/$(samtools idxstats Alignment_Germline_Aziz/$i.bam | grep "X" | cut -f 2)" | bc)
ycov=$(echo "scale=4; $(samtools idxstats Alignment_Germline_Aziz/$i.bam | grep "Y" | cut -f 3)/$(samtools idxstats Alignment_Germline_Aziz/$i.bam | grep "Y" | cut -f 2)" | bc)


rat=$(echo "scale=4; ${xcov}/${ycov}" | bc)



test1=$(echo "$rat"'>'"3.0" '&&' "$rat"'<'"11.1" | bc -l)
if [ $test1 -eq 1 ]
then
   
    test1R=Female
else
test1=$(echo "$rat"'>'"0.3" '&&' "$rat"'<'"1.6" | bc -l)
if [ $test1 -eq 1 ]
then
    
     test1R=Male
fi
fi

test2=$(echo "$xcov"'>'"$min" | bc -l)
if [ $test2 -eq 1 ]
then
   test2R=Female
    
else   
   
    test2R=Male
fi

test3=$(echo "$ycov"'<('"$min"/10\)| bc -l)
if [ $test3 -eq 1 ]
then
   
    test3R=Female
else 
      echo Male
   
fi

echo "$min" > Output_Germline_Aziz/$i-autosomalcovResult.txt
echo "$xcov" >> Output_Germline_Aziz/$i-autosomalcovResult.txt
echo "$ycov" >> Output_Germline_Aziz/$i-autosomalcovResult.txt
echo "$rat" >> Output_Germline_Aziz/$i-autosomalcovResult.txt
echo "$test1R" >> Output_Germline_Aziz/$i-autosomalcovResult.txt
echo "$test2R" >> Output_Germline_Aziz/$i-autosomalcovResult.txt


 echo ' <hr> <h2 align="left" style="color:blue;"> Bam Estimated Gender </h2> <table > <tr > <th>Minimmum Autosomal Coverage</th> <th>Chr X Coverage </th> <th>Chr Y Coverage </th><th> X:Y-coverage ratio: </th> <th>test1-result(X/Y cov ratio): </th> <th>test2-result (compare X cov to min autosomal cov: </th> </tr>'>> Alignment_Germline_Aziz/${i}_Report.html 
echo ' <tr> ' >>Alignment_Germline_Aziz/${i}_Report.html 
awk 'BEGIN {}{for(i=1;i<=NF;i++)print "<td>" $i"</td>";print "</th>"} END{print "</tr></table>"}' Output_Germline_Aziz/$i-autosomalcovResult.txt>> Alignment_Germline_Aziz/${i}_Report.html

#------------------------------------------------------    
#---------Enrichment Summary-----------------------------
 echo ' <hr> <h2 align="left" style="color:blue;"> Enrichment Summary </h2> <table > <tr > <th>Target Manifest </th> <td>Custom (nextera-dna-exome-targeted-reg ions-manifest-v1-2)</td></tr> <tr> <th> Total Length of Targeted Reference </th><td>45,326,818</td></tr> </table>'>>Alignment_Germline_Aziz/${i}_Report.html 

#------------------------------------------------------    
#---------Read Level Enrichment -----------------------------
echo ' <hr> <h2 align="left" style="color:blue;"> Read Level Enrichment</h2> <table > <tr > <th>Total PF Aligned Reads </th> <td>'>> Alignment_Germline_Aziz/${i}_Report.html 
AR="$(awk 'FNR == 3 {print $6}' Output_Germline_Aziz/${i}_output11.txt)" 
echo  $AR>>  Alignment_Germline_Aziz/${i}_Report.html 
echo ' </td> </tr> <tr> <th> PCT PF Aligned Reads</th> <td> '>> Alignment_Germline_Aziz/${i}_Report.html 
PCTAR="$(awk 'FNR == 3 {print $7}' Output_Germline_Aziz/${i}_output11.txt)"
echo  "$PCTAR*100" |bc -l >>   Alignment_Germline_Aziz/${i}_Report.html 
echo ' </td> </tr> <tr> <th> PF Unique Aligned Reads</th> <td> '>> Alignment_Germline_Aziz/${i}_Report.html 
PFUR="$(awk '{print $11}'  Output_Germline_Aziz/${i}_collect_hs_metrics4.txt)"
echo  $PFUR>>  Alignment_Germline_Aziz/${i}_Report.html 
echo ' </td> </tr> <tr> <th> PCT PF Unique Aligned Reads</th> <td> '>> Alignment_Germline_Aziz/${i}_Report.html 
PCTPFUR="$(awk '{print $12}'  Output_Germline_Aziz/${i}_collect_hs_metrics4.txt)"
echo  "$PCTPFUR*100" |bc -l>>  Alignment_Germline_Aziz/${i}_Report.html 
echo ' </td> </tr> <tr> <th> Targeted Aligned Reads</th> <td> '>> Alignment_Germline_Aziz/${i}_Report.html 
TAR="$(awk '{print $1}'   Output_Germline_Aziz/${i}_On_Target_Read.txt)"
echo  $TAR>>  Alignment_Germline_Aziz/${i}_Report.html 
echo ' </td> </tr> <tr> <th> Read Enrichment</th> <td> '>> Alignment_Germline_Aziz/${i}_Report.html 
echo  " $TAR / $AR *100" |bc -l  >> Alignment_Germline_Aziz/${i}_Report.html  
echo ' </td> </tr> </table> '>> Alignment_Germline_Aziz/${i}_Report.html 




#------------------------------------------------------    
#---------Base Level Enrichment -----------------------------
echo ' <hr> <h2 align="left" style="color:blue;"> Base Level Enrichment</h2> <table > <tr > <th>Total PF Aligned Bases </th> <td>'>> Alignment_Germline_Aziz/${i}_Report.html 
AB="$(awk '{print $13}'  Output_Germline_Aziz/${i}_collect_hs_metrics4.txt)" 
echo  $AB>>  Alignment_Germline_Aziz/${i}_Report.html 
echo ' </td> </tr> <tr> <th>  PF Unique Aligned Reads</th> <td> '>> Alignment_Germline_Aziz/${i}_Report.html 
AUB="$(awk '{print $14}'  Output_Germline_Aziz/${i}_collect_hs_metrics4.txt)" 
echo  $AUB >>   Alignment_Germline_Aziz/${i}_Report.html 
echo ' </td> </tr> <tr> <th>  Targeted Aligned Bases </th> <td> '>> Alignment_Germline_Aziz/${i}_Report.html 
TAB="$(awk -F"," '{s+=$4}END{print s}' Output_Germline_Aziz/${i}.coverage)"
echo  $TAB >>   Alignment_Germline_Aziz/${i}_Report.html 
echo ' </td> </tr> <tr> <th> Base Enrichment</th> <td> '>> Alignment_Germline_Aziz/${i}_Report.html 
echo  " $TAB / $AB *100" |bc -l  >> Alignment_Germline_Aziz/${i}_Report.html  
echo ' </td> </tr> </table> '>> Alignment_Germline_Aziz/${i}_Report.html 

#------------------------------------------------------    
#--------- Coverage Summary -----------------------------
  meanTC=`echo  $TAB/45326818 |bc -l `
awk '$4 > $meanTC * 2/10' Output_Germline_Aziz/${i}_base-coverage >Output_Germline_Aziz/${i}_Uniformity_RD.txt
     b=$(wc -l < Output_Germline_Aziz/${i}_Uniformity_RD.txt)
     c=$(wc -l < Output_Germline_Aziz/${i}_base-coverage)
     div=`echo $b / $c \*100 | bc -l`




echo ' <hr> <h2 align="left" style="color:blue;"> Coverage Summary</h2> <table > <tr > <th> Mean target coverage  </th> <td>'>> Alignment_Germline_Aziz/${i}_Report.html 
echo   " $TAB /45326818" |bc -l >>  Alignment_Germline_Aziz/${i}_Report.html 
echo ' </td> </tr> <tr> <th>  Uniformirty </th> <td> ' >> Alignment_Germline_Aziz/${i}_Report.html 
echo  $div >>  Alignment_Germline_Aziz/${i}_Report.html 
echo ' </td> </tr> </table> <br>' >> Alignment_Germline_Aziz/${i}_Report.html 

#------------------------------------------------------    
#---------sequencing Depth Coverage -----------------------------
echo '  <table > <tr > <th> Depth of ccoverage </th> <th> 1x </th><th> 10x </th><th> 20x </th><th> 30x </th><th> 50x </th> <th> 100x </th><th> 200x </th></tr>'>> Alignment_Germline_Aziz/${i}_Report.html 

echo ' <tr> <th>Number of Targeted Bases Covered at or Above Indicated Depth of Coverage</th> <td>'>> Alignment_Germline_Aziz/${i}_Report.html 
awk -F"," '$2 >=1'  Output_Germline_Aziz/${i}.coverage > Output_Germline_Aziz/${i}_tt1.txt
echo "$(wc -l < Output_Germline_Aziz/${i}_tt1.txt)" >> Alignment_Germline_Aziz/${i}_Report.html 
echo '</td><td>'>> Alignment_Germline_Aziz/${i}_Report.html

awk -F"," '$2 >=10'  Output_Germline_Aziz/${i}.coverage > Output_Germline_Aziz/${i}_tt10.txt
echo "$(wc -l < Output_Germline_Aziz/${i}_tt10.txt)">> Alignment_Germline_Aziz/${i}_Report.html 
echo '</td><td>'>> Alignment_Germline_Aziz/${i}_Report.html

awk -F"," '$2 >=20'  Output_Germline_Aziz/${i}.coverage > Output_Germline_Aziz/${i}_tt20.txt
echo "$(wc -l < Output_Germline_Aziz/${i}_tt20.txt)">> Alignment_Germline_Aziz/${i}_Report.html 
echo '</td><td>'>> Alignment_Germline_Aziz/${i}_Report.html

awk -F"," '$2 >=30'  Output_Germline_Aziz/${i}.coverage > Output_Germline_Aziz/${i}_tt30.txt
echo "$(wc -l < Output_Germline_Aziz/${i}_tt30.txt)">> Alignment_Germline_Aziz/${i}_Report.html 
echo '</td><td>'>> Alignment_Germline_Aziz/${i}_Report.html

awk -F"," '$2 >=50'  Output_Germline_Aziz/${i}.coverage > Output_Germline_Aziz/${i}_tt50.txt
echo "$(wc -l < Output_Germline_Aziz/${i}_tt50.txt)">> Alignment_Germline_Aziz/${i}_Report.html 
echo '</td><td>'>> Alignment_Germline_Aziz/${i}_Report.html

awk -F"," '$2 >=100'  Output_Germline_Aziz/${i}.coverage > Output_Germline_Aziz/${i}_tt100.txt
echo "$(wc -l < Output_Germline_Aziz/${i}_tt100.txt)">> Alignment_Germline_Aziz/${i}_Report.html 
echo '</td><td>'>> Alignment_Germline_Aziz/${i}_Report.html

awk -F"," '$2 >=200'  Output_Germline_Aziz/${i}.coverage > Output_Germline_Aziz/${i}_tt200.txt
echo "$(wc -l < Output_Germline_Aziz/${i}_tt200.txt)">> Alignment_Germline_Aziz/${i}_Report.html 
echo '</td></tr>'>> Alignment_Germline_Aziz/${i}_Report.html



echo ' <tr> <th>Target Coverage at or Above Indicated Depth of Coverage</th> <td>'>> Alignment_Germline_Aziz/${i}_Report.html 
x1= `echo  $TBx1/45326818 |bc -l`
echo "  $(wc -l < Output_Germline_Aziz/${i}_tt1.txt) / 45326818 *100" |bc -l >> Alignment_Germline_Aziz/${i}_Report.html 
echo '</td><td>'>> Alignment_Germline_Aziz/${i}_Report.html

x10= `echo  $TBx10/45326818 |bc -l`
echo "  $(wc -l < Output_Germline_Aziz/${i}_tt10.txt) / 45326818 *100" |bc -l >> Alignment_Germline_Aziz/${i}_Report.html 
echo '</td><td>'>> Alignment_Germline_Aziz/${i}_Report.html

x20= `echo  $TBx20/45326818 |bc -l`
echo "  $(wc -l < Output_Germline_Aziz/${i}_tt20.txt) / 45326818 *100" |bc -l >> Alignment_Germline_Aziz/${i}_Report.html 
echo '</td><td>'>> Alignment_Germline_Aziz/${i}_Report.html

x30= `echo  $TBx30/45326818 |bc -l `
echo "  $(wc -l < Output_Germline_Aziz/${i}_tt30.txt) / 45326818 *100" |bc -l >> Alignment_Germline_Aziz/${i}_Report.html 
echo '</td><td>'>> Alignment_Germline_Aziz/${i}_Report.html

x50= `echo  $TBx50/45326818 |bc -l`
echo "  $(wc -l < Output_Germline_Aziz/${i}_tt50.txt) / 45326818 *100" |bc -l >> Alignment_Germline_Aziz/${i}_Report.html 
echo '</td><td>'>> Alignment_Germline_Aziz/${i}_Report.html

x100= `echo  $TBx100/45326818 |bc -l`
echo "  $(wc -l < Output_Germline_Aziz/${i}_tt100.txt) / 45326818 *100" |bc -l >> Alignment_Germline_Aziz/${i}_Report.html 
echo '</td><td>'>> Alignment_Germline_Aziz/${i}_Report.html

x200= `echo  $TBx200/45326818 |bc -l`
echo "  $(wc -l < Output_Germline_Aziz/${i}_tt200.txt) / 45326818 *100" |bc -l >> Alignment_Germline_Aziz/${i}_Report.html 
echo '</td></tr></table>'>> Alignment_Germline_Aziz/${i}_Report.html

#====================================================
#------------Insert size---------------------------------
 
 
echo '<h2 align="left" style="color:blue;"> Insert Size Metrics</h2> <table> <tr> <th>MEDIAN_INSERT_SIZE</th><th>MEDIAN_ABSOLUTE_DEVIATION</th> <th>MIN_INSERT_SIZE</th> <th>MAX_INSERT_SIZE</th> <th>MEAN_INSERT_SIZE</th> <th>STANDARD_DEVIATION</th> <th>READ_PAIRS</th><th>PAIR_ORIENTATION</th><th>WIDTH_OF_10_PERCENT</th><th>WIDTH_OF_20_PERCENT</th><th>WIDTH_OF_30_PERCENT</th> <th>WIDTH_OF_40_PERCENT</th> <th>WIDTH_OF_50_PERCENT</th> <th>WIDTH_OF_60_PERCENT</th> <th>WIDTH_OF_70_PERCENT</th> <th>WIDTH_OF_80_PERCENT</th> <th>WIDTH_OF_90_PERCENT</th><th>WIDTH_OF_99_PERCENT</th> </tr>' >> Alignment_Germline_Aziz/${i}_Report.html
  
  awk 'BEGIN {}{print "<tr>";for(i=1;i<=NF;i++)print "<td>" $i"</td>";print "</tr>"} END{print "</table>"}' Output_Germline_Aziz/${i}_insert_size_metrics2.txt>>Alignment_Germline_Aziz/${i}_Report.html
 
 #====================================================
#------------Duplication Metrics---------------------------------
  
  echo '<h2 align="left" style="color:blue;"> Duplication Metrics</h2> <table> <tr> <th>LIBRARY</th><th>UNPAIRED_READS_EXAMINED</th> <th>READ_PAIRS_EXAMINED</th> <th>UNMAPPED_READS</th> <th>UNPAIRED_READ_DUPLICATES</th> <th>READ_PAIR_DUPLICATES</th> <th>READ_PAIR_OPTICAL_DUPLICATES</th><th>PERCENT_DUPLICATION</th><th>ESTIMATED_LIBRARY_SIZE</th> </tr>' >> Alignment_Germline_Aziz/${i}_Report.html
  
awk 'BEGIN {}{print "<tr>";for(i=2;i<=NF;i++)print "<td>" $i"</td>";print "</tr>"} END{print "</table>"}' Output_Germline_Aziz/${i}.rmdup.metrics2.txt>>Alignment_Germline_Aziz/${i}_Report.html
 
 #====================================================
#------------ Calling Summary Metrics --------------------------------- 
  
 echo '<h2 align="left" style="color:blue;"> Variant Calling Summary Metrics</h2> <table> <tr> <th>TOTAL_SNPS</th><th>NUM_IN_DB_SNP</th> <th>NOVEL_SNPS</th> <th>FILTERED_SNPS</th> <th>PCT_DBSNP</th> <th>DBSNP_TITV</th><th>NOVEL_TITV</th><th>TOTAL_INDELS</th><th>NOVEL_INDELS</th><th>FILTERED_INDELS</th> <th>PCT_DBSNP_INDELS</th> <th>NUM_IN_DB_SNP_INDELS</th> <th>DBSNP_INS_DEL_RATIO</th> <th>NOVEL_INS_DEL_RATIO</th> <th>TOTAL_MULTIALLELIC_SNPS</th> <th>NUM_IN_DB_SNP_MULTIALLELIC</th><th>TOTAL_COMPLEX_INDELS</th> <th>NUM_IN_DB_SNP_COMPLEX_INDELS</th> <th>SNP_REFERENCE_BIAS</th><th>NUM_SINGLETONS</th> </tr>' >> Alignment_Germline_Aziz/${i}_Report.html

 awk 'BEGIN {}{print "<tr>";for(i=1;i<=NF;i++)print "<td>" $i"</td>";print "</tr>"} END{print "</table>"}' Output_Germline_Aziz/${i}_output33.txt>> Alignment_Germline_Aziz/${i}_Report.html

echo '<h2 align="left" style="color:blue;"> Analysis Details</h2> <h3 align="left">  Setting </h3> <table> <tr> <th>Setting Name</th><th>Reference Genome</th> <th>Targeted Regions</th> <th>Depth Threshold</th> <th>Flag PCR Duplicates </th> </tr> <tr> <td> Value </td> <td> Human (UCSC hg19) </td> <td> Custom (nextera-dna-exome-targeted-regions-manifest-v1-2)</td> <td> 2 </td> <td> True </td> </tr> </table>'>> Alignment_Germline_Aziz/${i}_Report.html

echo '<h3 align="left"> Software Versions </h3> <table> <tr> <th>Software</th> <th>bcl2fastq</th><th>BWA-mem (Aligner)</th> <th>Samtool</th> <th>Picard</th> <th>GATK-HaplotypeCaller </th> </tr> <tr> <td> Version </td> <td>v2-20 </td> <td> v-0.7.12</td> <td> v-1.2 </td>  <td> v-2.2.1 </td> <td>v-4.2.2</td></tr></table>'>> Alignment_Germline_Aziz/${i}_Report.html

 echo ' </body> </html>' >> Alignment_Germline_Aziz/${i}_Report.html
 #==========================================================================
 #==========================================================================
 sed -i '52 i ##reference=hg19' Alignment_Germline_Aziz/${i}_HaplotypeCaller.vcf
sed -i '52 i ##reference=hg19' Alignment_Germline_Aziz/${i}_HaplotypeCaller_Filtered.vcf
rm Alignment_Germline_Aziz/FastQC/*.zip
 
 #==========================================================================

# Validation 
#java -jar /app/Genome/GATK-3.8/GenomeAnalysisTK.jar -T VariantEval -R ../../genome/genome.fa -eval Alignment_Germline_Aziz/${i}_HaplotypeCaller.vcf  -D ../../recalset/dbsnp_132.hg19.vcf -o Alignment_Germline_Aziz/${i}_HaplotypeCaller_Evaluation.eval.grp
#==========================================================================

# Generating Summary Report
#===============================================================================
###### QC tarck

echo ' <tr> <td>' >>Output_Germline_Aziz/${i}_trackQC.html

echo $i>> Output_Germline_Aziz/${i}_trackQC.html
 
echo '</td><td>'>> Output_Germline_Aziz/${i}_trackQC.html

 Q30="$(awk '{print $8}' Output_Germline_Aziz/${i}_output22.txt)"
 PB="$(awk '{print $4}' Output_Germline_Aziz/${i}_output22.txt)"
echo  " scale=6; $Q30 / $PB *100" |bc   >> Output_Germline_Aziz/${i}_trackQC.html 
echo ' </td> <td>'>> Output_Germline_Aziz/${i}_trackQC.html
DR="$(awk '{print $9}' Output_Germline_Aziz/${i}.rmdup.metrics2.txt)"

echo  " scale=6; $DR *100" |bc -l >> Output_Germline_Aziz/${i}_trackQC.html
echo '</td><td>'>> Output_Germline_Aziz/${i}_trackQC.html
AR="$(awk 'FNR == 3 {print $6}' Output_Germline_Aziz/${i}_output11.txt)" 
echo  $AR>>  Output_Germline_Aziz/${i}_trackQC.html 
echo ' </td> <td> '>> Output_Germline_Aziz/${i}_trackQC.html 
TAR="$(awk '{print $1}'   Output_Germline_Aziz/${i}_On_Target_Read.txt)"
echo  $TAR>>  Output_Germline_Aziz/${i}_trackQC.html
echo ' </td> <td> '>> Output_Germline_Aziz/${i}_trackQC.html  
echo  " scale=6;$TAR / $AR *100" |bc -l  >> Output_Germline_Aziz/${i}_trackQC.html  

echo ' </td> <td> '>> Output_Germline_Aziz/${i}_trackQC.html 

TAB="$(awk -F"," '{s+=$4}END{print s}' Output_Germline_Aziz/${i}.coverage)"
echo   " scale=6; $TAB /45326818" |bc -l >>  Output_Germline_Aziz/${i}_trackQC.html
  
echo ' </td> <td> '>> Output_Germline_Aziz/${i}_trackQC.html 

  b=$(wc -l < Output_Germline_Aziz/${i}_Uniformity_RD.txt)
     c=$(wc -l < Output_Germline_Aziz/${i}_base-coverage)
     div=`echo $b / $c \*100 | bc -l`

echo  " scale=6; $div" |bc -l >>  Output_Germline_Aziz/${i}_trackQC.html 

echo ' </td> <td> '>> Output_Germline_Aziz/${i}_trackQC.html 

echo " scale=6; $(wc -l < Output_Germline_Aziz/${i}_tt30.txt) / 45326818 *100" |bc -l >> Output_Germline_Aziz/${i}_trackQC.html 

echo '</td><td>'>> Output_Germline_Aziz/${i}_trackQC.html 

echo " scale=6; $(wc -l < Output_Germline_Aziz/${i}_tt100.txt) / 45326818 *100" |bc -l >> Output_Germline_Aziz/${i}_trackQC.html 
 
echo '</td><td>'>> Output_Germline_Aziz/${i}_trackQC.html 


PB="$(awk '{print $6}' Output_Germline_Aziz/${i}_output33.txt)"
echo " scale=6; $PB" |bc -l >> Output_Germline_Aziz/${i}_trackQC.html

echo ' </td> </tr>  '>> Output_Germline_Aziz/${i}_trackQC.html 
# echo '  </table> '>> Output_Germline_Aziz/${i}_trackQC.html 
 
 
cat Output_Germline_Aziz/${i}_trackQC.html  >> Output_Germline_Aziz/tests.html

 echo '<!DOCTYPE html> <html> <head>  <style>  { table#ReportTable {border-width: 1px 1px 1px 1px;border-collapse: collapse;}table#ReportTable th {   }  img {display: block; margin-left: auto;margin-right: auto;} img {display: block; margin-left: auto;margin-right: auto;}</style> </head>  <body> <img src="Logo.png" width="990"  height="130" > <br><br><br><br><br><br> <br><br><br><br><br><br> <h1 align="center">Bioinformatics Unit</h1> <h2 align="center">GenaTi</h2> <h2 align="center" style="color:orange;"> Summary Information For Experiment:</h2> <h2 align="center" > '> Alignment_Germline_Aziz/Summary.html
 
 echo "${PWD##*/}" >> Alignment_Germline_Aziz/Summary.html
 echo '</h2>  <h2 align="center" >Report Date:<script> document.write(new Date().toLocaleDateString()); </script></h2> <hr>'>>Alignment_Germline_Aziz/Summary.html

string="${PWD##*/}"
remainder="$string"
first="${remainder%%_*}"; remainder="${remainder#*_}"
second="${remainder%%_*}"; remainder="${remainder#*_}"
third="${remainder%%_*}"; remainder="${remainder#*_}"
fourth="${remainder%%_*}";
echo '<h2 align="left" style="color:black;"> Chip Summary</h2> <table border="1" ID="ReportTable"> <tr> <th>Machine</th> <td>'>>Alignment_Germline_Aziz/Summary.html
echo $second >> Alignment_Germline_Aziz/Summary.html
echo ' </td> </tr> <tr> <th>Run Folder</th> <td>'>> Alignment_Germline_Aziz/Summary.html
echo "${PWD##*/}">> Alignment_Germline_Aziz/Summary.html
echo ' </td> </tr> <tr> <th>Chip ID</th> <td>'>> Alignment_Germline_Aziz/Summary.html
echo $fourth >> Alignment_Germline_Aziz/Summary.html
echo ' </td> </tr> </table>'>> Alignment_Germline_Aziz/Summary.html


echo '<h2 align="left" style="color:black;"> Sample Information</h2> <table border="1" ID="ReportTable"> <tr> <th>Sample Name</th></tr> ' >> Alignment_Germline_Aziz/Summary.html
awk ' BEGIN {}{print "  <tr>";for(i=1;i<=NF;i++)print "<td>" $i"</td>";} END{print "</tr></table>"}' SampleList.txt>> Alignment_Germline_Aziz/Summary.html
echo '<br>' >> Alignment_Germline_Aziz/Summary.html

sed '1,11d' Data/Intensities/BaseCalls/Reports/html/*/all/all/all/lane.html >> Alignment_Germline_Aziz/Summary.html
echo ' <h2 align="left" style="color:black;"> QC Mertrics </h2> <table table border="1" ID="ReportTable"> <tr>  <th>Sample ID</th><th>Percent Q30</th><th>the duplication rate</th> <th>total aligned reads </th> <th>total aligned targeted reads </th><th>% read enrichment </th><th>Mean region coverage depth </th> <th>Uniformity of coverage (Pct &gt; 0.2*mean)</th> <th>Depth of Sequencing Coverage at 30x</th> <th>Depth of Sequencing Coverage at 100x</th> <th>SNV transition/transversion ratio (Ti/Tv)</th></tr>'>> Alignment_Germline_Aziz/Summary.html

 
 sed -n 1,437p Output_Germline_Aziz/tests.html>> Alignment_Germline_Aziz/Summary.html
echo '</table></body> </html>' >>Alignment_Germline_Aziz/Summary.html

#===============================================================================
timer2=$(date +'%s')
diff=$(($timer2-$timer1))
echo "Total excution time: $(($diff / 3600 )) hours $((($diff % 3600) / 60)) minutes $(($diff % 60)) seconds" 