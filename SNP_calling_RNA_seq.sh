#!/bin/bash
folder_bam=""
folder_SNP=""
cd ${folder_bam}
bam=`ls *-STARAligned.sortedByCoord.out.bam | head -n $SLURM_ARRAY_TASK_ID|tail -n 1`
sample=`echo ${bam} | cut -d "-" -f 1 `
reference="/project/uvm_mckay/bismark_genome/Bovine_genome/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa"
dbSNP="/project/uvm_mckay/WGBS_others/tissue-specific-MHL/SNPS/ARS1.2_BQSR.vcf" #vcf format run7 #wget https://sites.ualberta.ca/~stothard/1000_bull_genomes/ARS1.2PlusY_BQSR.vcf.gz
interval="/project/uvm_mckay/WGBS_others/tissue-specific-MHL/SNPS/ARS1.2_BQSR.bed" #generated from vcf? awk '{ print $1,$2-1,$2,$3}' OFS="\t" ARS1.2_BQSR.vcf | sed '/#/d; /-1/d' > ARS1.2_BQSR.bed
picard="/software/7/apps/picard/64/2.9.2" 

###########################################################################
##1. add read groups, sort, mark duplicates, and create index

cd ${folder_snp}

mkdir ${sample}
java -jar -Xmx15g -Djava.io.tmpdir=$TMPDIR ${picard}/picard.jar AddOrReplaceReadGroups I=${folder_bam}/${bam} O=./${sample}/rg_added_${bam} \
RGID=4 RGLB=lib1 RGPL=illumina RGPU=run RGSM=20 CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate

cd ${sample}
java -jar -Xmx15g -Djava.io.tmpdir=$TMPDIR ${picard}/picard.jar MarkDuplicates I=rg_added_${bam} O=dedupped_${bam}  \
CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=marked_dup_metrics.txt

##2. Split N trim and reassign mapping qualities
#reassignOneMapping Quality: reassign all good alignments (MAPQ=225) to the default value of 60. 
gatk SplitNCigarReads \
   --spark-runner LOCAL \
   --TMP_DIR ${TMPDIR} \
      -R ${reference} \
      -I dedupped_${bam} \
      -O SplitNCigarReads_${bam}

 rm dedupped_${bam} #remove

##3. Base recalibration (BQSR). ( used or not, determine based on the run time; the effect is marginal)
 gatk BaseRecalibrator \
   --spark-runner LOCAL \
   --TMP_DIR ${TMPDIR} \
   -I SplitNCigarReads_${bam} \
   -R ${reference} \
   --known-sites ${dbSNP} \
   -O ${sample}-recal.table

 gatk ApplyBQSR \
  --spark-runner LOCAL \
   --TMP_DIR ${TMPDIR} \
   -R ${reference} \
   -I SplitNCigarReads_${bam} \
   --bqsr-recal-file ${sample}-recal.table \
   -O BQSR_${bam}
  
 rm SplitNCigarReads_${bam} #remove
 
##4. run the haplotypecaller. 

 gatk HaplotypeCaller  \
   --spark-runner LOCAL \
   --TMP_DIR ${TMPDIR} \
   -R ${reference} \
   -I BQSR_${bam} \
   -O ${sample}.vcf.gz \
   -dbsnp ${dbSNP} \
    -L ${interval} \
   --dont-use-soft-clipped-bases \
   --output-mode EMIT_ALL_CONFIDENT_SITES \
   -stand-call-conf 0
   
rm BQSR_${bam} #remove

 ##5. filter the snps
 
 gatk VariantFiltration \
 --TMP_DIR ${TMPDIR} \
 -R ${reference} \
 -V ${sample}.vcf.gz \
 -O ${sample}.filtered.vcf.gz \
 -window 35 -cluster 3 \
 --filter-name one \
 --filter-expression "FS>30.0" \
 --filter-name two \
 --filter-expression "QD<2.0" 

##6. selected variants, only keep snps
gatk SelectVariants \
-R ${reference} \
-V  ${sample}.filtered.vcf.gz \
--select-type-to-include SNP \
--exclude-filtered \
-O selected.${sample}.filtered.vcf.gz 
