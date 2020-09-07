#!/bin/bash
module load java

subset="~/RNA-seq/run7_yak_indicus/Sample_id"
chrom="~/RNA-seq/SNP_calling/All_SNPs/chromsome.list.txt"
run7="~/RNA-seq/run7_yak_indicus/genotypes"
Beagle="~/RNA-seq/beagle.21Sep19.ec3.jar"

######################################################
cd ${subset}
sub_set=`ls *| head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
cd /home/shuli.liu/bull_scr/RNA-seq/run7_yak_indicus/imputation/before_imputation
for region in `cat ${chrom}`
do
chr=`echo ${region} | cut -d ":" -f 1`
java -jar -Xmx120g -Djava.io.tmpdir=$TMPDIR ${Beagle} gt=Chr${chr}.${sub_set}.genotyped.vcf.gz ref=${run7}/Chr${chr}-Run7-TAU-Beagle-toDistribute.vcf.gz chrom=${chr} nthreads=18 out=../after_imputation/Chr${chr}.${sub_set}-beagle.vcf.gz
done
