#!/bin/bash
tissue=`ls Adipose.vcf Blood.vcf Embryo.vcf Hypothalamus.vcf Ileum.vcf Intramuscular_fat.vcf Jejunum.vcf Leukocyte.vcf Liver.vcf Lung.vcf Lymph_node.vcf \
Macrophage.vcf Mammary.vcf Milk_cell.vcf Monocytes.vcf Muscle.vcf Ovary.vcf \
Oviduct.vcf Pituitary.vcf Rumen.vcf Salivary_gland.vcf Skin_fibroblast.vcf \
Testis.vcf Uterus.vcf | head -n $SLURM_ARRAY_TASK_ID|tail -n 1 | cut -d "." -f 1`
#mkdir ${tissue}
###minor allele frequencies >=0.01; with the minor allele observed in at least 4 samples.
bcftools view -q 0.01:minor -v snps -c 4:minor ${tissue}.vcf > ./${tissue}/${tissue}.filtered.vcf
cd ${tissue}
rm ccm.gds
bgzip -cf ${tissue}.filtered.vcf > ${tissue}.filtered.vcf.gz
tabix -p vcf ${tissue}.filtered.vcf.gz
rm ${tissue}.covariance.txt
Rscript /home/shuli.liu/commands/commands_for_RNA_seq_2/eQTL/prepare_for_eQTL.detection.r ${tissue} 10

sort -k1,1 -k2,2n ${tissue}.phenotype.bed  > ${tissue}.phenotype.s.bed 
bgzip -cf ${tissue}.phenotype.s.bed > ${tissue}.phenotype.bed.gz
tabix -p bed ${tissue}.phenotype.bed.gz
bgzip -cf ${tissue}.covariance.txt > ${tissue}.covariance.txt.gz
rm ${tissue}.phenotype.bed
rm ${tissue}.phenotype.s.bed

for j in $(seq 1 60); 
do
fastQTL.static --vcf ${tissue}.filtered.vcf.gz --bed ${tissue}.phenotype.bed.gz --cov ${tissue}.covariance.txt.gz  --out ${tissue}.nominals.chunk${j}.txt.gz --chunk $j 60& 
done
wait

zcat ${tissue}.nominals.chunk*.txt.gz | gzip -c > ${tissue}.nominals.2rd.txt.gz
rm ${tissue}.nominals.chunk*.txt.gz
#Rscript eQTL.p-value.nominal.BH.correction.r ${tissue} --no-save

for j in $(seq 1 60); 
do
fastQTL.static --vcf ${tissue}.filtered.vcf.gz --bed ${tissue}.phenotype.bed.gz --cov ${tissue}.covariance.txt.gz  --permute 1000 10000 --out ${tissue}.permutations.chunk${j}.txt.gz --chunk $j 60& 
done
wait

zcat ${tissue}.permutations.chunk*.txt.gz | gzip -c > ${tissue}.permutations.2rd.txt.gz
rm ${tissue}.permutations.chunk*.txt.gz
Rscript eQTL.p-value.permutation.correction.r ${tissue} --no-save

rm ${tissue}.filtered.vcf ${tissue}.filtered.vcf.gz ${tissue}.filtered.vcf.gz.tbi
