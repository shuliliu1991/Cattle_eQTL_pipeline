#!/bin/bash
tissue=`ls *.vcf | head -n $SLURM_ARRAY_TASK_ID|tail -n 1 | cut -d "." -f 1`
bcftools view -q 0.01:minor -v snps -c 4:minor ${tissue}.vcf > ./${tissue}/${tissue}.filtered.vcf
cd ${tissue}
bgzip -cf ${tissue}.filtered.vcf > ${tissue}.filtered.vcf.gz
tabix -p vcf ${tissue}.filtered.vcf.gz
rm ccm.gds

Rscript /home/shuli.liu/commands/commands_for_RNA_seq_2/eQTL/prepare_afc_phenotype.r ${tissue}

sort -k1,1 -k2,2n ${tissue}.phenotype.afc.bed  > ${tissue}.phenotype.afc.s.bed 
bgzip -cf ${tissue}.phenotype.afc.s.bed > ${tissue}.phenotype.afc.bed.gz
tabix -p bed ${tissue}.phenotype.afc.bed.gz
rm ${tissue}.phenotype.afc.bed
rm ${tissue}.phenotype.afc.s.bed
sed '2d' ${tissue}.covariance.txt > ${tissue}.covariance2.txt
bgzip -cf ${tissue}.covariance2.txt > ${tissue}.covariance2.txt.gz

rm ${tissue}.covariance2.txt
sed '1d' ${tissue}.permutations.2rd.storey.txt | awk -F "_|\t| " '{print $1"\t"$6"_"$7"_"$8"_"$9"\t"$6"\t"$7}' | sed '1 i pid\tsid\tsid_chr\tsid_pos' > ${tissue}.permutation.sig.txt

${aFC} --vcf ${tissue}.filtered.vcf.gz --pheno ${tissue}.phenotype.afc.bed.gz --cov ${tissue}.covariance2.txt.gz --qtl ${tissue}.permutation.sig.txt --log_xform 0 --output ${tissue}.permutation.sig.afc.txt

rm ${tissue}.covariance2.txt.gz
rm ${tissue}.filtered.vcf
rm ${tissue}.filtered.vcf.gz
rm ${tissue}.filtered.vcf.gz.tbi
rm ${tissue}.phenotype.afc.bed.gz
