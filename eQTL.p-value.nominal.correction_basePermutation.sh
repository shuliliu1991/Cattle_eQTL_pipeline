#!/bin/bash
tissue=`ls *.vcf | head -n $SLURM_ARRAY_TASK_ID|tail -n 1 | cut -d "." -f 1`
cd ${tissue}

Rscript eQTL.p-value.nominal.correction_basePermutation.r ${tissue}.permutations.2rd.txt.gz 0.05 ${tissue}.nominals.2rd.txt.gz ${tissue}.nominals.2rd.sig.txt
