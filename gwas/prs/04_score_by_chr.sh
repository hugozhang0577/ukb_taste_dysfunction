#!/bin/bash
# Clumping and thresholding with PRSice-2, per chromosome, in the 30% held-out subsample only: LD for clumping
# (250 kb, r2 0.1) is estimated from the held-out genotypes; scores at 14 P-value thresholds, no regression.
THREADS=16
for CHR in {1..22}; do
  PRSice_linux --base base_chr${CHR}.txt --snp SNP --chr CHR --bp BP --a1 A1 --a2 A2 --stat BETA --beta --pvalue P \
    --target ukb22828_c${CHR}_b0_v3_qc_filtered --type bgen --keep heldout_keep.txt \
    --pheno heldout_prsice_pheno.txt --pheno-col pheno --binary-target T \
    --clump-kb 250kb --clump-r2 0.1 --clump-p 1 \
    --bar-levels 5e-08,1e-07,1e-06,1e-05,0.0001,0.001,0.01,0.05,0.1,0.2,0.3,0.4,0.5,1 \
    --fastscore --all-score --print-snp --keep-ambig --no-regress --num-auto 22 --seed 42 \
    --thread $THREADS --out prs_chr${CHR}
done
