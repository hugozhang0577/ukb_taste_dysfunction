#!/bin/bash
# Variant QC in PLINK 2. Imputed variants (field 22828): INFO >= 0.8, MAF > 0.001, missingness < 0.01,
# HWE P > 1e-6, duplicates removed; written as 8-bit BGEN and indexed with bgenix.
# Array variants (field 22418): MAF > 0.01, missingness < 0.01, HWE P > 1e-6, then MAF > 0.05 and
# LD pruning (200 50 0.2) for the sparse genetic relationship matrix.
KEEP=taste_2w_strict_white_keep_ids.txt
PHENO=taste_2w_strict_white_gwas_pheno.txt

for CHR in {1..22}; do
  awk '$8 >= 0.8 {print $2}' ukb22828_c${CHR}_b0_v3.mfi.txt > c${CHR}_info08.snplist
  plink2 --bgen ukb22828_c${CHR}_b0_v3.bgen ref-first --sample ukb22828_c${CHR}_b0_v3.sample \
    --keep $KEEP --pheno $PHENO --pheno-name pheno --extract c${CHR}_info08.snplist \
    --maf 0.001 --geno 0.01 --hwe 1e-6 --rm-dup exclude-all \
    --export bgen-1.2 bits=8 --out ukb22828_c${CHR}_b0_v3_qc_filtered
  bgenix -g ukb22828_c${CHR}_b0_v3_qc_filtered.bgen -index -clobber
done

for CHR in {1..22}; do
  plink2 --bfile ukb22418_c${CHR}_b0_v2 --keep $KEEP --pheno $PHENO --pheno-name pheno --no-psam-pheno \
    --maf 0.01 --geno 0.01 --hwe 1e-6 --rm-dup exclude-all --make-bed --out ukb22418_c${CHR}_b0_v2_filt
done
for c in $(seq 2 22); do echo ukb22418_c${c}_b0_v2_filt; done > merge_list.txt
plink2 --bfile ukb22418_c1_b0_v2_filt --pmerge-list merge_list.txt bfile --make-bed --out maf_flt_22chroms

plink2 --bfile maf_flt_22chroms --maf 0.05 --indep-pairwise 200 50 0.2 --out grm_prune
plink2 --bfile maf_flt_22chroms --extract grm_prune.prune.in --keep $KEEP --make-bed \
  --out grm_for_saige_taste_2w_strict_white
