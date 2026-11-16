#!/bin/bash
# Additive dosages of the two APOE-defining variants (rs429358, rs7412) from the QC-filtered imputed data.
printf "19 45411941 45411941 rs429358\n19 45412079 45412079 rs7412\n" > apoe_snps.txt
plink2 --bgen ukb22828_c19_b0_v3_qc_filtered.bgen ref-first --sample ukb22828_c19_b0_v3_qc_filtered.sample \
  --extract range apoe_snps.txt --export A --out apoe_genotypes
