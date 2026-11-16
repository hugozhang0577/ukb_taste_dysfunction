#!/bin/bash
# Merge the 22 per-chromosome SAIGE outputs into one genome-wide file
# (used downstream as input/assoc_results/gwas_primary_sumstats.tsv.gz).
{ head -n 1 taste_change_chr1.txt; for CHR in $(seq 1 22); do tail -n +2 taste_change_chr${CHR}.txt; done; } |
  gzip -c > gwas_sumstats.tsv.gz
