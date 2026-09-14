#!/usr/bin/env bash
#> ship gwas_sumstats.tsv.gz -> input/assoc_results/gwas_primary_sumstats.tsv.gz
for CHR in $(seq 1 22); do test -f taste_change_chr${CHR}.txt || { echo "missing chromosome ${CHR}" >&2; exit 1; }; done && { head -n 1 taste_change_chr1.txt; for CHR in $(seq 1 22); do tail -n +2 taste_change_chr${CHR}.txt; done; } | gzip -c > gwas_sumstats.tsv.gz && echo "variants: $(( $(zcat gwas_sumstats.tsv.gz | wc -l) - 1 ))"
# Export and dx upload to RAP  (the merged genome-wide summary statistics)
