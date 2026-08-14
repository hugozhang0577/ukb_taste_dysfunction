#!/usr/bin/env bash
# =============================================================================
# SAIGE — concatenate the per-chromosome step-2 output into one summary file
# =============================================================================
#
# Step 2 writes one file per chromosome, each with its own header. Everything
# downstream — the diagnostics in 04, the FUMA input, the figure code — works on
# one genome-wide table, so the concatenation happens here rather than being
# repeated in each consumer.
#
# The header of chromosome 1 is kept and the rest are dropped. No filtering, no
# reformatting, no column changes.
#
# A missing chromosome stops the merge rather than producing a partial file that
# looks complete. A summary table silently short one chromosome is worse than no
# table at all, because nothing downstream can tell the difference.
#
# Inputs:  taste_change_chr{1..22}.txt
# Output:  gwas_sumstats.tsv.gz
# =============================================================================

for CHR in $(seq 1 22); do test -f taste_change_chr${CHR}.txt || { echo "missing chromosome ${CHR}" >&2; exit 1; }; done && { head -n 1 taste_change_chr1.txt; for CHR in $(seq 1 22); do tail -n +2 taste_change_chr${CHR}.txt; done; } | gzip -c > gwas_sumstats.tsv.gz && echo "variants: $(( $(zcat gwas_sumstats.tsv.gz | wc -l) - 1 ))"
# Export and dx upload to RAP  (the merged genome-wide summary statistics)
