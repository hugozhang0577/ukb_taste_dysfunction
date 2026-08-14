#!/bin/bash
# =============================================================================
# PRS — turn the discovery summary statistics into a PRSice base file
# =============================================================================
#
# The base file is the discovery GWAS: ../saige/ run with
# --phenoFile=discovery_70_gwas_pheno.txt (01_split_cohort.R writes it), then
# ../saige/03_merge_sumstats.sh. It must NOT be the full-cohort summary
# statistics — those were estimated in the same people the score is evaluated in.
#
# Column mapping. SAIGE reports the effect on Allele2, so A1 (the effect allele
# PRSice weights by) is Allele2 and A2 is Allele1. Getting this the wrong way
# round flips every weight and produces a score that is exactly as strong and
# points the wrong way, with nothing in the output to show it.
#
# The P pre-filter is a speed measure, not a threshold: variants above P = 0.5
# cannot enter any of the reported thresholds, so dropping them changes nothing
# and roughly halves the clumping work.
#
# Run these lines in the JupyterLab terminal.
#
# Input:  gwas_sumstats_discovery70.tsv.gz  (the 70% discovery scan)
# Output: gwas_base.txt, then gwas_base_chr{1..22}.txt
# =============================================================================

# Map SAIGE columns to what PRSice expects, dropping incomplete rows

zcat gwas_sumstats_discovery70.tsv.gz | awk 'BEGIN{OFS="\t"} NR==1{for(i=1;i<=NF;i++)h[$i]=i; print "SNP","CHR","BP","A1","A2","BETA","SE","P"; next} $h["p.value"]!="NA" && $h["BETA"]!="NA" && $h["SE"]!="NA" && $h["p.value"]+0<=0.5 {print $h["MarkerID"],$h["CHR"],$h["POS"],$h["Allele2"],$h["Allele1"],$h["BETA"],$h["SE"],$h["p.value"]}' > gwas_base.txt && wc -l gwas_base.txt


# Split by chromosome, so that scoring can go one chromosome at a time

for CHR in {1..22}; do awk -v c=${CHR} 'NR==1{print; next} $2==c' gwas_base.txt > gwas_base_chr${CHR}.txt; echo "chr${CHR}: $(( $(wc -l < gwas_base_chr${CHR}.txt) - 1 )) variants"; done
