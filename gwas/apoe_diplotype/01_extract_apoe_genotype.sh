#!/bin/bash
# =============================================================================
# APOE diplotype — extract rs429358 and rs7412 dosages
# =============================================================================
#
# The two variants that define the APOE alleles, exported as an additive dosage
# .raw file for the diplotype call in 02_apoe_association.R.
#
#   rs429358  19:45411941  (GRCh37/hg19)
#   rs7412    19:45412079
#
# Extraction is by position, not by rsID: the QC re-export in ../02_variant_qc.sh
# does not guarantee that variant identifiers survive as rsIDs, whereas the
# coordinates are stable.
#
# Run in the Swiss Army Knife terminal with the chromosome 19 QC-filtered bgen
# and its .sample attached as job inputs. They arrive in the working directory,
# so the command uses bare filenames. The range file is written by the command
# itself rather than shipped as a separate input.
#
# --export A writes one dosage column per variant, named <ID>_<counted allele>.
# The command prints that header so the counted allele is visible in the job log;
# 02_apoe_association.R reads it out of the header rather than assuming it,
# because a silently flipped allele would invert the diplotype call.
#
# Inputs:  ukb22828_c19_b0_v3_qc_filtered.{bgen,sample}
# Output:  apoe_genotypes.raw
# =============================================================================

printf "19 45411941 45411941 rs429358\n19 45412079 45412079 rs7412\n" > apoe_snps.txt && plink2 --bgen ukb22828_c19_b0_v3_qc_filtered.bgen ref-first --sample ukb22828_c19_b0_v3_qc_filtered.sample --extract range apoe_snps.txt --export A --threads "$(nproc)" --memory "$(( $(awk '/MemTotal/{print int($2*0.8/1024)}' /proc/meminfo) ))" --out apoe_genotypes && head -1 apoe_genotypes.raw | tr "\t" "\n" | tail -n +7
