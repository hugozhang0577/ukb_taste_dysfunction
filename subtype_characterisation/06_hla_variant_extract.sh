#!/bin/bash
# =============================================================================
# Subtype genetics — extract the chromosome 6 lead variant
# =============================================================================
#
# Additive dosages at the HLA class III lead variant, for the per-subtype
# association in 07_hla_persubtype.R (Figure 5D).
#
# The variant is taken as GIVEN. It was identified by the subtype-stratified
# genome-wide scan and the ASSET subset-based meta-analysis, which run separately
# from this package and are not reproduced by it: rs2071293, chr6:32,062,687
# (GRCh37), effect allele A, joint P = 3.0e-7. This script and 07 reproduce the
# panel's four per-subtype estimates at that locus, not the search that found it.
#
# Extraction is by position, not by rsID: the QC re-export in
# ../gwas/02_variant_qc.sh does not guarantee that variant identifiers survive as
# rsIDs, whereas the coordinates are stable.
#
# Run in the Swiss Army Knife terminal with the chromosome 6 QC-filtered bgen and
# its .sample attached as job inputs. They arrive in the working directory, so
# the command uses bare filenames and the whole-chromosome bgen is never read
# across the project mount. The range file is written by the command itself.
#
# --export A writes one dosage column per variant, named <ID>_<counted allele>.
# The command prints that header so the counted allele is visible in the job log,
# and 07_hla_persubtype.R stops if it is not the effect allele, rather than
# silently reporting a flipped odds ratio.
#
# Inputs:  ukb22828_c6_b0_v3_qc_filtered.{bgen,sample}
# Output:  hla_lead_genotypes.raw
# =============================================================================

printf "6 32062687 32062687 rs2071293\n" > hla_lead_snps.txt && plink2 --bgen ukb22828_c6_b0_v3_qc_filtered.bgen ref-first --sample ukb22828_c6_b0_v3_qc_filtered.sample --extract range hla_lead_snps.txt --export A --threads "$(nproc)" --memory "$(( $(awk '/MemTotal/{print int($2*0.8/1024)}' /proc/meminfo) ))" --out hla_lead_genotypes && head -1 hla_lead_genotypes.raw | tr "\t" "\n" | tail -n +7
