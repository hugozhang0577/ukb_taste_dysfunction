#!/bin/bash
# =============================================================================
# PRS — clumping-and-thresholding scores, one chromosome at a time
# =============================================================================
#
# PRSice-2 computes the score at every P-value threshold. Clumping keeps the most
# significant variant in each 250 kb window and drops anything correlated with it
# above r2 = 0.1, which is what stops a single strong locus contributing once per
# variant in linkage disequilibrium with it.
#
# --no-regress is deliberate: this step only writes scores. The association
# between score and outcome is fitted in 06_regression.R, in the held-out subset
# alone. Letting PRSice regress here would fit it in whoever happens to be in the
# target file, which is everyone, which is the leakage this design exists to
# avoid.
#
# Scoring the whole cohort is fine and is what the merge step expects: the
# weights come from the discovery 70% only, so the held-out rows are genuinely
# out of sample. It is the *evaluation* that has to be restricted, not the
# scoring.
#
# The phenotype file is passed because PRSice requires one; with --no-regress it
# is used only to decide which samples are scored.
#
# Reading a whole-chromosome bgen across the read-only project mount is slow, so
# download the chromosome to the workspace first. Sizes are large; delete each
# one after scoring if disk is tight.
#
# Run these lines in the JupyterLab terminal.
#
# Inputs:  gwas_base_chr{1..22}.txt, ukb22828_c{1..22}_b0_v3_qc_filtered.{bgen,sample}
#          cohort3_2w_strict_gwas_pheno.txt
# Output:  prs_chr{1..22}.all_score
# =============================================================================

# One chromosome (CHR=22 shown), including the download and cleanup

dx download "/gwas/cohort_primary/genotype/ukb22828_c22_b0_v3_qc_filtered.bgen" -o ukb22828_c22_b0_v3_qc_filtered.bgen && dx download "/gwas/cohort_primary/genotype/ukb22828_c22_b0_v3_qc_filtered.sample" -o ukb22828_c22_b0_v3_qc_filtered.sample && ./PRSice_linux --base gwas_base_chr22.txt --snp SNP --chr CHR --bp BP --a1 A1 --a2 A2 --stat BETA --pvalue P --target ukb22828_c22_b0_v3_qc_filtered --type bgen --pheno cohort3_2w_strict_gwas_pheno.txt --pheno-col pheno --binary-target T --clump-kb 250 --clump-r2 0.1 --clump-p 1 --bar-levels 1e-5,1e-4,1e-3,0.01,0.05,0.1,0.2,0.3,0.4,0.5,1 --fastscore --all-score --print-snp --keep-ambig --no-regress --seed 42 --thread "$(nproc)" --out prs_chr22 && rm -f ukb22828_c22_b0_v3_qc_filtered.bgen


# All 22 in sequence. This is long-running; start it with nohup so that it
# survives the notebook session being closed, and watch prs_all_chr.log.

nohup bash -c 'for CHR in {1..22}; do dx download "/gwas/cohort_primary/genotype/ukb22828_c${CHR}_b0_v3_qc_filtered.bgen" -o ukb22828_c${CHR}_b0_v3_qc_filtered.bgen && dx download "/gwas/cohort_primary/genotype/ukb22828_c${CHR}_b0_v3_qc_filtered.sample" -o ukb22828_c${CHR}_b0_v3_qc_filtered.sample && ./PRSice_linux --base gwas_base_chr${CHR}.txt --snp SNP --chr CHR --bp BP --a1 A1 --a2 A2 --stat BETA --pvalue P --target ukb22828_c${CHR}_b0_v3_qc_filtered --type bgen --pheno cohort3_2w_strict_gwas_pheno.txt --pheno-col pheno --binary-target T --clump-kb 250 --clump-r2 0.1 --clump-p 1 --bar-levels 1e-5,1e-4,1e-3,0.01,0.05,0.1,0.2,0.3,0.4,0.5,1 --fastscore --all-score --print-snp --keep-ambig --no-regress --seed 42 --thread "$(nproc)" --out prs_chr${CHR} && rm -f ukb22828_c${CHR}_b0_v3_qc_filtered.bgen; done' > prs_all_chr.log 2>&1 &


# Check what finished

ls prs_chr*.all_score | wc -l
