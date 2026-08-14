#!/bin/bash
# =============================================================================
# GWAS — genotype preparation and variant-level QC
# =============================================================================
#
# Produces the two genotype inputs the association engines consume:
#   (A) analysis-ready per-chromosome IMPUTED bgen (UK Biobank field 22828, v3)
#   (B) the LD-pruned HARD-CALL PLINK set (from field 22418 array calls) that
#       SAIGE step 1 uses for the GRM and the variance-ratio markers
#
# Each block below is a command line to run in the Swiss Army Knife terminal,
# with the inputs for that step selected as the job's inputs. The files arrive in
# the working directory, so the commands use bare filenames and no paths. This is
# the one place in the package where paths are not written relative to
# $PROJECT_DIR: there is no project mount to be relative to inside a SAK job.
#
# Variant-level QC only. Individual-level sample QC is 01_sample_qc.R, which
# works from the UK Biobank central QC fields rather than PLINK --mind.
#
# Thresholds:
#   imputed  INFO >= 0.8 (.mfi column 8), MAF >= 0.001, call rate >= 0.99,
#            Hardy-Weinberg P > 1e-6 in controls, duplicates dropped, exported
#            as 8-bit bgen (what SAIGE reads) and bgenix-indexed
#   array    MAF >= 0.01, call rate >= 0.99, Hardy-Weinberg P > 1e-6, duplicates
#            dropped
#   GRM      MAF >= 0.05, --indep-pairwise 200 50 0.2
#
# INFO is read from the .mfi file because PLINK cannot see it. The Hardy-Weinberg
# test is control-only, which is why the phenotype is supplied at this stage
# rather than left to the association step.
#
# The association-level minor-allele-count filter is applied later by the engine
# (SAIGE --minMAC 20); the PRS base uses its own stricter filter.
#
# Two different sample files are used and they are NOT interchangeable:
# cohort3_2w_strict_keep_ids.txt is headerless FID/IID and is what --keep reads,
# while cohort3_2w_strict_gwas_pheno.txt carries the header and the `pheno`
# column that --pheno-name needs. Both come from
# ../preprocessing/04_gwas_sample_qc.R.
#
# --threads and --memory are resolved from the instance the job is running on, so
# the same line works on any instance type without editing.
# =============================================================================


# -----------------------------------------------------------------------------
# (A) imputed bgen — one job per chromosome
#     inputs: ukb22828_c${CHR}_b0_v3.{bgen,sample,mfi.txt}
#             cohort3_2w_strict_{keep_ids,gwas_pheno}.txt
# -----------------------------------------------------------------------------
# Single chromosome (CHR=22 shown; change the number for each job):

awk '$8 >= 0.8 {print $2}' ukb22828_c22_b0_v3.mfi.txt > c22_info08.snplist && plink2 --bgen ukb22828_c22_b0_v3.bgen ref-first --sample ukb22828_c22_b0_v3.sample --keep cohort3_2w_strict_keep_ids.txt --pheno cohort3_2w_strict_gwas_pheno.txt --pheno-name pheno --extract c22_info08.snplist --maf 0.001 --geno 0.01 --hwe 1e-6 --rm-dup exclude-all --export bgen-1.2 bits=8 --threads "$(nproc)" --memory "$(( $(awk '/MemTotal/{print int($2*0.8/1024)}' /proc/meminfo) ))" --out ukb22828_c22_b0_v3_qc_filtered && bgenix -g ukb22828_c22_b0_v3_qc_filtered.bgen -index -clobber

# All 22 in one job (only if every chromosome is attached to that job):

for CHR in {1..22}; do awk '$8 >= 0.8 {print $2}' ukb22828_c${CHR}_b0_v3.mfi.txt > c${CHR}_info08.snplist && plink2 --bgen ukb22828_c${CHR}_b0_v3.bgen ref-first --sample ukb22828_c${CHR}_b0_v3.sample --keep cohort3_2w_strict_keep_ids.txt --pheno cohort3_2w_strict_gwas_pheno.txt --pheno-name pheno --extract c${CHR}_info08.snplist --maf 0.001 --geno 0.01 --hwe 1e-6 --rm-dup exclude-all --export bgen-1.2 bits=8 --threads "$(nproc)" --memory "$(( $(awk '/MemTotal/{print int($2*0.8/1024)}' /proc/meminfo) ))" --out ukb22828_c${CHR}_b0_v3_qc_filtered && bgenix -g ukb22828_c${CHR}_b0_v3_qc_filtered.bgen -index -clobber; done


# -----------------------------------------------------------------------------
# (B) hard-call array — one job per chromosome
#     inputs: ukb22418_c${CHR}_b0_v2.{bed,bim,fam}
#             cohort3_2w_strict_{keep_ids,gwas_pheno}.txt
# -----------------------------------------------------------------------------

plink2 --bfile ukb22418_c22_b0_v2 --keep cohort3_2w_strict_keep_ids.txt --pheno cohort3_2w_strict_gwas_pheno.txt --pheno-name pheno --no-psam-pheno --maf 0.01 --geno 0.01 --hwe 1e-6 --rm-dup exclude-all --make-bed --threads "$(nproc)" --memory "$(( $(awk '/MemTotal/{print int($2*0.8/1024)}' /proc/meminfo) ))" --out ukb22418_c22_b0_v2_filt

# All 22 in one job:

for CHR in {1..22}; do plink2 --bfile ukb22418_c${CHR}_b0_v2 --keep cohort3_2w_strict_keep_ids.txt --pheno cohort3_2w_strict_gwas_pheno.txt --pheno-name pheno --no-psam-pheno --maf 0.01 --geno 0.01 --hwe 1e-6 --rm-dup exclude-all --make-bed --threads "$(nproc)" --memory "$(( $(awk '/MemTotal/{print int($2*0.8/1024)}' /proc/meminfo) ))" --out ukb22418_c${CHR}_b0_v2_filt; done


# -----------------------------------------------------------------------------
# Merge the 22 filtered array sets
#     inputs: ukb22418_c{1..22}_b0_v2_filt.{bed,bim,fam}  (from the step above)
# The merge list is written here so it names the files as they arrive.
# -----------------------------------------------------------------------------

for c in $(seq 2 22); do echo ukb22418_c${c}_b0_v2_filt; done > merge_list.txt && plink2 --bfile ukb22418_c1_b0_v2_filt --merge-list merge_list.txt --make-bed --threads "$(nproc)" --memory "$(( $(awk '/MemTotal/{print int($2*0.8/1024)}' /proc/meminfo) ))" --out maf_flt_22chroms


# -----------------------------------------------------------------------------
# LD-prune for the GRM
#     inputs: maf_flt_22chroms.{bed,bim,fam}, cohort3_2w_strict_keep_ids.txt
# -----------------------------------------------------------------------------

plink2 --bfile maf_flt_22chroms --maf 0.05 --indep-pairwise 200 50 0.2 --threads "$(nproc)" --memory "$(( $(awk '/MemTotal/{print int($2*0.8/1024)}' /proc/meminfo) ))" --out grm_prune && plink2 --bfile maf_flt_22chroms --extract grm_prune.prune.in --keep cohort3_2w_strict_keep_ids.txt --make-bed --threads "$(nproc)" --memory "$(( $(awk '/MemTotal/{print int($2*0.8/1024)}' /proc/meminfo) ))" --out grm_for_saige_cohort3_2w_strict
