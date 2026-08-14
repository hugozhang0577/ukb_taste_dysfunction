#!/bin/bash
# =============================================================================
# SAIGE — the primary genome-wide scan
# =============================================================================
#
# Two command lines to run in the Swiss Army Knife terminal, with SAIGE selected
# as the job's application. SAIGE is installed there, so these call
# step1_fitNULLGLMM.R and step2_SPAtests.R directly; the job's inputs arrive in
# the working directory, so the commands use bare filenames and no paths.
#
# Model: the primary specification reported in the Methods. Sex, the first ten
# genetic principal components, age at the taste assessment, genotyping batch,
# non-current smoking, alcohol risk and taste-affecting surgery, with relatedness
# carried by the sparse GRM.
#
# Which covariates are declared categorical matters and is easy to get wrong.
# batch_number is a batch label, not a quantity: entering it as a number would
# impose an ordering the assay does not have, so it goes in --qCovarColList along
# with the other non-numeric terms. Note also that sex_plink is genetic sex
# (field 22001, PLINK coding 1 male / 2 female), which is a different field and a
# different coding from the self-reported `sex` the panel scans use.
#
# The Supplementary Methods report three covariate variants of this same fit.
# Each is this command with one term changed in --covarColList (and, where the
# term is binary, in --qCovarColList), and a different --outputPrefix:
#   - drop  surg_taste_affecting_full      behaviour-free baseline
#   - add   BMI            (continuous)    metabolic sensitivity
#   - add   smell_any      (binary)        smell-adjusted sensitivity
# The smell-adjusted fit is a lower bound on the taste-specific signal rather
# than an alternative primary estimate: adjusting for concurrent smell change
# removes the chemosensory component the two share.
#
# Step 2 is per chromosome. --minMAC 20 sets the association-level variant
# filter; Firth refinement is applied to variants reaching P < 0.05.
#
# Inputs (step 1): sparseGRM_relatednessCutoff_0.125_5000_randomMarkersUsed.sparseGRM.mtx
#                  ...sparseGRM.mtx.sampleIDs.txt
#                  random_50000_markers.{bed,bim,fam}
#                  cohort3_2w_strict_gwas_pheno.txt
# Inputs (step 2): step1_taste_change.rda, step1_taste_change.varianceRatio.txt
#                  ukb22828_c{1..22}_b0_v3_qc_filtered.{bgen,bgen.bgi,sample}
# Output:          taste_change_chr{1..22}.txt  -> 03_merge_sumstats.sh
# =============================================================================


# -----------------------------------------------------------------------------
# Step 1 — fit the null logistic mixed model
# -----------------------------------------------------------------------------

step1_fitNULLGLMM.R --sparseGRMFile=sparseGRM_relatednessCutoff_0.125_5000_randomMarkersUsed.sparseGRM.mtx --sparseGRMSampleIDFile=sparseGRM_relatednessCutoff_0.125_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt --useSparseGRMtoFitNULL=TRUE --plinkFile=random_50000_markers --phenoFile=cohort3_2w_strict_gwas_pheno.txt --phenoCol=pheno --sampleIDColinphenoFile=IID --traitType=binary --covarColList=sex_plink,batch_number,smoking,drink,surg_taste_affecting_full,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,age --qCovarColList=sex_plink,batch_number,smoking,drink,surg_taste_affecting_full --skipVarianceRatioEstimation=FALSE --IsOverwriteVarianceRatioFile=TRUE --nThreads="$(nproc)" --outputPrefix=step1_taste_change


# -----------------------------------------------------------------------------
# Step 2 — single-variant association, chromosomes 1 to 22
# -----------------------------------------------------------------------------

for CHR in {1..22}; do step2_SPAtests.R --bgenFile=ukb22828_c${CHR}_b0_v3_qc_filtered.bgen --bgenFileIndex=ukb22828_c${CHR}_b0_v3_qc_filtered.bgen.bgi --sampleFile=ukb22828_c${CHR}_b0_v3_qc_filtered.sample --GMMATmodelFile=step1_taste_change.rda --varianceRatioFile=step1_taste_change.varianceRatio.txt --SAIGEOutputFile=taste_change_chr${CHR}.txt --chrom=${CHR} --AlleleOrder=ref-first --minMAC=20 --LOCO=FALSE --is_Firth_beta=TRUE --pCutoffforFirth=0.05 --is_fastTest=FALSE --SPAcutoff=0.5; done


# One chromosome at a time, if step 2 is split across jobs (CHR=22 shown):

step2_SPAtests.R --bgenFile=ukb22828_c22_b0_v3_qc_filtered.bgen --bgenFileIndex=ukb22828_c22_b0_v3_qc_filtered.bgen.bgi --sampleFile=ukb22828_c22_b0_v3_qc_filtered.sample --GMMATmodelFile=step1_taste_change.rda --varianceRatioFile=step1_taste_change.varianceRatio.txt --SAIGEOutputFile=taste_change_chr22.txt --chrom=22 --AlleleOrder=ref-first --minMAC=20 --LOCO=FALSE --is_Firth_beta=TRUE --pCutoffforFirth=0.05 --is_fastTest=FALSE --SPAcutoff=0.5
