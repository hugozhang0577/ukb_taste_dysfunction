# Shared settings for the metabolome-wide scans: inputs, outcome, factor coding, the
# 19 base covariates, the 37 CE-IVD certified Nightingale biomarkers and the model runner.
library(data.table)
source("mwas/mwas_glm_analysis.R")

NMR_FILES   <- c(group1 = "input/analysis_ready/metabolomics_group1.csv",
                 group2 = "input/analysis_ready/metabolomics_group2.csv",
                 group3 = "input/analysis_ready/metabolomics_group3.csv")
PHENO_FILES <- c(group1 = "input/analysis_ready/phenotype_group1.csv",
                 group2 = "input/analysis_ready/phenotype_group2.csv",
                 group3 = "input/analysis_ready/phenotype_group3.csv")
PHENO_SUPP  <- "input/analysis_ready/nmr_extra_covariates.csv"   # fasting_hours, townsend, blood_hour, assess_month

OUTCOME <- "taste_2w_strict"
N_JOBS  <- 10   # parallel workers

FACTOR_VARS <- c(
  "sex", "smoking", "drink", "assess_centre_id", "diabetes", "hypertension",
  "taste_drug_user_baseline_numeric", "surg_taste_affecting_baseline", "surg_taste_affecting_full",
  "any_oral_problem_baseline", "any_oral_problem", "periodontal_indicator_baseline",
  "periodontal_indicator", "neuro_disease_baseline", "head_neck_cancer_baseline",
  "psychiatric_baseline", "depression_sr", "covid_ever_total")

# Age at baseline, sex, PC1-10, smoking, alcohol, assessment center, fasting time, Townsend,
# time of blood draw, month of assessment
BASE_COV <- c("age_baseline", "sex", paste0("PC", 1:10), "smoking", "drink", "assess_centre_id",
              "fasting_hours", "townsend", "blood_hour", "assess_month")

# Nightingale biomarkers with CE-IVD certification
CE_37 <- c(
  "Total_C", "VLDL_C", "LDL_C", "HDL_C", "Total_TG", "ApoB", "ApoA1", "ApoB_by_ApoA1",
  "Total_FA", "Omega_3", "Omega_6", "PUFA", "MUFA", "SFA", "DHA",
  "Omega_3_pct", "Omega_6_pct", "PUFA_pct", "MUFA_pct", "SFA_pct", "DHA_pct",
  "PUFA_by_MUFA", "Omega_6_by_Omega_3", "Ala", "Gly", "His", "Ile", "Leu", "Val",
  "Phe", "Tyr", "Total_BCAA", "Glucose", "Lactate", "Creatinine", "Albumin", "GlycA")

# Writes <out_dir>/<name>_<number of covariates>cov.csv
run_model <- function(name, covariates, out_dir, nmr_file, pheno_files)
  mwas_glm(nmr_file, pheno_files, OUTCOME, covariates, FACTOR_VARS,
           file.path(out_dir, sprintf("%s_%dcov.csv", name, length(covariates))), N_JOBS)
