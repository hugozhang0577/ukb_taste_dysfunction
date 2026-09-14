#!/usr/bin/env Rscript
source(file.path(Sys.getenv("CODE_DIR", unset = "."), "mwas_config.R"))

NMR_FILE   <- NMR_FILES[["group1"]]
PHENO_FILE <- paste(c(PHENO_FILES[["group1"]], PHENO_SUPP), collapse = ",")
OUTPUT_DIR <- file.path(MWAS_RESULTS_DIR, "batch2_sensitivity")
BATCH_NAME <- "MWAS batch 2: sensitivity and oral stratification"

BASE_COV <- paste0("age_baseline,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,",
                   "smoking,drink,assess_centre_id,fasting_hours,townsend,",
                   "blood_hour,assess_month")
stopifnot(length(strsplit(BASE_COV, ",", fixed = TRUE)[[1]]) == 19L)

MODELS <- data.table(
  name = c(
    "O1a_oral_any_baseline", "O1b_oral_any_cumulative",
    "O2a_periodontal_baseline", "O2b_periodontal_cumulative",
    "O3_oral_all_items_baseline",
    "M1_bmi_only", "M2_diabetes_only", "M3_hypertension_only", "M4_metabolic_all",
    "D19_neuro_baseline", "D22_cancer_baseline", "D24_neuro_cancer_baseline",
    "P1_psychiatric_any", "P4_depression_only",
    "S1_minimal_clean", "S2_full_baseline",
    "N4a_drug_baseline", "N4b_surgery_baseline"
  ),
  covariates = c(
    paste0(BASE_COV, ",any_oral_problem_baseline"),
    paste0(BASE_COV, ",any_oral_problem"),
    paste0(BASE_COV, ",periodontal_indicator_baseline"),
    paste0(BASE_COV, ",periodontal_indicator"),
    paste0(BASE_COV, ",oral_bleeding_gums_baseline,oral_painful_gums_baseline,",
                     "oral_loose_teeth_baseline,oral_toothache_baseline,",
                     "oral_mouth_ulcers_baseline,oral_dentures_baseline"),
    paste0(BASE_COV, ",BMI"),
    paste0(BASE_COV, ",diabetes"),
    paste0(BASE_COV, ",hypertension"),
    paste0(BASE_COV, ",BMI,diabetes,hypertension"),
    paste0(BASE_COV, ",neuro_disease_baseline"),
    paste0(BASE_COV, ",head_neck_cancer_baseline"),
    paste0(BASE_COV, ",neuro_disease_baseline,head_neck_cancer_baseline"),
    paste0(BASE_COV, ",psychiatric_baseline"),
    paste0(BASE_COV, ",depression_sr"),
    "age_baseline,sex",
    paste0(BASE_COV, ",taste_drug_user_baseline_numeric,surg_taste_affecting_baseline,",
                     "any_oral_problem_baseline,BMI,diabetes,hypertension,",
                     "neuro_disease_baseline,head_neck_cancer_baseline,",
                     "psychiatric_baseline,covid_ever_total"),
    paste0(BASE_COV, ",taste_drug_user_baseline_numeric"),
    paste0(BASE_COV, ",surg_taste_affecting_baseline")
  ),
  description = c(
    "O1a: any oral problem, baseline definition",
    "O1b: any oral problem, cumulative definition",
    "O2a: periodontal indicator, baseline definition",
    "O2b: periodontal indicator, cumulative definition",
    "O3: each oral item separately, baseline definitions",
    "M1: +BMI",
    "M2: +diabetes",
    "M3: +hypertension",
    "M4: +BMI, diabetes, hypertension - signal loss expected",
    "D19: +neurological disease, baseline",
    "D22: +head and neck cancer, baseline",
    "D24: +neurological disease and cancer",
    "P1: +any psychiatric history - reporting-bias check",
    "P4: +self-reported depression",
    "S1: minimal adjustment - upper bound on the effect size",
    "S2: fullest adjustment - lower bound on the effect size",
    "N4a: +taste-affecting medication, baseline",
    "N4b: +taste-affecting surgery, baseline"
  )
)

t0 <- start_batch(BATCH_NAME, OUTPUT_DIR, nrow(MODELS))

for (i in seq_len(nrow(MODELS))) {
  run_model(MODELS$name[i], MODELS$covariates[i], MODELS$description[i],
            OUTPUT_DIR, NMR_FILE, PHENO_FILE, i, nrow(MODELS))
}

summarise_batch(OUTPUT_DIR, BATCH_NAME)
end_batch(BATCH_NAME, t0, OUTPUT_DIR)
