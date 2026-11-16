# Sensitivity models in the discovery cohort (Group 1): the 19 base covariates plus oral health,
# metabolic, neurological, cancer, psychiatric, medication or surgery terms; minimal and full models.
source("mwas/mwas_config.R")
setwd(Sys.getenv("PROJECT_DIR", "."))

ORAL_ITEMS <- paste0("oral_", c("bleeding_gums", "painful_gums", "loose_teeth", "toothache",
                                "mouth_ulcers", "dentures"), "_baseline")

MODELS <- list(
  O1a_oral_any_baseline      = c(BASE_COV, "any_oral_problem_baseline"),
  O1b_oral_any_cumulative    = c(BASE_COV, "any_oral_problem"),
  O2a_periodontal_baseline   = c(BASE_COV, "periodontal_indicator_baseline"),
  O2b_periodontal_cumulative = c(BASE_COV, "periodontal_indicator"),
  O3_oral_all_items_baseline = c(BASE_COV, ORAL_ITEMS),
  M1_bmi_only                = c(BASE_COV, "BMI"),
  M2_diabetes_only           = c(BASE_COV, "diabetes"),
  M3_hypertension_only       = c(BASE_COV, "hypertension"),
  M4_metabolic_all           = c(BASE_COV, "BMI", "diabetes", "hypertension"),
  D19_neuro_baseline         = c(BASE_COV, "neuro_disease_baseline"),
  D22_cancer_baseline        = c(BASE_COV, "head_neck_cancer_baseline"),
  D24_neuro_cancer_baseline  = c(BASE_COV, "neuro_disease_baseline", "head_neck_cancer_baseline"),
  P1_psychiatric_any         = c(BASE_COV, "psychiatric_baseline"),
  P4_depression_only         = c(BASE_COV, "depression_sr"),
  S1_minimal_clean           = c("age_baseline", "sex"),
  S2_full_baseline           = c(BASE_COV, "taste_drug_user_baseline_numeric",
                                 "surg_taste_affecting_baseline", "any_oral_problem_baseline",
                                 "BMI", "diabetes", "hypertension", "neuro_disease_baseline",
                                 "head_neck_cancer_baseline", "psychiatric_baseline",
                                 "covid_ever_total"),
  N4a_drug_baseline          = c(BASE_COV, "taste_drug_user_baseline_numeric"),
  N4b_surgery_baseline       = c(BASE_COV, "surg_taste_affecting_baseline")
)

for (m in names(MODELS))
  run_model(m, MODELS[[m]], "output/mwas/batch2_sensitivity",
            NMR_FILES[["group1"]], c(PHENO_FILES[["group1"]], PHENO_SUPP))
