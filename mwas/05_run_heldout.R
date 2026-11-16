# Primary model (and the primary model plus BMI, diabetes and hypertension) refitted in the
# held-out cohorts: Group 2 (other White) and Group 3 (non-White).
source("mwas/mwas_config.R")
setwd(Sys.getenv("PROJECT_DIR", "."))

PRIMARY_COV <- c(BASE_COV, "surg_taste_affecting_full")
OUT_DIR     <- "output/mwas/batch3_validation"

# V1 and V3 are copied to input/assoc_results/mwas_group2_primary.csv and mwas_group3_primary.csv
run_model("V1_group2_primary", PRIMARY_COV, OUT_DIR,
          NMR_FILES[["group2"]], c(PHENO_FILES[["group2"]], PHENO_SUPP))
run_model("V2_group2_metabolic", c(PRIMARY_COV, "BMI", "diabetes", "hypertension"), OUT_DIR,
          NMR_FILES[["group2"]], c(PHENO_FILES[["group2"]], PHENO_SUPP))
run_model("V3_group3_primary", PRIMARY_COV, OUT_DIR,
          NMR_FILES[["group3"]], c(PHENO_FILES[["group3"]], PHENO_SUPP))
