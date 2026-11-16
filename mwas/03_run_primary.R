# Primary metabolome-wide scan in the discovery cohort (Group 1): 20 pre-specified covariates;
# BMI and the metabolic block are left out as potential mediators.
source("mwas/mwas_config.R")
setwd(Sys.getenv("PROJECT_DIR", "."))

PRIMARY_COV <- c("age_baseline", "sex", "assess_centre_id", "fasting_hours", "townsend",
                 "blood_hour", "assess_month", paste0("PC", 1:10), "smoking", "drink",
                 "surg_taste_affecting_full")

# primary_20cov.csv is copied to input/assoc_results/mwas_primary.csv for the figures
run_model("primary", PRIMARY_COV, "output/mwas/primary",
          NMR_FILES[["group1"]], c(PHENO_FILES[["group1"]], PHENO_SUPP))
