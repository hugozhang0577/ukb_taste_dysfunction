#!/usr/bin/env Rscript
#> ship output/mwas/primary/primary_20cov.csv -> input/assoc_results/mwas_primary.csv
source(file.path(Sys.getenv("CODE_DIR", unset = "."), "mwas_config.R"))

NMR_FILE   <- NMR_FILES[["group1"]]
PHENO_FILE <- paste(c(PHENO_FILES[["group1"]], PHENO_SUPP), collapse = ",")
OUTPUT_DIR <- file.path(MWAS_RESULTS_DIR, "primary")
BATCH_NAME <- "MWAS primary model (discovery cohort)"

# --- the pre-specified covariate set -----------------------------------------
DESIGN_PREANALYTIC <- paste("age_baseline", "sex", "assess_centre_id",
                            "fasting_hours", "townsend", "blood_hour",
                            "assess_month", sep = ",")
ANCESTRY   <- paste0("PC", 1:10, collapse = ",")
LIFESTYLE  <- "smoking,drink"
TASTE_SURG <- "surg_taste_affecting_full"

PRIMARY_COV <- paste(DESIGN_PREANALYTIC, ANCESTRY, LIFESTYLE, TASTE_SURG, sep = ",")
stopifnot(length(strsplit(PRIMARY_COV, ",", fixed = TRUE)[[1]]) == 20L)

t0 <- start_batch(BATCH_NAME, OUTPUT_DIR, 1L,
                  extra = c(paste("NMR      :", NMR_FILE),
                            paste("phenotype:", PHENO_FILE)))

run_model("primary", PRIMARY_COV,
          "pre-specified covariate set; BMI and the metabolic block excluded as mediators",
          OUTPUT_DIR, NMR_FILE, PHENO_FILE, 1L, 1L)

summarise_batch(OUTPUT_DIR, BATCH_NAME)
end_batch(BATCH_NAME, t0, OUTPUT_DIR)
