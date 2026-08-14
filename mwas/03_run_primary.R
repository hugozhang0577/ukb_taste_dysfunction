#!/usr/bin/env Rscript
# =============================================================================
# MWAS — primary model
# =============================================================================
#
# The metabolome-wide scan of the discovery cohort under the pre-specified
# covariate set. The covariates are fixed in advance from the study design, not
# chosen by comparing candidate models:
#
#   age at baseline, sex, assessment centre
#   PC1-PC10                              population structure
#   baseline smoking status, alcohol-risk category
#   Townsend deprivation index
#   cumulative taste-affecting surgery    the exposure-side confounder the
#                                         protocol specifies for a taste outcome
#   fasting hours, blood-draw hour, assessment month
#                                         pre-analytical, specific to the NMR assay
#
# BMI and the wider metabolic block are deliberately absent. They lie on the
# causal path between the metabolite panel and the outcome rather than
# confounding it, so adjusting for them would remove signal by construction.
# Their effect is quantified in the sensitivity analyses (04), not here.
#
# The same set is used for the held-out cohorts in 05, so discovery and
# validation estimates are directly comparable.
#
# Runtime is long (roughly 15-30 h on ~28 cores: 327 measures x ~91k
# participants), so run it detached:
#   CODE_DIR=. PROJECT_DIR=/path/to/project \
#     nohup Rscript 03_run_primary.R > mwas_primary.log 2>&1 &
# =============================================================================

source(file.path(Sys.getenv("CODE_DIR", unset = "."), "mwas_config.R"))

NMR_FILE   <- NMR_FILES[["group1"]]
PHENO_FILE <- paste(c(PHENO_FILES[["group1"]], PHENO_SUPP), collapse = ",")
OUTPUT_DIR <- file.path(MWAS_RESULTS_DIR, "primary")
BATCH_NAME <- "MWAS primary model (discovery cohort)"

# --- the pre-specified covariate set -----------------------------------------
# Design and pre-analytical terms first, then ancestry, lifestyle, and the
# taste-specific surgery term. Term order does not affect the metabolite
# coefficient; it is kept as written so the string matches the one the reported
# results were produced under.
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
