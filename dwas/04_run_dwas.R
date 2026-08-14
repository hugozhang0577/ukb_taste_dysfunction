#!/usr/bin/env Rscript
# =============================================================================
# Disease phenome-wide association scan — batch runner
# =============================================================================
#
# PheCodes are treated as exposures and tested against the taste outcome. The
# discovery cohort is run under five covariate models; both held-out
# cross-population cohorts are run under the smell-adjusted model (see below).
#
#   primary        age, sex, assessment centre, Townsend, smoking, alcohol, BMI,
#                  and the participant's total distinct-diagnosis count
#   primary+smell  primary plus current smell change. Bounds how much of the
#                  disease signal is chemosensory rather than taste-specific
#   no_bmi         primary without BMI, since BMI may sit on the causal path
#                  between several of these diagnoses and taste change
#   no_centre      drops assessment centre and the diagnosis count together
#   minimal        age, sex and the diagnosis count only
#
# The total diagnosis count is load-bearing rather than cosmetic. Participants
# with more recorded contact with health services accumulate more PheCodes for
# reasons unrelated to taste, which inflates associations across the whole scan.
# Adjusting for it removes that ascertainment gradient.
#
# Multiple testing: Benjamini-Hochberg and Bonferroni across the full analysable
# PheCode set, computed inside dwas_regression.R.
#
# The held-out cohorts use a lower minimum-case threshold: they are much smaller
# than the discovery cohort, so the discovery threshold would leave almost
# nothing testable.
#
# Usage (locally, or as a RAP job):
#   CODE_DIR=. PROJECT_DIR=/path/to/project \
#     nohup Rscript 04_run_dwas.R > dwas_run.log 2>&1 &
# =============================================================================

suppressPackageStartupMessages(library(data.table))

CODE_DIR    <- Sys.getenv("CODE_DIR", unset = ".")
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)

ENGINE    <- file.path(CODE_DIR, "dwas_regression.R")
SUMMARISE <- file.path(CODE_DIR, "summarise_results.R")
if (!file.exists(ENGINE))
  stop("regression engine not found: ", ENGINE,
       "\n  set CODE_DIR to the directory containing dwas_regression.R")

# PheCode matrices from 02; phenotype files carrying total_dx_count from 03.
# PHECODE_DIR can point straight at 02's output directory instead of copying the
# matrices into input/analysis_ready/.
READY_DIR   <- file.path(PROJECT_DIR, "input", "analysis_ready")
PHECODE_DIR <- Sys.getenv("PHECODE_DIR", unset = READY_DIR)
RESULTS_DIR <- file.path(PROJECT_DIR, "output", "dwas", "results")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

phecode_file <- function(g) file.path(PHECODE_DIR, sprintf("phecode_matrix_%s.rds", g))
pheno_file   <- function(g) file.path(READY_DIR, sprintf("phenotype_%s.csv", g))

OUTCOME              <- "taste_2w_strict"
N_JOBS               <- as.integer(Sys.getenv("N_JOBS", unset = "10"))
MIN_CASES_DISCOVERY  <- 100
MIN_CASES_HELDOUT    <- 50
FACTOR_VARS          <- "sex,assess_centre_id,smoking,drink,smell_any"

MODELS <- data.table(
  name = c("primary+smell", "primary", "no_bmi", "no_centre", "minimal"),
  covariates = c(
    "age_baseline,sex,assess_centre_id,townsend,smoking,drink,BMI,total_dx_count,smell_any",
    "age_baseline,sex,assess_centre_id,townsend,smoking,drink,BMI,total_dx_count",
    "age_baseline,sex,assess_centre_id,townsend,smoking,drink,total_dx_count",
    "age_baseline,sex,townsend,smoking,drink,BMI",
    "age_baseline,sex,total_dx_count"
  ),
  description = c(
    "primary covariates plus current smell change",
    "age, sex, centre, Townsend, smoking, alcohol, BMI, diagnosis count",
    "primary without BMI",
    "primary without assessment centre and without the diagnosis count",
    "age, sex and the diagnosis count only"
  )
)

run_single_model <- function(name, covariates, description,
                             phecode, pheno, min_cases, output_file) {
  n_cov <- length(strsplit(covariates, ",", fixed = TRUE)[[1]])
  cat(strrep("-", 60), "\n>>> ", name, "\n", sep = "")
  cat(sprintf("    covariates (%d): %s\n", n_cov, covariates))
  cat(sprintf("    %s | min_cases = %d\n", description, min_cases))
  cat(strrep("-", 60), "\n", sep = "")
  status <- system2("Rscript", c(
    shQuote(ENGINE),
    "--phecode",     shQuote(phecode),
    "--phenotype",   shQuote(pheno),
    "--outcome",     shQuote(OUTCOME),
    "--covariates",  shQuote(covariates),
    "--factor-vars", shQuote(FACTOR_VARS),
    "--min-cases",   min_cases,
    "--output",      shQuote(output_file),
    "--n-jobs",      N_JOBS
  ))
  if (status != 0) stop("model ", name, " failed")
  invisible(status)
}

# ---- discovery cohort: all models -------------------------------------------
cat(sprintf("=== discovery cohort (%d models, min_cases = %d) ===\n",
            nrow(MODELS), MIN_CASES_DISCOVERY))
for (i in seq_len(nrow(MODELS))) {
  run_single_model(MODELS$name[i], MODELS$covariates[i], MODELS$description[i],
                   phecode_file("group1"), pheno_file("group1"), MIN_CASES_DISCOVERY,
                   file.path(RESULTS_DIR, sprintf("group1_%s.csv", MODELS$name[i])))
}

# ---- held-out cohorts --------------------------------------------------------
# Both held-out cohorts are run under the smell-adjusted covariate set, not the
# discovery primary set. The held-out cohorts exist to test whether a disease
# association is specific to taste, and smell change is the one variable that
# most readily produces a taste-like report without a taste-specific mechanism;
# adjusting for it makes the held-out estimate the conservative one.
heldout_cov <- MODELS[name == "primary+smell", covariates]
for (g in c("group2", "group3")) {
  if (!file.exists(phecode_file(g))) {
    cat("[skip]", g, "- PheCode matrix not found\n"); next
  }
  cat(sprintf("\n=== %s held-out validation (min_cases = %d) ===\n", g, MIN_CASES_HELDOUT))
  run_single_model("primary+smell", heldout_cov,
                   sprintf("%s held-out validation, smell-adjusted", g),
                   phecode_file(g), pheno_file(g), MIN_CASES_HELDOUT,
                   file.path(RESULTS_DIR, sprintf("%s_primary.csv", g)))
}

# ---- comparison table -------------------------------------------------------
system2("Rscript", c(shQuote(SUMMARISE), shQuote(RESULTS_DIR)))
# Export and dx upload to RAP  (per-model result CSVs + the comparison table)

cat("[DONE] results in:", RESULTS_DIR, "\n")
