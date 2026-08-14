#!/usr/bin/env Rscript
# =============================================================================
# Proteome-wide association scan — primary model
# =============================================================================
#
# The proteome-wide scan of the discovery cohort under the pre-specified
# covariate set. The covariates are fixed in advance from the study design, not
# chosen by comparing candidate models:
#
#   age at baseline, sex                design
#   Townsend deprivation index          socioeconomic position
#   Olink batch                         assay-side technical variation
#   PC1-PC10                            population structure
#   baseline smoking status, alcohol-risk category
#   cumulative taste-affecting surgery  the exposure-side confounder the
#                                       protocol specifies for a taste outcome
#
# The metabolomic primary model is the same set with the assay-specific terms
# swapped: no Olink batch (the NMR assay has no equivalent), plus assessment
# centre and three pre-analytical covariates the NMR measurements need. The
# shared core makes the two omics results directly comparable.
#
# BMI and the wider metabolic block are deliberately absent: they lie on the
# causal path between the protein panel and the outcome rather than confounding
# it, so adjusting for them would remove signal by construction.
#
# Every protein is tested by Firth penalised-likelihood logistic regression, not
# as a fallback but for all proteins: with ~17k assayed participants and few
# cases per protein, maximum likelihood is biased and can fail to converge.
# Benjamini-Hochberg and Bonferroni are applied across all proteins tested.
#
# The held-out cohorts use this same engine and covariate set with their own
# phenotype inputs.
#
# Usage:
#   CODE_DIR=. PROJECT_DIR=/path/to/project \
#     nohup Rscript 02_run_pwas.R > pwas_run.log 2>&1 &
# =============================================================================

CODE_DIR    <- Sys.getenv("CODE_DIR", unset = ".")
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)

ENGINE <- file.path(CODE_DIR, "pwas_firth_analysis.R")
if (!file.exists(ENGINE))
  stop("regression engine not found: ", ENGINE,
       "\n  set CODE_DIR to the directory containing pwas_firth_analysis.R")

READY_DIR   <- file.path(PROJECT_DIR, "input", "analysis_ready")
RESULTS_DIR <- file.path(PROJECT_DIR, "output", "pwas", "results")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# QC'd protein matrix from 01. Either copy it into input/analysis_ready/ or
# point OLINK_DIR at 01's output directory.
OLINK_DIR  <- Sys.getenv("OLINK_DIR", unset = READY_DIR)
OLINK_FILE <- file.path(OLINK_DIR, "proteomics_group1.csv")   # eid + protein NPX columns
if (!file.exists(OLINK_FILE))
  stop("QC'd protein matrix not found: ", OLINK_FILE,
       "\n  run 01_olink_qc.R first, then copy its output here or set OLINK_DIR")
PHENO_FILE <- paste(c(file.path(READY_DIR, "phenotype_group1.csv"),
                      file.path(READY_DIR, "extra_covariates.csv")), collapse = ",")

OUTCOME     <- "taste_2w_strict"
N_JOBS      <- as.integer(Sys.getenv("N_JOBS", unset = "16"))
FACTOR_VARS <- "sex,olink_batch_number,smoking,drink"

# --- the pre-specified covariate set -----------------------------------------
PCS         <- paste0("PC", 1:10, collapse = ",")
PRIMARY_COV <- paste0("age_baseline,sex,olink_batch_number,", PCS,
                      ",smoking,drink,surg_taste_affecting_full,townsend")
stopifnot(length(strsplit(PRIMARY_COV, ",", fixed = TRUE)[[1]]) == 17L)

OUTPUT_FILE <- file.path(RESULTS_DIR, "primary.csv")

cat(strrep("-", 60), "\n>>> primary model\n    covariates (",
    length(strsplit(PRIMARY_COV, ",", fixed = TRUE)[[1]]), "): ", PRIMARY_COV, "\n",
    strrep("-", 60), "\n", sep = "")

status <- system2("Rscript", c(
  shQuote(ENGINE),
  "--olink",       shQuote(OLINK_FILE),
  "--phenotype",   shQuote(PHENO_FILE),
  "--outcome",     shQuote(OUTCOME),
  "--covariates",  shQuote(PRIMARY_COV),
  "--factor-vars", shQuote(FACTOR_VARS),
  "--eid-col",     shQuote("eid"),
  "--output",      shQuote(OUTPUT_FILE),
  "--n-jobs",      N_JOBS
))
if (status != 0) stop("the primary model failed")
# Export and dx upload to RAP  (the result table and its summary)

cat("[DONE] results in:", RESULTS_DIR, "\n")
