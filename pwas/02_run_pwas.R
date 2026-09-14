#!/usr/bin/env Rscript
#> in   input/analysis_ready/proteomics_group1.csv
#> ship output/pwas/results/primary.csv -> input/assoc_results/pwas_primary.csv
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
