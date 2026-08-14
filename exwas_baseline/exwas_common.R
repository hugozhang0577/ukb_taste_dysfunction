#!/usr/bin/env Rscript
# =============================================================================
# Shared configuration and helpers for the exposure-wide association scans
# =============================================================================
#
# Sourced by 02_run_exwas.R and by the follow-up scan's runner. Both scans use
# the same regression engine (exwas_regression.R) and the same covariate
# strategy; they differ only in which exposures they test, which dictionary
# defines the FDR families, and which run variants they sweep.
#
# Covariate strategy, shared with the protein-, metabolite- and disease-wide
# scans: a six-variable primary set, plus an ethnicity term in the ancestrally
# heterogeneous cohort. BMI and smell change are deliberately NOT primary
# covariates for questionnaire exposures, because for many of them BMI is a
# consequence rather than a confounder and smell change is part of what is being
# measured; both are added back as explicit sensitivity variants instead.
#
# All paths resolve under $PROJECT_DIR; nothing here is machine-specific.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
CODE_DIR <- Sys.getenv("CODE_DIR", unset = ".")

INPUT_DIR  <- file.path(PROJECT_DIR, "input")
OUTPUT_DIR <- file.path(PROJECT_DIR, "output")
READY_DIR  <- file.path(INPUT_DIR, "analysis_ready")

# Per-cohort phenotype files, shared with the other scans.
# group1 = discovery; group2 and group3 = held-out cross-population cohorts.
PHENO_FILES <- setNames(
  file.path(READY_DIR, sprintf("phenotype_%s.csv", c("group1", "group2", "group3"))),
  c("1", "2", "3"))

# --- analysis parameters -----------------------------------------------------
OUTCOME <- "taste_2w_strict"
EID_COL <- "eid"

COVARIATES_BASE     <- "age_baseline,sex,townsend,smoking,drink,assess_centre_id"
COVARIATES_BMI      <- paste0(COVARIATES_BASE, ",BMI")
COVARIATES_SMELL    <- paste0(COVARIATES_BASE, ",smell_any")
COVARIATES_MINIMAL  <- "age_baseline,sex"
# Female-only sub-batch: sex is dropped because the data are subset to one sex
# and the term would be single-level.
COVARIATES_FEMALE   <- "age_baseline,townsend,smoking,drink,assess_centre_id"

GROUP3_ETHNICITY_VAR <- "ethnicity_code"
FACTOR_VARS <- "sex,smoking,drink,assess_centre_id,smell_any,ethnicity_code"

MIN_CASES_DEFAULT  <- 100   # minimum outcome-complete cases for a variable to be tested
CELL_DROP_DEFAULT  <- 10    # a binary exposure with a smaller 2x2 cell is dropped
CELL_FIRTH_DEFAULT <- 20    # below this, the fit switches to Firth
METHOD_DEFAULT     <- "hybrid"

# UK Biobank Field 31 coding: 0 = female, 1 = male
SEX_FEMALE_CODE <- 0

ENGINE <- file.path(CODE_DIR, "exwas_regression.R")
if (!file.exists(ENGINE)) ENGINE <- file.path(CODE_DIR, "..", "exwas_baseline", "exwas_regression.R")
if (!file.exists(ENGINE))
  stop("regression engine not found; set CODE_DIR to the directory holding ",
       "exwas_regression.R (or its sibling exwas_baseline/)")

# =============================================================================
# Helpers
# =============================================================================

# The ancestrally heterogeneous cohort carries an extra ethnicity term, since
# residual population structure there would otherwise be absorbed by the
# exposures.
add_ethnicity <- function(covars, group_num)
  if (group_num == "3") paste0(covars, ",", GROUP3_ETHNICITY_VAR) else covars

# Subset a phenotype file to one sex. The engine's inner join then drops the
# other sex on the exposure side automatically.
prep_sex_subset <- function(pheno_in, pheno_out, sex_code = SEX_FEMALE_CODE) {
  p <- fread(pheno_in)
  if (!"sex" %in% names(p)) stop("no sex column in phenotype file: ", pheno_in)
  n0 <- nrow(p)
  p_sub <- p[sex == sex_code]
  fwrite(p_sub, pheno_out)
  cat(sprintf("[PREP] %s: %d -> %d rows (sex == %d)\n",
              basename(pheno_out), n0, nrow(p_sub), sex_code))
  invisible(pheno_out)
}

# Run one variant of the scan and tee the engine's log to disk.
call_exwas <- function(omics, pheno, covars, metadata, output, log_file,
                       method = METHOD_DEFAULT, min_cases = MIN_CASES_DEFAULT,
                       cell_drop = CELL_DROP_DEFAULT, cell_firth = CELL_FIRTH_DEFAULT) {
  out <- system2("Rscript", c(
    shQuote(ENGINE),
    "--omics",                 shQuote(omics),
    "--phenotype",             shQuote(pheno),
    "--outcome",               shQuote(OUTCOME),
    "--covariates",            shQuote(covars),
    "--factor-vars",           shQuote(FACTOR_VARS),
    "--variable-metadata",     shQuote(metadata),
    "--eid-col",               shQuote(EID_COL),
    "--output",                shQuote(output),
    "--method",                shQuote(method),
    "--min-cases",             min_cases,
    "--cell-drop-threshold",   cell_drop,
    "--cell-firth-threshold",  cell_firth
  ), stdout = TRUE, stderr = TRUE)
  writeLines(out, log_file)
  cat(out, sep = "\n")
  invisible(out)
}

# The engine prints one SUMMARY line per FDR family; collect them into the
# run-level report so the whole sweep can be read in one place.
append_summary <- function(log_lines, label, report_file) {
  cat(c("", strrep("=", 62), paste(">>>", label), strrep("=", 62),
        sub("^SUMMARY ", "", grep("^SUMMARY ", log_lines, value = TRUE))),
      sep = "\n", file = report_file, append = TRUE)
}
