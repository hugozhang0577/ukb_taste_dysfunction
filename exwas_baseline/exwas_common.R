# Settings shared by the baseline and follow-up ExWAS runners: cohort phenotype files,
# outcome, covariate sets, filtering thresholds, and the call into the regression engine.
library(data.table)
ENGINE <- "exwas_baseline/exwas_regression.R"
source(ENGINE)

PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."))
INPUT_DIR   <- file.path(PROJECT_DIR, "input")
OUTPUT_DIR  <- file.path(PROJECT_DIR, "output")
READY_DIR   <- file.path(INPUT_DIR, "analysis_ready")

# group1 = discovery; group2, group3 = held-out cohorts
PHENO_FILES <- setNames(
  file.path(READY_DIR, sprintf("phenotype_%s.csv", c("group1", "group2", "group3"))),
  c("1", "2", "3"))

OUTCOME <- "taste_2w_strict"
EID_COL <- "eid"

COVARIATES_BASE    <- "age_baseline,sex,townsend,smoking,drink,assess_centre_id"
COVARIATES_BMI     <- paste0(COVARIATES_BASE, ",BMI")
COVARIATES_SMELL   <- paste0(COVARIATES_BASE, ",smell_any")
COVARIATES_MINIMAL <- "age_baseline,sex"
COVARIATES_FEMALE  <- "age_baseline,townsend,smoking,drink,assess_centre_id"

GROUP3_ETHNICITY_VAR <- "ethnicity_code"   # added to every model in group3
FACTOR_VARS <- "sex,smoking,drink,assess_centre_id,smell_any,ethnicity_code"

MIN_CASES_DEFAULT  <- 100     # minimum cases with complete data per exposure
CELL_DROP_DEFAULT  <- 10      # binary exposures with a smaller 2x2 cell are not tested
CELL_FIRTH_DEFAULT <- 20      # binary exposures with a smaller 2x2 cell are fitted with Firth
METHOD_DEFAULT     <- "hybrid"

SEX_FEMALE_CODE <- 0          # UK Biobank field 31: 0 = female, 1 = male

add_ethnicity <- function(covars, group_num)
  if (group_num == "3") paste0(covars, ",", GROUP3_ETHNICITY_VAR) else covars

prep_sex_subset <- function(pheno_in, pheno_out, sex_code = SEX_FEMALE_CODE) {
  p <- fread(pheno_in)
  fwrite(p[sex == sex_code], pheno_out)
  invisible(pheno_out)
}

# Results are written to `output`; log_file and append_summary() are kept for the runners' call interface.
call_exwas <- function(omics, pheno, covars, metadata, output, log_file,
                       method = METHOD_DEFAULT, min_cases = MIN_CASES_DEFAULT,
                       cell_drop = CELL_DROP_DEFAULT, cell_firth = CELL_FIRTH_DEFAULT)
  run_exwas(omics, pheno, OUTCOME, covars, FACTOR_VARS, metadata, output,
            method = method, min_cases = min_cases, cell_drop = cell_drop,
            cell_firth = cell_firth, eid_col = EID_COL)

append_summary <- function(log_lines, label, report_file)
  cat(label, "\n", file = report_file, append = TRUE)
