#!/usr/bin/env Rscript
# =============================================================================
# MWAS shared configuration and batch runner
# =============================================================================
#
# Sourced by every 03-05_run_batch*.R script:
#   source(file.path(dirname(sys.frame(1)$ofile), "mwas_config.R"))
#
# Holds the paths, the analysis-wide parameters, and the two helpers that the
# batch scripts use: run_model() for a single covariate model, and
# summarise_batch() for the per-batch comparison table.
#
# All I/O is rooted at $PROJECT_DIR; set it before running (defaults to ".").
# =============================================================================

suppressPackageStartupMessages(library(data.table))

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR))
  stop("PROJECT_DIR does not exist: ", PROJECT_DIR,
       "\n  set it to the project root, e.g. Sys.setenv(PROJECT_DIR = '/mnt/project')")

# --- inputs ------------------------------------------------------------------
# QC'd, preprocessed NMR matrices (output of 02_preprocess_cohorts.R)
NMR_FILES <- c(
  group1 = file.path(PROJECT_DIR, "input", "analysis_ready", "metabolomics_group1.csv"),
  group2 = file.path(PROJECT_DIR, "input", "analysis_ready", "metabolomics_group2.csv"),
  group3 = file.path(PROJECT_DIR, "input", "analysis_ready", "metabolomics_group3.csv")
)

# Phenotype/covariate files, shared with the PWAS pipeline. PHENO_SUPP carries
# the pre-analytical covariates that only the NMR models need (fasting hours,
# blood-draw hour, assessment month).
PHENO_FILES <- c(
  group1 = file.path(PROJECT_DIR, "input", "analysis_ready", "phenotype_group1.csv"),
  group2 = file.path(PROJECT_DIR, "input", "analysis_ready", "phenotype_group2.csv"),
  group3 = file.path(PROJECT_DIR, "input", "analysis_ready", "phenotype_group3.csv")
)
PHENO_SUPP <- file.path(PROJECT_DIR, "input", "analysis_ready", "nmr_extra_covariates.csv")

# --- engine ------------------------------------------------------------------
# The per-feature regression engine is shared with the PWAS pipeline; MWAS
# selects standard maximum likelihood with --method glm. With ~91k samples and
# >3k cases Firth penalisation is unnecessary, and there is no glm->Firth
# fallback: the estimator is fixed for every feature.
CODE_DIR <- Sys.getenv("CODE_DIR", unset = ".")   # directory holding these scripts
ENGINE   <- file.path(CODE_DIR, "mwas_glm_analysis.R")
if (!file.exists(ENGINE))
  stop("regression engine not found: ", ENGINE,
       "\n  set CODE_DIR to the directory containing mwas_glm_analysis.R")
METHOD  <- "glm"
N_JOBS  <- 10

# --- analysis-wide parameters ------------------------------------------------
OUTCOME <- "taste_2w_strict"
EID_COL <- "eid"

# NMR values are already standardised per SD in 02_preprocess_cohorts.R, so the
# engine is called with its default --standardize FALSE and the effect sizes are
# odds ratios per SD.

# Categorical covariates. Relative to the PWAS pipeline, olink_batch_number is
# absent: the NMR assay has no equivalent batch variable.
FACTOR_VARS <- paste(c(
  "sex", "smoking", "drink", "assess_centre_id", "diabetes", "hypertension",
  "taste_drug_user_baseline_numeric", "taste_drug_user_numeric",
  "surg_taste_affecting_baseline", "surg_taste_affecting_full",
  "any_oral_problem_baseline", "any_oral_problem",
  "periodontal_indicator_baseline", "periodontal_indicator", "ckd_binary",
  "neuro_disease_baseline", "neuro_disease_ever", "parkinsons_baseline",
  "ms_baseline", "epilepsy_baseline", "head_neck_cancer_baseline",
  "head_neck_cancer_ever", "psychiatric_baseline", "psychiatric_severe",
  "psychiatric_common", "depression_sr", "anxiety_sr", "covid_ever_total"
), collapse = ",")

# --- results root ------------------------------------------------------------
MWAS_RESULTS_DIR <- file.path(PROJECT_DIR, "output", "mwas")

# =============================================================================
# Run one covariate model
# =============================================================================
# The output file name carries the covariate count so that models are
# distinguishable at a glance and the summary step can parse them back out.

run_model <- function(model_name, covariates, description, output_dir,
                      nmr_file, pheno_file, index, n_models) {
  n_cov       <- length(strsplit(covariates, ",", fixed = TRUE)[[1]])
  output_file <- file.path(output_dir, sprintf("%s_%dcov.csv", model_name, n_cov))

  cat("\n", strrep("-", 72), "\n", sep = "")
  cat(sprintf(" Model %d/%d: %s\n", index, n_models, model_name))
  cat(sprintf("   covariates (%d): %s\n", n_cov, covariates))
  cat(sprintf("   description    : %s\n", description))
  cat(sprintf("   output         : %s\n", output_file))
  cat(strrep("-", 72), "\n", sep = "")

  t0 <- Sys.time()
  status <- system2("Rscript", c(
    shQuote(ENGINE),
    "--olink",        shQuote(nmr_file),
    "--phenotype",    shQuote(pheno_file),
    "--outcome",      shQuote(OUTCOME),
    "--covariates",   shQuote(covariates),
    "--factor-vars",  shQuote(FACTOR_VARS),
    "--eid-col",      shQuote(EID_COL),
    "--output",       shQuote(output_file),
    "--method",       shQuote(METHOD),
    "--n-jobs",       N_JOBS
  ))
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  result  <- if (status == 0) "SUCCESS" else "FAILED"
  if (status != 0) cat("[ERROR] model", model_name, "failed\n")

  # Prepend the model specification to the engine's own summary, so each
  # result file documents the covariate set it was produced under.
  summary_file <- sub("\\.csv$", "_summary.txt", output_file)
  if (status == 0 && file.exists(summary_file)) {
    writeLines(c(
      strrep("=", 60),
      "Pipeline: MWAS (Nightingale NMR metabolomics)",
      paste("Model:", model_name),
      paste("Description:", description),
      sprintf("Covariates (%d): %s", n_cov, covariates),
      paste("Factor variables:", FACTOR_VARS),
      strrep("=", 60), "",
      readLines(summary_file)
    ), summary_file)
  }

  cat(sprintf("Status: %s | %.1f min\n", result, elapsed))
  cat(sprintf("[%s] %s: %.1f min\n", result, model_name, elapsed),
      file = file.path(output_dir, "batch_run_log.txt"), append = TRUE)
  invisible(status)
}

# =============================================================================
# Per-batch comparison table
# =============================================================================
# One row per model run: convergence, hit counts at each threshold, the
# inflation factor, and the strongest metabolite. The inflation factor is a run
# diagnostic only; it is not reported for the scan, where an outcome-anchored
# panel makes a high value expected rather than informative.

summarise_batch <- function(output_dir, batch_name) {
  files <- list.files(output_dir, pattern = ".*cov\\.csv$", full.names = TRUE)
  if (length(files) == 0) { cat("No result files found in", output_dir, "\n"); return(invisible(NULL)) }

  results <- rbindlist(lapply(files, function(f) {
    dt <- tryCatch(fread(f), error = function(e) NULL)
    if (is.null(dt)) return(NULL)
    fname   <- basename(f)
    valid_p <- dt$pval[!is.na(dt$pval)]
    lambda  <- if (length(valid_p) > 0)
      median(qchisq(1 - valid_p, df = 1)) / qchisq(0.5, df = 1) else NA_real_
    top_row <- dt[which.min(pval)]
    data.table(
      Model          = gsub("_[0-9]+cov\\.csv$", "", fname),
      N_Covariates   = as.integer(gsub(".*_([0-9]+)cov\\.csv$", "\\1", fname)),
      N_Metabolites  = nrow(dt),
      N_Converged    = sum(dt$converged, na.rm = TRUE),
      Converge_Rate  = round(mean(dt$converged, na.rm = TRUE) * 100, 1),
      N_Nominal      = sum(dt$pval < 0.05, na.rm = TRUE),
      N_Suggestive   = sum(dt$pval < 0.001, na.rm = TRUE),
      N_FDR          = sum(dt$pval_fdr < 0.05, na.rm = TRUE),
      N_Bonferroni   = sum(dt$pval_bonf < 0.05, na.rm = TRUE),
      Lambda         = round(lambda, 3),
      Top_Metabolite = top_row$protein,
      Top_OR         = round(top_row$or, 3),
      Top_P          = format(top_row$pval, scientific = TRUE, digits = 2)
    )
  }), fill = TRUE)

  setorder(results, Model)
  cat("\n", strrep("=", 90), "\n  ", batch_name, " - summary\n", strrep("=", 90), "\n\n", sep = "")
  print(results[, .(Model, N_Covariates, Lambda, N_Nominal, N_FDR, N_Bonferroni,
                    Top_Metabolite, Top_P)])
  fwrite(results, file.path(output_dir, "batch_models_comparison.csv"))
  # Export and dx upload to RAP  (per-model result CSVs + batch_models_comparison.csv)
  cat("\nwritten:", file.path(output_dir, "batch_models_comparison.csv"), "\n")
  invisible(results)
}

# =============================================================================
# Batch header / footer
# =============================================================================

start_batch <- function(batch_name, output_dir, n_models, extra = character(0)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  cat("\n", strrep("=", 72), "\n  ", batch_name, "\n", strrep("=", 72), "\n", sep = "")
  cat("  started : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
  cat("  models  : ", n_models, "\n", sep = "")
  cat("  workers : ", N_JOBS, "\n", sep = "")
  for (e in extra) cat("  ", e, "\n", sep = "")
  cat("  output  : ", output_dir, "\n", sep = "")
  cat(strrep("=", 72), "\n", sep = "")
  writeLines(c(paste(batch_name, "- run log"),
               paste("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
               paste("Models:", n_models)),
             file.path(output_dir, "batch_run_log.txt"))
  Sys.time()
}

end_batch <- function(batch_name, t0, output_dir) {
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat("\n", strrep("=", 72), "\n", sep = "")
  cat("  ", batch_name, " complete\n", sep = "")
  cat(sprintf("  total elapsed: %.0fh %.0fm\n", elapsed %/% 60, elapsed %% 60))
  cat("  results: ", output_dir, "\n", sep = "")
  cat(strrep("=", 72), "\n", sep = "")
  cat(sprintf("Completed: %s\nTotal: %.0fh %.0fm\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), elapsed %/% 60, elapsed %% 60),
      file = file.path(output_dir, "batch_run_log.txt"), append = TRUE)
}

cat("[CONFIG] MWAS configuration loaded\n")
cat("  NMR (discovery): ", NMR_FILES[["group1"]], "\n", sep = "")
cat("  phenotype      : ", PHENO_FILES[["group1"]], "\n", sep = "")
cat("  engine         : ", ENGINE, " (--method ", METHOD, ")\n", sep = "")
cat("  workers        : ", N_JOBS, "\n", sep = "")
