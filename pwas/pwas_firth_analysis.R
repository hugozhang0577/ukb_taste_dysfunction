#!/usr/bin/env Rscript
# =============================================================================
# pwas_firth_analysis.R
# Proteome-wide association analysis (Olink) by Firth penalised-likelihood
# logistic regression — one model per protein, parallelised, with BH-FDR and
# Bonferroni correction.
#
# Generic engine: outcome, covariates, factor variables, and I/O paths are all
# supplied on the command line; no study-specific configuration here.
#
# Usage (see 02_run_pwas.R for the model matrix):
#   Rscript pwas_firth_analysis.R \
#     --olink      $READY_DIR/proteomics_group1.csv \
#     --phenotype  $PHENO_MAIN,$PHENO_SUPP \
#     --outcome    taste_2w_strict \
#     --covariates "age_baseline,sex,olink_batch_number,PC1,...,PC10,smoking,drink,surg_taste_affecting_full,townsend" \
#     --factor-vars "sex,olink_batch_number,smoking,drink" \
#     --output     $RESULTS_DIR/primary.csv \
#     --n-jobs     16
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(logistf)
  library(parallel)
})

# =============================================================================
# Command-line interface
# =============================================================================
option_list <- list(
  make_option(c("-o", "--olink"),       type = "character", help = "Olink protein matrix (.parquet/.csv/.rds)"),
  make_option(c("-p", "--phenotype"),   type = "character", help = "Phenotype file(s), comma-separated for multi-file merge"),
  make_option(c("-y", "--outcome"),     type = "character", help = "Binary outcome column (0/1)"),
  make_option(c("-c", "--covariates"),  type = "character", help = "Comma-separated covariate list"),
  make_option(c("-f", "--factor-vars"), type = "character", default = NULL, help = "Comma-separated factor covariates"),
  make_option(c("--output"),            type = "character", help = "Output results CSV path"),
  make_option(c("--eid-col"),           type = "character", default = "eid", help = "Sample ID column [default %default]"),
  make_option(c("--protein-cols"),      type = "character", default = NULL, help = "Explicit protein columns (default: auto-detect)"),
  make_option(c("-j", "--n-jobs"),      type = "integer", default = 1, help = "Parallel cores [default %default]"),
  make_option(c("--chunk-size"),        type = "integer", default = 100, help = "Proteins per parallel chunk [default %default]"),
  make_option(c("--standardize"),       action = "store_true", default = FALSE, help = "Z-score the protein values"),
  make_option(c("--resume"),            action = "store_true", default = FALSE, help = "Resume from an existing output")
)
args <- parse_args(OptionParser(option_list = option_list,
                                description = "PWAS Firth logistic regression"))

# =============================================================================
# Data loading (parquet / csv / rds / tsv)
# =============================================================================
load_data <- function(filepath) {
  ext <- tools::file_ext(filepath)
  if (ext == "parquet") {
    if (!requireNamespace("arrow", quietly = TRUE)) stop("arrow needed for parquet")
    as.data.table(arrow::read_parquet(filepath))
  } else if (ext == "csv") {
    fread(filepath)
  } else if (ext == "rds") {
    as.data.table(readRDS(filepath))
  } else if (ext %in% c("tsv", "txt")) {
    fread(filepath, sep = "\t")
  } else stop(paste("unsupported file format:", ext))
}

# =============================================================================
# Single-protein Firth regression
#   skip if n_total < 50 or n_case < 5; otherwise outcome ~ protein + covariates
# =============================================================================
run_firth_single <- function(protein_name, data, outcome, covariates) {
  vars_needed <- c(protein_name, outcome, covariates)
  analysis_data <- na.omit(data[, ..vars_needed])
  n_total <- nrow(analysis_data)
  n_case  <- sum(analysis_data[[outcome]])

  if (n_total < 50 || n_case < 5) {
    return(data.table(protein = protein_name, n_total = n_total, n_case = n_case,
                      beta = NA_real_, se = NA_real_, or = NA_real_,
                      or_lower = NA_real_, or_upper = NA_real_, z = NA_real_,
                      pval = NA_real_, converged = FALSE,
                      error = "Insufficient samples or cases"))
  }

  formula_obj <- as.formula(paste0(outcome, " ~ ", protein_name, " + ",
                                   paste(covariates, collapse = " + ")))
  tryCatch({
    fit <- logistf(formula_obj, data = analysis_data,
                   control = logistf.control(maxit = 250))
    idx <- which(names(fit$coefficients) == protein_name)
    if (length(idx) == 0) idx <- 2
    beta <- fit$coefficients[idx]
    se   <- sqrt(diag(vcov(fit)))[idx]
    data.table(
      protein = protein_name, n_total = n_total, n_case = n_case,
      beta = beta, se = se,
      or = exp(beta),
      or_lower = exp(fit$ci.lower[idx]),   # profile-likelihood CI
      or_upper = exp(fit$ci.upper[idx]),
      z = beta / se,
      pval = fit$prob[idx],
      converged = TRUE, error = NA_character_
    )
  }, error = function(e) {
    data.table(protein = protein_name, n_total = n_total, n_case = n_case,
               beta = NA_real_, se = NA_real_, or = NA_real_,
               or_lower = NA_real_, or_upper = NA_real_, z = NA_real_,
               pval = NA_real_, converged = FALSE, error = as.character(e$message))
  })
}

# =============================================================================
# Main pipeline
# =============================================================================
main <- function() {
  for (req in c("olink", "phenotype", "outcome", "covariates", "output"))
    if (is.null(args[[req]])) stop("missing required arg: --", req)

  covariates <- trimws(strsplit(args$covariates, ",")[[1]])
  cat("[INFO] covariates:", paste(covariates, collapse = ", "), "\n")

  # ---- Olink matrix + protein-column detection ----
  olink_dt <- load_data(args$olink)
  cat(sprintf("[INFO] Olink: %d x %d\n", nrow(olink_dt), ncol(olink_dt)))
  if (!is.null(args$`protein-cols`)) {
    protein_list <- trimws(strsplit(args$`protein-cols`, ",")[[1]])
  } else {
    exclude_patterns <- c(args$`eid-col`, "batch", "plate", "well", "qc", "index")
    protein_list <- setdiff(names(olink_dt),
                            names(olink_dt)[grepl(paste(exclude_patterns, collapse = "|"),
                                                  names(olink_dt), ignore.case = TRUE)])
    cat(sprintf("[INFO] auto-detected %d protein columns\n", length(protein_list)))
  }

  # ---- phenotype (inner-join multiple files) ----
  pheno_files <- trimws(strsplit(args$phenotype, ",")[[1]])
  pheno_list <- lapply(pheno_files, load_data)
  pheno_dt <- if (length(pheno_list) == 1) pheno_list[[1]] else
    Reduce(function(x, y) merge(x, y, by = args$`eid-col`, all = FALSE), pheno_list)

  needed_cols  <- c(args$`eid-col`, args$outcome, covariates)
  analysis_dt  <- merge(olink_dt, pheno_dt[, ..needed_cols], by = args$`eid-col`)
  cat(sprintf("[INFO] merged sample size: %d\n", nrow(analysis_dt)))
  rm(olink_dt, pheno_dt, pheno_list); gc()

  # ---- outcome to 0/1 ----
  analysis_dt[[args$outcome]] <- as.numeric(analysis_dt[[args$outcome]])
  analysis_dt <- analysis_dt[analysis_dt[[args$outcome]] %in% c(0, 1), ]
  cat(sprintf("[INFO] cases: %d, controls: %d\n",
              sum(analysis_dt[[args$outcome]]), sum(1 - analysis_dt[[args$outcome]])))

  # ---- factor coercion (user-specified, integer64, or character) ----
  user_factor_vars <- if (!is.null(args$`factor-vars`))
    trimws(strsplit(args$`factor-vars`, ",")[[1]]) else character(0)
  for (cov in covariates) {
    if (cov %in% user_factor_vars || inherits(analysis_dt[[cov]], "integer64") ||
        is.character(analysis_dt[[cov]]) || is.factor(analysis_dt[[cov]])) {
      analysis_dt[[cov]] <- as.factor(as.character(analysis_dt[[cov]]))
    }
  }

  if (args$standardize) {
    for (p in protein_list) if (p %in% names(analysis_dt))
      analysis_dt[[p]] <- scale(analysis_dt[[p]])[, 1]
  }

  protein_list <- protein_list[protein_list %in% names(analysis_dt)]
  cat(sprintf("[INFO] proteins to analyse: %d\n", length(protein_list)))

  # ---- resume support ----
  completed <- character(0)
  if (args$resume && file.exists(args$output)) completed <- fread(args$output)$protein
  proteins_to_run <- setdiff(protein_list, completed)
  if (length(proteins_to_run) == 0) { cat("[INFO] nothing to run\n"); return(invisible()) }

  # ---- run (parallel chunks or serial) ----
  if (args$`n-jobs` > 1) {
    cat(sprintf("[INFO] %d cores\n", args$`n-jobs`))
    n_chunks <- ceiling(length(proteins_to_run) / args$`chunk-size`)
    cl <- makeCluster(args$`n-jobs`)
    clusterExport(cl, c("analysis_dt", "args", "covariates", "run_firth_single"),
                  envir = environment())
    clusterEvalQ(cl, { library(data.table); library(logistf) })
    results_list <- list()
    for (ci in seq_len(n_chunks)) {
      cs <- (ci - 1) * args$`chunk-size` + 1
      ce <- min(ci * args$`chunk-size`, length(proteins_to_run))
      results_list <- c(results_list, parLapply(cl, proteins_to_run[cs:ce], function(p)
        run_firth_single(p, analysis_dt, args$outcome, covariates)))
      cat(sprintf("  %d / %d proteins\n", ce, length(proteins_to_run)))
    }
    stopCluster(cl)
  } else {
    results_list <- lapply(seq_along(proteins_to_run), function(i) {
      if (i %% 100 == 0) cat(sprintf("  %d / %d\n", i, length(proteins_to_run)))
      run_firth_single(proteins_to_run[i], analysis_dt, args$outcome, covariates)
    })
  }

  results_dt <- rbindlist(results_list)
  if (length(completed) > 0) results_dt <- rbind(fread(args$output), results_dt)

  # ---- multiple testing ----
  ok <- !is.na(results_dt$pval)
  results_dt$pval_fdr <- NA_real_
  results_dt$pval_fdr[ok] <- p.adjust(results_dt$pval[ok], method = "BH")
  results_dt$pval_bonf <- pmin(results_dt$pval * sum(ok), 1)
  results_dt[, sig_nominal    := pval < 0.05]
  results_dt[, sig_suggestive := pval < 0.001]
  results_dt[, sig_fdr        := pval_fdr < 0.05]
  results_dt[, sig_bonf       := pval_bonf < 0.05]
  setorder(results_dt, pval)

  dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
  fwrite(results_dt, args$output)
  # Export and dx upload to RAP  (this model's result table and its summary)

  # ---- summary ----
  vp <- results_dt$pval[ok]
  lambda_gc <- median(qchisq(1 - vp, df = 1)) / qchisq(0.5, df = 1)
  cat(sprintf("\n[SUMMARY] proteins=%d converged=%.1f%% medN=%d medCases=%d lambda=%.3f FDR=%d Bonf=%d\n",
              nrow(results_dt), 100 * mean(results_dt$converged),
              as.integer(median(results_dt$n_total, na.rm = TRUE)),
              as.integer(median(results_dt$n_case, na.rm = TRUE)),
              lambda_gc, sum(results_dt$sig_fdr, na.rm = TRUE),
              sum(results_dt$sig_bonf, na.rm = TRUE)))

  summary_path <- sub("\\.csv$", "_summary.txt", args$output)
  writeLines(c(
    "PWAS Firth logistic regression summary",
    paste("Date:", as.character(Sys.time())),
    paste("Proteins analysed:", nrow(results_dt)),
    paste("Lambda:", round(lambda_gc, 3)),
    paste("Nominal significant:", sum(results_dt$sig_nominal, na.rm = TRUE)),
    paste("FDR significant:", sum(results_dt$sig_fdr, na.rm = TRUE)),
    paste("Bonferroni significant:", sum(results_dt$sig_bonf, na.rm = TRUE))
  ), summary_path)
  cat("[INFO] saved:", args$output, "\n")
}

main()
