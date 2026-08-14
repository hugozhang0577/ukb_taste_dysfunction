#!/usr/bin/env Rscript
#' =============================================================================
#' PWAS / MWAS per-feature logistic regression engine
#' =============================================================================
#'
#' Shared engine for the UK Biobank Olink proteome-wide (PWAS) and Nightingale
#' NMR metabolome-wide (MWAS) association scans. One feature (protein or
#' metabolite) is tested against a binary outcome at a time.
#'
#' Two interchangeable estimators are selected with --method:
#'   - firth : logistf penalised-likelihood logistic regression (PWAS; small N,
#'             rare events). NOT a fallback — it is used for every feature.
#'   - glm   : standard maximum-likelihood logistic regression (MWAS; large N,
#'             ~91k samples). No Firth penalisation is applied.
#' The two methods are mutually exclusive; there is no automatic glm->firth
#' fallback. Multiple testing: Benjamini-Hochberg + Bonferroni across all
#' tested features.
#'
#' Parallelised over features (chunked makeCluster); memory-aware.
#' All I/O is via command-line arguments — no hard-coded paths.
#' =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(parallel)
  library(arrow)  # for parquet support
})

# =============================================================================
# Command-line arguments
# =============================================================================

option_list <- list(
  make_option(c("-o", "--olink"), type = "character", 
              help = "Feature matrix: proteins or metabolites (.parquet, .csv, .rds)"),
  
  make_option(c("-p", "--phenotype"), type = "character",
              help = "Phenotype/covariate file(s), comma-separated"),
  
  make_option(c("-y", "--outcome"), type = "character",
              help = "Outcome column name (binary 0/1)"),
  
  make_option(c("-c", "--covariates"), type = "character",
              help = "Covariate column names, comma-separated"),
  
  make_option(c("-f", "--factor-vars"), type = "character", default = NULL,
              help = "Variables to treat as factors, comma-separated (e.g. sex,Batch,centre)"),
  
  make_option(c("--output"), type = "character",
              help = "Output file path"),
  
  make_option(c("--eid-col"), type = "character", default = "eid",
              help = "Sample-ID column [default: %default]"),
  
  make_option(c("--protein-cols"), type = "character", default = NULL,
              help = "Feature columns to analyse, comma-separated; auto-detected by default"),
  
  make_option(c("-j", "--n-jobs"), type = "integer", default = 1,
              help = "Number of parallel workers [default: %default]"),
  
  make_option(c("--chunk-size"), type = "integer", default = 100,
              help = "Features per chunk [default: %default]"),
  
  make_option(c("--standardize"), action = "store_true", default = FALSE,
              help = "Z-score standardising the feature values"),
  
  make_option(c("--resume"), action = "store_true", default = FALSE,
              help = "Resume from an existing partial output"),
  
  make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
              help = "Verbose logging"),
  
  make_option(c("-m", "--method"), type = "character", default = "firth",
              help = "Estimator: 'firth' (logistf, small samples) or 'glm' (standard logistic, faster on large samples) [default: %default]")
)

parser <- OptionParser(
  usage = "%prog [options]",
  option_list = option_list,
  description = "PWAS Firth Logistic Regression Analysis"
)

args <- parse_args(parser)

# Validate --method and load the matching package
args$method <- tolower(args$method)
if (!args$method %in% c("firth", "glm")) {
  stop("--method must be either 'firth' or 'glm'")
}
if (args$method == "firth") {
  suppressPackageStartupMessages(library(logistf))
}

# =============================================================================
# Data loading
# =============================================================================

load_data <- function(filepath) {
  ext <- tools::file_ext(filepath)
  
  if (ext == "parquet") {
    dt <- as.data.table(read_parquet(filepath))
  } else if (ext == "csv") {
    dt <- fread(filepath)
  } else if (ext == "rds") {
    dt <- as.data.table(readRDS(filepath))
  } else if (ext %in% c("tsv", "txt")) {
    dt <- fread(filepath, sep = "\t")
  } else {
    stop(paste("unsupported file format:", ext))
  }
  
  return(dt)
}

# =============================================================================
# Empty-result template shared by both estimators
# =============================================================================

empty_result <- function(protein_name, n_total, n_case, error_msg) {
  data.table(
    protein = protein_name,
    n_total = n_total,
    n_case = n_case,
    beta = NA_real_,
    se = NA_real_,
    or = NA_real_,
    or_lower = NA_real_,
    or_upper = NA_real_,
    z = NA_real_,
    pval = NA_real_,
    converged = FALSE,
    error = error_msg
  )
}

# =============================================================================
# Firth penalised-likelihood regression for one feature
# =============================================================================

run_firth_single <- function(protein_name, data, outcome, covariates) {
  
  vars_needed <- c(protein_name, outcome, covariates)
  analysis_data <- na.omit(data[, ..vars_needed])
  
  n_total <- nrow(analysis_data)
  n_case <- sum(analysis_data[[outcome]])
  
  if (n_total < 50 || n_case < 5) {
    return(empty_result(protein_name, n_total, n_case, "Insufficient samples or cases"))
  }
  
  formula_str <- paste0(outcome, " ~ ", protein_name, " + ", 
                        paste(covariates, collapse = " + "))
  formula_obj <- as.formula(formula_str)
  
  result <- tryCatch({
    fit <- logistf(formula_obj, data = analysis_data, 
                   control = logistf.control(maxit = 250))
    
    idx <- which(names(fit$coefficients) == protein_name)
    if (length(idx) == 0) idx <- 2
    
    beta <- fit$coefficients[idx]
    se <- sqrt(diag(vcov(fit)))[idx]
    
    data.table(
      protein = protein_name,
      n_total = n_total,
      n_case = n_case,
      beta = beta,
      se = se,
      or = exp(beta),
      or_lower = exp(fit$ci.lower[idx]),
      or_upper = exp(fit$ci.upper[idx]),
      z = beta / se,
      pval = fit$prob[idx],
      converged = TRUE,
      error = NA_character_
    )
    
  }, error = function(e) {
    empty_result(protein_name, n_total, n_case, as.character(e$message))
  })
  
  return(result)
}

# =============================================================================
# Standard GLM logistic regression for one feature (large samples; 20-50x faster)
# =============================================================================

run_glm_single <- function(protein_name, data, outcome, covariates) {
  
  vars_needed <- c(protein_name, outcome, covariates)
  analysis_data <- na.omit(data[, ..vars_needed])
  
  n_total <- nrow(analysis_data)
  n_case <- sum(analysis_data[[outcome]])
  
  if (n_total < 50 || n_case < 5) {
    return(empty_result(protein_name, n_total, n_case, "Insufficient samples or cases"))
  }
  
  formula_str <- paste0(outcome, " ~ ", protein_name, " + ", 
                        paste(covariates, collapse = " + "))
  formula_obj <- as.formula(formula_str)
  
  result <- tryCatch({
    fit <- glm(formula_obj, data = analysis_data, family = binomial(link = "logit"))
    
    coef_summary <- summary(fit)$coefficients
    idx <- which(rownames(coef_summary) == protein_name)
    if (length(idx) == 0) idx <- 2
    
    beta <- coef_summary[idx, "Estimate"]
    se <- coef_summary[idx, "Std. Error"]
    z <- coef_summary[idx, "z value"]
    pval <- coef_summary[idx, "Pr(>|z|)"]
    
    ci <- confint.default(fit)[protein_name, ]
    
    data.table(
      protein = protein_name,
      n_total = n_total,
      n_case = n_case,
      beta = beta,
      se = se,
      or = exp(beta),
      or_lower = exp(ci[1]),
      or_upper = exp(ci[2]),
      z = z,
      pval = pval,
      converged = fit$converged,
      error = NA_character_
    )
    
  }, error = function(e) {
    empty_result(protein_name, n_total, n_case, as.character(e$message))
  })
  
  return(result)
}

# =============================================================================
# Main analysis
# =============================================================================

main <- function() {
  
  cat("\n")
  cat("============================================================\n")
  cat("PWAS/MWAS Logistic Regression Analysis (R Version)\n")
  cat("============================================================\n")
  cat("Method:", toupper(args$method), 
      ifelse(args$method == "glm", "(standard, fast)", "(Firth penalized, robust)"), "\n\n")
  
  # Check the required arguments
  if (is.null(args$olink) || is.null(args$phenotype) || 
      is.null(args$outcome) || is.null(args$covariates) || 
      is.null(args$output)) {
    stop("missing required arguments; see --help")
  }
  
  # Parse the covariate list
  covariates <- trimws(strsplit(args$covariates, ",")[[1]])
  cat("covariates:", paste(covariates, collapse = ", "), "\n")
  
  # Load the feature matrix
  cat("[INFO] loading feature matrix:", args$olink, "\n")
  olink_dt <- load_data(args$olink)
  cat("[INFO] feature matrix dimensions:", nrow(olink_dt), "x", ncol(olink_dt), "\n")
  
  # Determine which columns are features
  if (!is.null(args$`protein-cols`)) {
    protein_list <- trimws(strsplit(args$`protein-cols`, ",")[[1]])
  } else {
    # Auto-detect: drop eid and the usual non-feature columns
    exclude_patterns <- c(args$`eid-col`, "batch", "plate", "well", "qc", "index")
    protein_list <- setdiff(names(olink_dt), 
                            names(olink_dt)[grepl(paste(exclude_patterns, collapse = "|"), 
                                                   names(olink_dt), ignore.case = TRUE)])
    cat("[INFO] auto-detected", length(protein_list), "feature columns\n")
  }
  
  # Load the phenotype/covariate data
  pheno_files <- trimws(strsplit(args$phenotype, ",")[[1]])
  cat("[INFO] loading", length(pheno_files), "phenotype file(s)\n")
  
  pheno_list <- lapply(pheno_files, function(f) {
    dt <- load_data(f)
    cat("  -", f, ":", nrow(dt), "rows,", ncol(dt), "columns\n")
    return(dt)
  })
  
  # Merge the phenotype files
  if (length(pheno_list) == 1) {
    pheno_dt <- pheno_list[[1]]
  } else {
    pheno_dt <- Reduce(function(x, y) {
      merge(x, y, by = args$`eid-col`, all = FALSE)
    }, pheno_list)
  }
  cat("[INFO] merged phenotype data:", nrow(pheno_dt), "x", ncol(pheno_dt), "\n")
  
  # Join features to phenotype
  needed_cols <- c(args$`eid-col`, args$outcome, covariates)
  pheno_subset <- pheno_dt[, ..needed_cols]
  
  analysis_dt <- merge(olink_dt, pheno_subset, by = args$`eid-col`)
  cat("[INFO] samples after join:", nrow(analysis_dt), "\n")
  
  # Release memory
  rm(olink_dt, pheno_dt, pheno_list)
  gc()
  
  # Check the outcome variable
  analysis_dt[[args$outcome]] <- as.numeric(analysis_dt[[args$outcome]])
  valid_outcome <- analysis_dt[[args$outcome]] %in% c(0, 1)
  if (sum(!valid_outcome) > 0) {
    cat("[WARNING] removed samples with a non-0/1 outcome:", sum(!valid_outcome), "\n")
    analysis_dt <- analysis_dt[valid_outcome, ]
  }
  
  cat("[INFO] cases:", sum(analysis_dt[[args$outcome]]), 
      ", controls:", sum(1 - analysis_dt[[args$outcome]]), "\n")
  
  # Handle categorical covariates
  # Parse the user-specified factor variables
  if (!is.null(args$`factor-vars`)) {
    user_factor_vars <- trimws(strsplit(args$`factor-vars`, ",")[[1]])
    cat("[INFO] user-specified factor variables:", paste(user_factor_vars, collapse = ", "), "\n")
  } else {
    user_factor_vars <- character(0)
  }
  
  for (cov in covariates) {
    # Explicitly requested as a factor
    if (cov %in% user_factor_vars) {
      cat("[INFO] converting", cov, "converted to a factor (user-specified)\n")
      analysis_dt[[cov]] <- as.factor(as.character(analysis_dt[[cov]]))
      next
    }
    
    # integer64 columns (from arrow/parquet) need explicit conversion
    if (inherits(analysis_dt[[cov]], "integer64")) {
      cat("[INFO] converting", cov, "converted from integer64 to a factor\n")
      analysis_dt[[cov]] <- as.factor(as.character(analysis_dt[[cov]]))
      next
    }
    
    # Character, or already a factor
    if (is.character(analysis_dt[[cov]]) || is.factor(analysis_dt[[cov]])) {
      cat("[INFO] converting", cov, "converted to a factor (was character/factor)\n")
      analysis_dt[[cov]] <- as.factor(analysis_dt[[cov]])
    }
  }
  
  # Standardise the feature values
  if (args$standardize) {
    cat("[INFO] Z-score standardising the feature values\n")
    for (p in protein_list) {
      if (p %in% names(analysis_dt)) {
        analysis_dt[[p]] <- scale(analysis_dt[[p]])[, 1]
      }
    }
  }
  
  # Keep only features present in the data
  protein_list <- protein_list[protein_list %in% names(analysis_dt)]
  cat("[INFO] features to analyse:", length(protein_list), "\n")
  
  # Resume check
  completed_proteins <- character(0)
  if (args$resume && file.exists(args$output)) {
    existing <- fread(args$output)
    completed_proteins <- existing$protein
    cat("[INFO] resuming; already complete:", length(completed_proteins), "features\n")
  }
  
  proteins_to_run <- setdiff(protein_list, completed_proteins)
  cat("[INFO] features remaining:", length(proteins_to_run), "\n\n")
  
  if (length(proteins_to_run) == 0) {
    cat("[INFO] all features already analysed\n")
    return(invisible(NULL))
  }
  
  # Run the analysis
  start_time <- Sys.time()
  
  # Progress display
  show_progress <- function(current, total, start_time, converged_count, width = 50) {
    percent <- current / total
    filled <- round(width * percent)
    bar <- paste0(
      "\r[",
      paste(rep("█", filled), collapse = ""),
      paste(rep("░", width - filled), collapse = ""),
      "] "
    )
    
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    if (current > 0) {
      eta_secs <- elapsed / current * (total - current)
      eta_str <- sprintf("%02d:%02d:%02d", 
                         as.integer(eta_secs / 3600),
                         as.integer((eta_secs %% 3600) / 60),
                         as.integer(eta_secs %% 60))
      elapsed_str <- sprintf("%02d:%02d:%02d",
                            as.integer(elapsed / 3600),
                            as.integer((elapsed %% 3600) / 60),
                            as.integer(elapsed %% 60))
      speed <- current / elapsed * 60  # features per minute
    } else {
      eta_str <- "--:--:--"
      elapsed_str <- "00:00:00"
      speed <- 0
    }
    
    status <- sprintf("%d/%d (%.1f%%) | elapsed: %s | remaining: %s | %.1f/min | converged: %d", 
                      current, total, percent * 100, 
                      elapsed_str, eta_str, speed, converged_count)
    
    cat(paste0(bar, status))
    flush.console()
  }
  
  # Select the estimator
  run_func <- if (args$method == "firth") run_firth_single else run_glm_single
  
  if (args$`n-jobs` > 1) {
    # Parallel execution with progress reporting
    cat("[INFO] using", args$`n-jobs`, "parallel workers\n")
    cat("[INFO] progress is reported per chunk in parallel mode\n\n")
    
    # Chunked so that progress can be reported
    chunk_size <- args$`chunk-size`
    n_proteins <- length(proteins_to_run)
    n_chunks <- ceiling(n_proteins / chunk_size)
    
    results_list <- list()
    converged_count <- 0
    
    cl <- makeCluster(args$`n-jobs`)
    clusterExport(cl, c("analysis_dt", "args", "covariates", 
                        "run_func", "run_firth_single", "run_glm_single", "empty_result"),
                  envir = environment())
    clusterEvalQ(cl, {
      library(data.table)
      if (exists("args") && args$method == "firth") library(logistf)
    })
    
    for (chunk_idx in seq_len(n_chunks)) {
      chunk_start <- (chunk_idx - 1) * chunk_size + 1
      chunk_end <- min(chunk_idx * chunk_size, n_proteins)
      chunk_proteins <- proteins_to_run[chunk_start:chunk_end]
      
      chunk_results <- parLapply(cl, chunk_proteins, function(p) {
        run_func(p, analysis_dt, args$outcome, covariates)
      })
      
      results_list <- c(results_list, chunk_results)
      converged_count <- sum(sapply(results_list, function(x) x$converged))
      
      show_progress(chunk_end, n_proteins, start_time, converged_count)
    }
    
    stopCluster(cl)
    cat("\n")
    
  } else {
    # Serial execution with a progress bar
    results_list <- list()
    n_proteins <- length(proteins_to_run)
    converged_count <- 0
    
    cat("\n")
    for (i in seq_along(proteins_to_run)) {
      p <- proteins_to_run[i]
      
      results_list[[i]] <- run_func(p, analysis_dt, args$outcome, covariates)
      
      if (results_list[[i]]$converged) {
        converged_count <- converged_count + 1
      }
      
      # Update the progress bar
      show_progress(i, n_proteins, start_time, converged_count)
      
      # Verbose mode: print each feature's result on its own line
      if (args$verbose) {
        res <- results_list[[i]]
        status <- ifelse(res$converged, "[OK]", "[FAIL]")
        cat(sprintf("\n  %s %s: OR=%.3f, P=%.2e", 
                    status, res$protein, res$or, res$pval))
      }
    }
    cat("\n\n")
  }
  
  # Combine results
  results_dt <- rbindlist(results_list)
  
  # Merge with the existing partial output when resuming
  if (length(completed_proteins) > 0) {
    existing <- fread(args$output)
    results_dt <- rbind(existing, results_dt)
  }
  
  # Benjamini-Hochberg FDR
  valid_pvals <- !is.na(results_dt$pval)
  results_dt$pval_fdr <- NA_real_
  results_dt$pval_fdr[valid_pvals] <- p.adjust(results_dt$pval[valid_pvals], method = "BH")
  
  # Bonferroni correction
  n_tests <- sum(valid_pvals)
  results_dt$pval_bonf <- pmin(results_dt$pval * n_tests, 1)
  
  # Significance flags
  results_dt[, sig_nominal := pval < 0.05]
  results_dt[, sig_suggestive := pval < 0.001]
  results_dt[, sig_fdr := pval_fdr < 0.05]
  results_dt[, sig_bonf := pval_bonf < 0.05]
  
  # Order by p-value
  setorder(results_dt, pval)
  
  # Write results
  dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
  fwrite(results_dt, args$output)
  cat("\n[INFO] results written:", args$output, "\n")
  
  # Print the summary
  cat("\n")
  cat("============================================================\n")
  cat("Analysis Summary (", toupper(args$method), ")\n")
  cat("============================================================\n")
  cat("Method:", args$method, "\n")
  cat("features tested:", nrow(results_dt), "\n")
  cat("converged:", sum(results_dt$converged), 
      sprintf("(%.1f%%)", 100 * mean(results_dt$converged)), "\n")
  cat("analytic N:", median(results_dt$n_total, na.rm = TRUE), "(median)\n")
  cat("cases:", median(results_dt$n_case, na.rm = TRUE), "(median)\n")
  cat("\nsignificant results:\n")
  cat("  nominal (P < 0.05):", sum(results_dt$sig_nominal, na.rm = TRUE), "\n")
  cat("  suggestive (P < 0.001):", sum(results_dt$sig_suggestive, na.rm = TRUE), "\n")
  cat("  FDR-significant (FDR < 0.05):", sum(results_dt$sig_fdr, na.rm = TRUE), "\n")
  cat("  Bonferroni-significant:", sum(results_dt$sig_bonf, na.rm = TRUE), "\n")
  
  # Lambda (genomic inflation factor)
  valid_p <- results_dt$pval[valid_pvals]
  chi2 <- qchisq(1 - valid_p, df = 1)
  lambda_gc <- median(chi2) / qchisq(0.5, df = 1)
  cat("\nGenomic inflation factor (λ):", round(lambda_gc, 3), "\n")
  
  # Top 10
  cat("\nTop 10 features:\n")
  top10 <- head(results_dt, 10)
  for (i in 1:nrow(top10)) {
    row <- top10[i, ]
    cat(sprintf("  %s: OR=%.3f (%.3f-%.3f), P=%.2e\n",
                row$protein, row$or, row$or_lower, row$or_upper, row$pval))
  }
  
  elapsed_total <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  cat(sprintf("\ntotal elapsed: %.1f min\n", elapsed_total))
  cat("============================================================\n")
  
  # Write the summary
  summary_path <- sub("\\.csv$", "_summary.txt", args$output)
  sink(summary_path)
  cat("PWAS/MWAS Logistic Regression Summary\n")
  cat("Date:", as.character(Sys.time()), "\n")
  cat("Method:", args$method, "\n")
  cat("Proteins analyzed:", nrow(results_dt), "\n")
  cat("Lambda:", round(lambda_gc, 3), "\n")
  cat("Nominal significant:", sum(results_dt$sig_nominal, na.rm = TRUE), "\n")
  cat("FDR significant:", sum(results_dt$sig_fdr, na.rm = TRUE), "\n")
  cat("Bonferroni significant:", sum(results_dt$sig_bonf, na.rm = TRUE), "\n")
  sink()
  cat("[INFO] summary written:", summary_path, "\n")
}

# Entry point
main()