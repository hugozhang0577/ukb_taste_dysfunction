#!/usr/bin/env Rscript
# =============================================================================
# exwas_regression.R
# Generic hybrid glm / Firth exposure-wide regression engine with grouped
# Benjamini-Hochberg FDR control.
#
# Design:
#   - Pure algorithm + I/O; contains NO study-specific configuration.
#   - All study-specific parameters (outcome, covariates, factor variables,
#     FDR grouping, filter thresholds) are supplied via the command line and a
#     variable-metadata CSV.
#   - Supports three fit modes: 'glm', 'firth', 'hybrid'.
#   - FDR is controlled by Benjamini-Hochberg INDEPENDENTLY within each level
#     of metadata$fdr_group (i.e. one family per pre-specified source).
#
# This engine is invoked by 02_run_exwas.R, which assembles its arguments
# from exwas_common.R. It is identical to the engine used for the
# baseline and follow-up exposure-wide scans and for the disease phenome-wide
# scan, differing only in the metadata and covariates passed in.
#
# Reporting discipline: every filter step prints "N before -> N after" style
# counts and a per-variable exclusion log is written alongside the results.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(logistf)
})

# -----------------------------------------------------------------------------
# Command-line interface
# -----------------------------------------------------------------------------
option_list <- list(
  make_option("--omics",       type="character",
              help="Exposure matrix CSV (eid + one column per exposure)"),
  make_option("--phenotype",   type="character",
              help="Phenotype CSV path(s); comma-separated for multi-file merge"),
  make_option("--outcome",     type="character",
              help="Outcome column name (binary 0/1)"),
  make_option("--covariates",  type="character",
              help="Covariate column names, comma-separated"),
  make_option("--factor-vars", type="character", default="", dest="factor_vars",
              help="Columns to coerce to factor, comma-separated [default none]"),
  make_option("--variable-metadata", type="character", dest="var_meta",
              help="Metadata CSV; required columns: variable, var_type, fdr_group"),
  make_option("--min-cases",            type="integer", default=100, dest="min_cases",
              help="L1 filter: minimum outcome-complete cases per variable [default %default]"),
  make_option("--cell-drop-threshold",  type="integer", default=10,  dest="cell_drop",
              help="L2 filter: min 2x2 cell for a binary exposure; below this = DROP [default %default]"),
  make_option("--cell-firth-threshold", type="integer", default=20,  dest="cell_firth",
              help="L3 trigger: min 2x2 cell below which Firth is used in hybrid mode [default %default]"),
  make_option("--method", type="character", default="hybrid",
              help="Fit method: 'glm' | 'firth' | 'hybrid' [default %default]"),
  make_option("--eid-col", type="character", default="eid", dest="eid_col",
              help="Sample ID column [default %default]"),
  make_option("--output", type="character",
              help="Output results CSV path (filter log written alongside)")
)
opt <- parse_args(OptionParser(option_list=option_list))

for (f in c("omics","phenotype","outcome","covariates","var_meta","output"))
  if (is.null(opt[[f]])) stop(sprintf("missing required arg: --%s",
                                      gsub("_","-",f)))
stopifnot(opt$cell_drop <= opt$cell_firth,
          opt$method %in% c("glm","firth","hybrid"))

COV <- strsplit(opt$covariates, ",", fixed=TRUE)[[1]]
FAC <- if (nzchar(opt$factor_vars))
  strsplit(opt$factor_vars, ",", fixed=TRUE)[[1]] else character(0)

# -----------------------------------------------------------------------------
# Load
# -----------------------------------------------------------------------------
t0 <- Sys.time()

pheno_files <- strsplit(opt$phenotype, ",", fixed=TRUE)[[1]]
pheno <- fread(pheno_files[1])
for (pf in pheno_files[-1]) {
  p2 <- fread(pf)
  keep <- setdiff(names(p2), setdiff(names(pheno), opt$eid_col))
  pheno <- merge(pheno, p2[, ..keep], by=opt$eid_col, all.x=TRUE)
}
omics <- fread(opt$omics)
meta  <- fread(opt$var_meta)

req_meta <- c("variable","var_type","fdr_group")
if (!all(req_meta %in% names(meta)))
  stop(sprintf("metadata missing required columns: %s",
               paste(setdiff(req_meta, names(meta)), collapse=",")))

# --- Drop phenotype columns that overlap declared exposures ---
# Some exposures (e.g. oral-health, treatment-history, eGFR items) were
# originally stored in shared phenotype files and later migrated into the
# exposure pool. Without this drop, merge() creates .x/.y suffixes and the bare
# variable name in metadata fails to match, causing a silent
# 'not_in_omics_data' exclusion. Exposure-side values take precedence because
# the metadata describes the exposure-side cleaning.
overlap_cols <- intersect(names(pheno), meta$variable)
overlap_cols <- setdiff(overlap_cols, opt$eid_col)
if (length(overlap_cols) > 0) {
  cat(sprintf("[MERGE] dropping %d pheno columns overlapping exposures: %s\n",
              length(overlap_cols), paste(overlap_cols, collapse=",")))
  pheno[, (overlap_cols) := NULL]
}

dat <- merge(pheno, omics, by=opt$eid_col, all=FALSE)

missing_cov <- setdiff(COV, names(dat))
if (length(missing_cov) > 0)
  stop(sprintf("missing covariates in merged data: %s",
               paste(missing_cov, collapse=",")))
if (!(opt$outcome %in% names(dat)))
  stop(sprintf("outcome '%s' not in merged data", opt$outcome))

for (cv in intersect(FAC, names(dat))) dat[[cv]] <- factor(dat[[cv]])
dat[[".y"]] <- as.integer(dat[[opt$outcome]])

N_total <- nrow(dat)
N_cases <- sum(dat$.y == 1, na.rm=TRUE)
cat(sprintf("[LOAD] merged n=%d cases=%d covariates=%d vars_in_meta=%d\n",
            N_total, N_cases, length(COV), nrow(meta)))

# -----------------------------------------------------------------------------
# Hybrid fit for a single exposure
#   - hybrid: glm first; fall back to Firth on a pre-emptive sparse-cell
#     trigger (min 2x2 cell < cell_firth), non-convergence, or a separation
#     warning ("fitted probabilities numerically 0 or 1").
#   - firth : always Firth penalised likelihood (profile-likelihood CIs).
#   - glm   : standard maximum-likelihood only, never falls back.
# -----------------------------------------------------------------------------
fit_one <- function(sub, v, covars, var_type, method, cell_firth) {
  fml <- as.formula(paste(".y ~", paste(c(v, covars), collapse=" + ")))
  out <- list(beta=NA_real_, se=NA_real_, OR=NA_real_,
              OR_lower=NA_real_, OR_upper=NA_real_, p_value=NA_real_,
              method_used=NA_character_, converged=FALSE,
              error_msg=NA_character_)

  is_binary <- var_type %in% c("binary","derived_binary")
  force_firth <- (method == "firth")

  # pre-emptive Firth trigger (hybrid mode, binary exposures only)
  if (method == "hybrid" && is_binary) {
    tab <- table(sub$.y, sub[[v]])
    if (nrow(tab) < 2 || ncol(tab) < 2) {
      out$error_msg <- "degenerate_2x2"; return(out)
    }
    if (min(tab) < cell_firth) force_firth <- TRUE
  }

  # ---- glm path ----
  sep_warn <- FALSE; glm_err <- NULL; fit <- NULL
  if (!force_firth) {
    fit <- tryCatch(
      withCallingHandlers(
        glm(fml, data=sub, family=binomial()),
        warning = function(w) {
          if (grepl("fitted probabilities numerically 0 or 1",
                    conditionMessage(w), fixed=TRUE))
            sep_warn <<- TRUE
          invokeRestart("muffleWarning")
        }),
      error = function(e) { glm_err <<- conditionMessage(e); NULL }
    )
  }

  need_firth <- force_firth || is.null(fit) ||
                !isTRUE(fit$converged) || isTRUE(sep_warn)
  if (method == "glm") need_firth <- FALSE  # glm-only mode never falls back

  if (!need_firth && !is.null(fit)) {
    co <- tryCatch(summary(fit)$coefficients, error=function(e) NULL)
    if (is.null(co) || !(v %in% rownames(co))) {
      out$error_msg <- "glm_coef_missing"; return(out)
    }
    b  <- co[v, "Estimate"]
    se <- co[v, "Std. Error"]
    p  <- co[v, "Pr(>|z|)"]
    out$beta <- b; out$se <- se; out$p_value <- p
    out$OR   <- exp(b)
    out$OR_lower <- exp(b - 1.96*se)
    out$OR_upper <- exp(b + 1.96*se)
    out$method_used <- "glm"
    out$converged   <- TRUE
    return(out)
  }
  if (method == "glm") {
    out$error_msg <- if (is.null(fit)) paste0("glm_failed:", glm_err)
                     else "glm_nonconverged"
    out$method_used <- "glm_failed"
    return(out)
  }

  # ---- Firth path ----
  fit2 <- tryCatch(
    logistf::logistf(fml, data=sub,
                     control=logistf::logistf.control(maxit=100)),
    error = function(e) {
      out$error_msg <<- paste0("firth_error:", conditionMessage(e)); NULL }
  )
  if (is.null(fit2)) { out$method_used <- "firth_failed"; return(out) }
  idx <- which(names(fit2$coefficients) == v)
  if (length(idx) == 0) {
    out$error_msg <- "firth_coef_missing"; return(out)
  }
  b  <- unname(fit2$coefficients[idx])
  se <- unname(sqrt(diag(vcov(fit2)))[idx])
  out$beta      <- b
  out$se        <- se
  out$p_value   <- unname(fit2$prob[idx])
  out$OR        <- exp(b)
  out$OR_lower  <- exp(unname(fit2$ci.lower[idx]))    # profile-likelihood CI
  out$OR_upper  <- exp(unname(fit2$ci.upper[idx]))
  out$method_used <- if (method == "firth") "firth" else "firth_fallback"
  out$converged   <- TRUE
  return(out)
}

# -----------------------------------------------------------------------------
# Filter + fit loop
#   L1: minimum outcome-complete cases (--min-cases)
#   variance guard: drop zero-variance exposures
#   L2: sparse-cell DROP for binary exposures (min 2x2 cell < --cell-drop)
#   L3: sparse-cell Firth trigger handled inside fit_one (--cell-firth)
# -----------------------------------------------------------------------------
log_rows <- list()
log_add  <- function(vars, reason) {
  if (length(vars) == 0) return(invisible())
  log_rows[[length(log_rows)+1]] <<- data.table(variable=vars, reason=reason)
}

meta_in_data <- meta[variable %in% names(dat)]
log_add(setdiff(meta$variable, names(dat)), "not_in_omics_data")

res_rows <- vector("list", nrow(meta_in_data))
n_l1 <- 0L; n_l2 <- 0L; n_var0 <- 0L

for (i in seq_len(nrow(meta_in_data))) {
  v  <- meta_in_data$variable[i]
  vt <- meta_in_data$var_type[i]
  fg <- meta_in_data$fdr_group[i]

  cols <- c(".y", v, COV)
  sub  <- dat[, ..cols]
  sub  <- sub[complete.cases(sub)]
  n_c  <- nrow(sub)
  n_ca <- sum(sub$.y == 1)
  vr   <- if (n_c >= 2) suppressWarnings(var(as.numeric(sub[[v]]))) else 0

  # L1: outcome-complete cases
  if (n_ca < opt$min_cases) {
    log_add(v, sprintf("L1_n_cases_complete<%d", opt$min_cases))
    n_l1 <- n_l1 + 1L; next
  }
  # variance guard
  if (is.na(vr) || vr <= 0) {
    log_add(v, "zero_variance"); n_var0 <- n_var0 + 1L; next
  }
  # L2: sparse-cell drop (binary only)
  if (vt %in% c("binary","derived_binary")) {
    tab <- table(sub$.y, sub[[v]])
    if (nrow(tab) < 2 || ncol(tab) < 2 || min(tab) < opt$cell_drop) {
      log_add(v, sprintf("L2_min_cell<%d", opt$cell_drop))
      n_l2 <- n_l2 + 1L; next
    }
  }

  r <- fit_one(sub, v, COV, vt, opt$method, opt$cell_firth)
  res_rows[[i]] <- data.table(
    variable=v, var_type=vt, fdr_group=fg,
    n_complete=n_c, n_cases_complete=n_ca,
    beta=r$beta, se=r$se,
    OR=r$OR, OR_lower=r$OR_lower, OR_upper=r$OR_upper,
    p_value=r$p_value, method_used=r$method_used,
    converged=r$converged,
    error_msg=if (is.null(r$error_msg)) NA_character_ else r$error_msg)

  if (i %% 25 == 0)
    cat(sprintf("  ... %d/%d\n", i, nrow(meta_in_data)))
}
results <- rbindlist(res_rows, fill=TRUE)

# -----------------------------------------------------------------------------
# Benjamini-Hochberg FDR, independently within each fdr_group
# -----------------------------------------------------------------------------
results[, p_adj := NA_real_]
for (g in unique(results$fdr_group)) {
  idx <- which(results$fdr_group == g & !is.na(results$p_value))
  if (length(idx) > 0)
    results$p_adj[idx] <- p.adjust(results$p_value[idx], method="BH")
}
results[, significant := !is.na(p_adj) & p_adj < 0.05]

# -----------------------------------------------------------------------------
# Genomic inflation factor (1-df chi-square), reported globally and per source
# -----------------------------------------------------------------------------
compute_lambda <- function(p) {
  p <- p[!is.na(p) & p > 0 & p < 1]
  if (length(p) < 5) return(NA_real_)
  chisq <- qchisq(1 - p, df = 1)
  median(chisq) / qchisq(0.5, df = 1)
}

# -----------------------------------------------------------------------------
# Write results + per-variable exclusion log
# -----------------------------------------------------------------------------
dir.create(dirname(opt$output), showWarnings=FALSE, recursive=TRUE)
fwrite(results, opt$output)
log_path <- sub("\\.csv$", ".filter_log.csv", opt$output)
fwrite(rbindlist(log_rows, fill=TRUE), log_path)

# -----------------------------------------------------------------------------
# Structured stdout summary (harvested by the batch runner via grep "^SUMMARY")
# -----------------------------------------------------------------------------
elapsed <- as.numeric(difftime(Sys.time(), t0, units="secs"))
cat("\n========== SUMMARY ==========\n")
cat(sprintf("SUMMARY n_samples: %d\n", N_total))
cat(sprintf("SUMMARY n_cases: %d\n", N_cases))
cat(sprintf("SUMMARY vars_in_metadata: %d\n", nrow(meta)))
cat(sprintf("SUMMARY vars_in_data: %d\n", nrow(meta_in_data)))
cat(sprintf("SUMMARY excluded_L1: %d\n", n_l1))
cat(sprintf("SUMMARY excluded_zero_variance: %d\n", n_var0))
cat(sprintf("SUMMARY excluded_L2: %d\n", n_l2))
n_fit <- sum(!is.na(results$method_used))
cat(sprintf("SUMMARY models_fit: %d\n", n_fit))
cat(sprintf("SUMMARY method_glm: %d\n",
            sum(results$method_used=="glm", na.rm=TRUE)))
cat(sprintf("SUMMARY method_firth: %d\n",
            sum(results$method_used=="firth", na.rm=TRUE)))
cat(sprintf("SUMMARY method_firth_fallback: %d\n",
            sum(results$method_used=="firth_fallback", na.rm=TRUE)))
cat(sprintf("SUMMARY errors: %d\n", sum(!is.na(results$error_msg))))
for (g in sort(unique(results$fdr_group))) {
  sub <- results[fdr_group == g]
  cat(sprintf("SUMMARY sig_%s: %d/%d\n", g,
              sum(sub$significant, na.rm=TRUE), nrow(sub)))
}
# Global inflation factor (all fitted models in this cohort)
lambda_global <- compute_lambda(results$p_value)
cat(sprintf("SUMMARY lambda_global: %.3f\n", lambda_global))

# Per-source inflation factor
for (g in sort(unique(results$fdr_group))) {
  lg <- compute_lambda(results[fdr_group == g, p_value])
  cat(sprintf("SUMMARY lambda_%s: %.3f\n", g,
              if (is.na(lg)) NaN else lg))
}

top <- results[!is.na(p_adj)][order(p_adj)][1:min(10, .N)]
for (j in seq_len(nrow(top))) {
  cat(sprintf("SUMMARY top%d: %s | %s | OR=%.3f (%.3f-%.3f) | p=%.2e | q=%.2e | %s\n",
              j, top$fdr_group[j], top$variable[j],
              top$OR[j], top$OR_lower[j], top$OR_upper[j],
              top$p_value[j], top$p_adj[j], top$method_used[j]))
}
cat(sprintf("SUMMARY elapsed_sec: %.1f\n", elapsed))
cat(sprintf("SUMMARY output_results: %s\n", opt$output))
cat(sprintf("SUMMARY output_filter_log: %s\n", log_path))
cat("=============================\n")
