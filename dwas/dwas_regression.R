#!/usr/bin/env Rscript
# =============================================================================
# Disease-wide association scan — regression engine
# =============================================================================
#
# Tests each PheCode as a binary exposure against the taste outcome, with glm,
# a Firth penalised-likelihood fallback, and multiple-testing correction across
# the full analysable set.
#
# Design (consistent with the protein-, metabolite- and exposure-wide engines):
#   - The PheCode matrix (RDS) and the phenotype file are kept separate and
#     merged at run time by eid; no pre-merged subject-level file is stored.
#   - Covariates are supplied via --covariates, so the same engine runs the
#     primary model and every sensitivity variant.
#   - One PheCode at a time: rows the PheCode marks as excluded controls (NA)
#     are dropped for that test only, so each code keeps its own control set.
#
# Invoked by 04_run_dwas.R, once per cohort x model; see that script for the
# full model matrix.
#
# Usage:
#   Rscript dwas_regression.R \
#     --phecode    $READY_DIR/phecode_matrix_group1.rds \
#     --phenotype  $READY_DIR/phenotype_group1.csv \
#     --outcome    taste_2w_strict \
#     --covariates "age_baseline,sex,assess_centre_id,townsend,smoking,drink,BMI,total_dx_count" \
#     --factor-vars "sex,assess_centre_id,smoking,drink" \
#     --min-cases  100 \
#     --output     $RESULTS_DIR/group1_primary.csv \
#     --n-jobs     10
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

# =============================================================================
# Command-line interface
# =============================================================================
option_list <- list(
  make_option("--phecode",     type = "character", help = "PheCode matrix RDS file"),
  make_option("--phenotype",   type = "character", help = "Phenotype CSV"),
  make_option("--outcome",     type = "character", default = "taste_2w_strict", help = "Outcome column name"),
  make_option("--covariates",  type = "character", help = "Comma-separated covariate list"),
  make_option("--factor-vars", type = "character", default = "", help = "Comma-separated factor variables"),
  make_option("--min-cases",   type = "integer",   default = 100, help = "Minimum PheCode case count to test [default %default]"),
  make_option("--output",      type = "character", help = "Output CSV path"),
  make_option("--n-jobs",      type = "integer",   default = 1, help = "Parallel cores [default %default]"),
  make_option("--phecodes",    type = "character", default = "", help = "Optional: restrict to these PheCodes (comma-separated)")
)

opt <- parse_args(OptionParser(option_list = option_list))

# =============================================================================
# Step 1: load data
# =============================================================================
cat("=== Loading data ===\n")

# PheCode matrix (eid + one logical column per PheCode: TRUE/FALSE/NA)
phecode_mat <- as.data.table(readRDS(opt$phecode))
setnames(phecode_mat, "id", "eid", skip_absent = TRUE)
cat(sprintf("  PheCode matrix: %d people x %d PheCodes\n",
            nrow(phecode_mat), ncol(phecode_mat) - 1))

pheno <- fread(opt$phenotype)
cat(sprintf("  Phenotype file: %d people x %d columns\n", nrow(pheno), ncol(pheno)))

# =============================================================================
# Step 2: run-time merge
# =============================================================================
cat("\n=== Merging ===\n")

# Inner join by eid: keep only people present in both
merged <- merge(phecode_mat, pheno, by = "eid")
cat(sprintf("  After merge: %d people\n", nrow(merged)))

# Keep only people with a non-missing outcome
merged <- merged[!is.na(get(opt$outcome))]
cat(sprintf("  With outcome: %d people (case=%d, control=%d)\n",
            nrow(merged),
            sum(merged[[opt$outcome]] == 1, na.rm = TRUE),
            sum(merged[[opt$outcome]] == 0, na.rm = TRUE)))

# =============================================================================
# Step 3: parse covariates and factor variables
# =============================================================================
covariates  <- strsplit(opt$covariates, ",")[[1]]
factor_vars <- if (opt$`factor-vars` != "") strsplit(opt$`factor-vars`, ",")[[1]] else character(0)

missing_vars <- setdiff(c(covariates, opt$outcome), names(merged))
if (length(missing_vars) > 0) {
  stop(sprintf("variables not found: %s", paste(missing_vars, collapse = ", ")))
}

for (fv in factor_vars) {
  if (fv %in% names(merged)) {
    merged[[fv]] <- as.factor(merged[[fv]])
  }
}

cat(sprintf("  Covariates (%d): %s\n", length(covariates), paste(covariates, collapse = ", ")))
cat(sprintf("  Factor vars (%d): %s\n", length(factor_vars), paste(factor_vars, collapse = ", ")))

# =============================================================================
# Step 4: determine which PheCodes to analyse
# =============================================================================
all_phecodes <- setdiff(names(phecode_mat), "eid")

if (opt$phecodes != "") {
  target_phecodes <- strsplit(opt$phecodes, ",")[[1]]
  all_phecodes <- intersect(all_phecodes, target_phecodes)
  cat(sprintf("  Restricted to %d specified PheCodes\n", length(all_phecodes)))
}

# Drop low-frequency PheCodes (case count below threshold)
case_counts <- sapply(all_phecodes, function(pc) sum(merged[[pc]] == TRUE, na.rm = TRUE))
analyzable  <- all_phecodes[case_counts >= opt$`min-cases`]

cat(sprintf("  Analysable PheCodes (case >= %d): %d / %d\n",
            opt$`min-cases`, length(analyzable), length(all_phecodes)))

# =============================================================================
# Step 5: regression loop (glm; Firth fallback on warning)
# =============================================================================
cat(sprintf("\n=== Regression (%d PheCodes) ===\n", length(analyzable)))

formula_base <- sprintf("%s ~ phecode_exposure + %s", opt$outcome, paste(covariates, collapse = " + "))
cat(sprintf("  Formula: %s\n\n", formula_base))

run_one_phecode <- function(pc) {
  # cases + controls for this PheCode (NA = PheCode-specific exclusion, dropped)
  sub_dt <- merged[!is.na(get(pc)), ]
  sub_dt[, phecode_exposure := as.integer(get(pc) == TRUE)]

  n_total   <- nrow(sub_dt)
  n_case_pc <- sum(sub_dt$phecode_exposure == 1, na.rm = TRUE)
  n_ctrl_pc <- sum(sub_dt$phecode_exposure == 0, na.rm = TRUE)

  # complete cases on all model terms
  all_vars <- c(opt$outcome, "phecode_exposure", covariates)
  sub_dt <- sub_dt[complete.cases(sub_dt[, ..all_vars])]

  result <- tryCatch({
    fit <- glm(as.formula(formula_base), data = sub_dt, family = binomial())
    coef_row <- summary(fit)$coefficients["phecode_exposure", ]

    data.table(
      phecode    = pc,
      beta       = coef_row["Estimate"],
      se         = coef_row["Std. Error"],
      z          = coef_row["z value"],
      pval       = coef_row["Pr(>|z|)"],
      or         = exp(coef_row["Estimate"]),
      or_lower   = exp(coef_row["Estimate"] - 1.96 * coef_row["Std. Error"]),
      or_upper   = exp(coef_row["Estimate"] + 1.96 * coef_row["Std. Error"]),
      n_total    = nrow(sub_dt),
      n_case_taste = sum(sub_dt[[opt$outcome]] == 1),
      n_case_phecode = n_case_pc,
      n_ctrl_phecode = n_ctrl_pc,
      converged  = TRUE,
      method     = "glm"
    )
  }, warning = function(w) {
    # glm convergence / separation warning -> Firth penalised likelihood
    tryCatch({
      if (!requireNamespace("logistf", quietly = TRUE)) stop("logistf not installed")
      fit_f <- logistf::logistf(as.formula(formula_base), data = sub_dt)
      idx <- which(names(fit_f$coefficients) == "phecode_exposure")

      data.table(
        phecode    = pc,
        beta       = fit_f$coefficients[idx],
        se         = sqrt(diag(vcov(fit_f)))[idx],
        z          = NA_real_,
        pval       = fit_f$prob[idx],
        or         = exp(fit_f$coefficients[idx]),
        or_lower   = exp(fit_f$ci.lower[idx]),
        or_upper   = exp(fit_f$ci.upper[idx]),
        n_total    = nrow(sub_dt),
        n_case_taste = sum(sub_dt[[opt$outcome]] == 1),
        n_case_phecode = n_case_pc,
        n_ctrl_phecode = n_ctrl_pc,
        converged  = TRUE,
        method     = "firth"
      )
    }, error = function(e2) {
      data.table(
        phecode = pc, beta = NA, se = NA, z = NA, pval = NA,
        or = NA, or_lower = NA, or_upper = NA,
        n_total = nrow(sub_dt), n_case_taste = NA,
        n_case_phecode = n_case_pc, n_ctrl_phecode = n_ctrl_pc,
        converged = FALSE, method = "failed"
      )
    })
  }, error = function(e) {
    data.table(
      phecode = pc, beta = NA, se = NA, z = NA, pval = NA,
      or = NA, or_lower = NA, or_upper = NA,
      n_total = 0, n_case_taste = NA,
      n_case_phecode = n_case_pc, n_ctrl_phecode = n_ctrl_pc,
      converged = FALSE, method = "error"
    )
  })

  return(result)
}

if (opt$`n-jobs` > 1) {
  results <- parallel::mclapply(analyzable, run_one_phecode, mc.cores = opt$`n-jobs`)
} else {
  results <- lapply(analyzable, function(pc) {
    if (which(analyzable == pc) %% 50 == 0) {
      cat(sprintf("  progress: %d / %d\n", which(analyzable == pc), length(analyzable)))
    }
    run_one_phecode(pc)
  })
}

results_dt <- rbindlist(results, fill = TRUE)

# =============================================================================
# Step 6: multiple-testing correction + save
# =============================================================================
cat("\n=== Multiple-testing correction ===\n")

results_dt[, pval_fdr  := p.adjust(pval, method = "fdr")]          # Benjamini-Hochberg
results_dt[, pval_bonf := p.adjust(pval, method = "bonferroni")]   # Bonferroni

# Attach PheCode descriptions (PheWAS catalogue)
if (requireNamespace("PheWAS", quietly = TRUE)) {
  desc <- as.data.table(PheWAS::pheinfo)
  desc_cols <- intersect(c("phecode", "description", "group", "groupnum"), names(desc))
  if (length(desc_cols) > 1) {
    results_dt <- merge(results_dt, desc[, ..desc_cols], by = "phecode", all.x = TRUE)
  }
}

# Inflation factor: a run diagnostic only. This scan is outcome-anchored over
# health-relevant PheCodes, so a large share of the tests are genuinely
# non-null and this number is expected to exceed 1 without indicating
# confounding. It is not reported in the manuscript.
valid_p <- results_dt$pval[!is.na(results_dt$pval)]
lambda_gc <- median(qchisq(1 - valid_p, df = 1)) / qchisq(0.5, df = 1)

cat(sprintf("  Total PheCodes: %d\n", nrow(results_dt)))
cat(sprintf("  Converged: %d (%.1f%%)\n", sum(results_dt$converged), mean(results_dt$converged) * 100))
cat(sprintf("  Lambda GC: %.3f\n", lambda_gc))
cat(sprintf("  Nominal (p<0.05): %d\n", sum(results_dt$pval < 0.05, na.rm = TRUE)))
cat(sprintf("  FDR < 0.05: %d\n", sum(results_dt$pval_fdr < 0.05, na.rm = TRUE)))
cat(sprintf("  Bonferroni < 0.05: %d\n", sum(results_dt$pval_bonf < 0.05, na.rm = TRUE)))

setorder(results_dt, pval)
dir.create(dirname(opt$output), showWarnings = FALSE, recursive = TRUE)
fwrite(results_dt, opt$output)
# Export and dx upload to RAP  (this model's result table)
cat(sprintf("\n  Saved: %s\n", opt$output))
