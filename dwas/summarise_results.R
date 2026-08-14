#!/usr/bin/env Rscript
# =============================================================================
# Disease-wide association scan — per-model comparison table
# =============================================================================
#
# Aggregates the per-model result CSVs into one comparison table: genomic
# inflation, nominal / suggestive / FDR / Bonferroni hit counts, and the top
# PheCode per model. Invoked by 04_run_dwas.R once every model has finished;
# it can also be run on its own against an existing results directory.
#
# The inflation factor is a run diagnostic only. An outcome-anchored scan over
# health-relevant PheCodes has a genuinely non-null test distribution, so this
# number is not evidence of confounding and is not reported in the manuscript.
#
# Usage:
#   Rscript summarise_results.R [results_dir]
#   (default: $PROJECT_DIR/output/dwas/results)
# =============================================================================
suppressPackageStartupMessages(library(data.table))

args        <- commandArgs(trailingOnly = TRUE)
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
RESULTS_DIR <- if (length(args) >= 1) args[1] else
  file.path(PROJECT_DIR, "output", "dwas", "results")
if (!dir.exists(RESULTS_DIR)) stop("results directory not found: ", RESULTS_DIR)
cat("Reading results from:", RESULTS_DIR, "\n")

files <- list.files(RESULTS_DIR, pattern = "^group[123]_.*\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("no per-model result files in ", RESULTS_DIR)

results <- rbindlist(lapply(files, function(f) {
  dt <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)

  model_name <- gsub("\\.csv$", "", basename(f))

  valid_p <- dt$pval[!is.na(dt$pval)]
  lambda_gc <- if (length(valid_p) > 0) {
    median(qchisq(1 - valid_p, df = 1)) / qchisq(0.5, df = 1)
  } else NA

  top_row <- dt[which.min(pval)]

  data.table(
    Model           = model_name,
    N_PheCodes      = nrow(dt),
    N_Converged     = sum(dt$converged, na.rm = TRUE),
    Converge_Rate   = round(mean(dt$converged, na.rm = TRUE) * 100, 1),
    N_Nominal       = sum(dt$pval < 0.05, na.rm = TRUE),
    N_Suggestive    = sum(dt$pval < 0.001, na.rm = TRUE),
    N_FDR           = sum(dt$pval_fdr < 0.05, na.rm = TRUE),
    N_Bonferroni    = sum(dt$pval_bonf < 0.05, na.rm = TRUE),
    Lambda          = round(lambda_gc, 3),
    Top_PheCode     = top_row$phecode,
    Top_Description = if ("description" %in% names(top_row)) top_row$description else NA,
    Top_OR          = round(top_row$or, 3),
    Top_P           = format(top_row$pval, scientific = TRUE, digits = 2)
  )
}), fill = TRUE)

setorder(results, Model)

cat("\n", strrep("=", 90), "\n", sep = "")
cat("  DISEASE-WIDE SCAN - BATCH SUMMARY\n")
cat(strrep("=", 90), "\n\n")

for (grp in c("group1", "group2", "group3")) {
  grp_results <- results[grepl(grp, Model)]
  if (nrow(grp_results) > 0) {
    cat(sprintf("\n--- %s ---\n", toupper(grp)))
    print_cols <- c("Model", "Lambda", "N_Nominal", "N_FDR", "N_Bonferroni",
                    "Top_PheCode", "Top_OR", "Top_P")
    print_cols <- intersect(print_cols, names(grp_results))
    print(grp_results[, ..print_cols])
  }
}

# Group 1 model comparison (primary vs no-BMI)
cat("\n\n--- Group 1 model comparison ---\n")
g1 <- results[grepl("group1", Model)]
if (nrow(g1) > 0) {
  primary <- g1[grepl("primary", Model)]
  no_bmi  <- g1[grepl("no_bmi", Model)]
  if (nrow(primary) > 0 & nrow(no_bmi) > 0) {
    cat(sprintf("  Primary:  Lambda=%.3f  FDR=%d\n", primary$Lambda[1], primary$N_FDR[1]))
    cat(sprintf("  No BMI:   Lambda=%.3f  FDR=%d\n", no_bmi$Lambda, no_bmi$N_FDR))
    cat(sprintf("  FDR change after dropping BMI: %d\n", abs(primary$N_FDR[1] - no_bmi$N_FDR)))
  }
}

# Group 1 vs Group 2 held-out validation
cat("\n--- Group 1 vs Group 2 validation ---\n")
g1_primary <- results[grepl("group1_primary$", Model)]
g2_primary <- results[grepl("group2_primary", Model)]
if (nrow(g1_primary) > 0 & nrow(g2_primary) > 0) {
  cat(sprintf("  Group 1: %d FDR-sig PheCodes, Lambda=%.3f\n",
              g1_primary$N_FDR, g1_primary$Lambda))
  cat(sprintf("  Group 2: %d FDR-sig PheCodes, Lambda=%.3f\n",
              g2_primary$N_FDR, g2_primary$Lambda))
}

out_file <- file.path(RESULTS_DIR, "dwas_models_comparison.csv")
fwrite(results, out_file)
cat("\n\nSaved:", out_file, "\n")

# =============================================================================
# FDR-significant table (reporting)
# =============================================================================
#
# The labelled table of FDR-significant hits, built from the smell-adjusted
# discovery model — the conservative one, and the model the reported disease
# hits are taken from. PheCode descriptions and clinical categories come from
# the PheWAS catalogue; the ICD-10 column is the alphabetically first
# three-character code the PheCode maps back to, as a reading aid only (most
# PheCodes cover several).
#
# Rollup parent codes are kept, as everywhere else in this scan: dropping them
# here would silently change which hits are reported.

SOURCE_MODEL <- "group1_primary+smell.csv"
src <- file.path(RESULTS_DIR, SOURCE_MODEL)

if (!file.exists(src)) {
  cat("\n[skip] FDR-significant table:", SOURCE_MODEL, "not found\n")
} else {
  r <- fread(src)
  r[, phecode := as.character(phecode)]

  if (!"description" %in% names(r) && requireNamespace("PheWAS", quietly = TRUE)) {
    desc <- as.data.table(PheWAS::pheinfo)
    desc[, phecode := as.character(phecode)]
    r <- merge(r, desc[, .(phecode, description, group)], by = "phecode", all.x = TRUE)
  }
  # PheCodes the catalogue does not describe keep an empty Description: the
  # figure code falls back to the raw identifier, and inventing a label here
  # would put a different string on the published figure.
  if (!"description" %in% names(r)) r[, description := NA_character_]
  if (!"group" %in% names(r))       r[, group := NA_character_]

  # PheCode -> first three-character ICD-10 code
  r[, icd10_primary := NA_character_]
  if (requireNamespace("PheWAS", quietly = TRUE)) {
    pm <- as.data.table(PheWAS::phecode_map)
    pm <- pm[vocabulary_id == "ICD10CM"]
    pm[, `:=`(phecode = as.character(phecode), icd_3char = substr(code, 1, 3))]
    back <- pm[, .(icd10 = sort(unique(icd_3char))[1]), by = phecode]
    r[back, icd10_primary := i.icd10, on = "phecode"]
  }

  hits <- r[pval_fdr < 0.05][order(pval)]
  pub <- hits[, .(
    ICD10       = icd10_primary,
    PheCode     = phecode,
    Description = description,
    Category    = group,
    N_Case      = n_case_phecode,
    N_Control   = n_ctrl_phecode,
    OR          = sprintf("%.2f", or),
    CI_95       = sprintf("%.2f-%.2f", or_lower, or_upper),
    P_value     = formatC(pval, format = "e", digits = 2),
    FDR         = formatC(pval_fdr, format = "e", digits = 2),
    Direction   = fifelse(beta > 0, "Risk", "Protective")
  )]

  fdr_file <- file.path(RESULTS_DIR, "dwas_fdr_significant_table.csv")
  fwrite(pub, fdr_file)
  cat(sprintf("FDR-significant (%s): %d / %d PheCodes -> %s\n",
              sub("\\.csv$", "", SOURCE_MODEL), nrow(pub), nrow(r), basename(fdr_file)))
  cat(sprintf("  risk %d | protective %d | without an ICD-10 back-map %d\n",
              sum(pub$Direction == "Risk"), sum(pub$Direction == "Protective"),
              sum(is.na(pub$ICD10))))
}
# Export and dx upload to RAP  (the model comparison table and the
# FDR-significant table)
