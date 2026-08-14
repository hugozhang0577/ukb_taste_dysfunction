#!/usr/bin/env Rscript
# =============================================================================
# 08_build_tables.R  (table assembly)
#
# Consolidate the modelling outputs into publication-ready tables for the
# manuscript. Reads the existing eval CSVs + ml_ready RDS, reformats with
# human-readable names, writes to model_reports/fig_tables/.
#
# Main: cohort characteristics, model performance (G1 CV + G2/G3), top-20
#       features (M1_TierA), fairness subgroups (M1_TierA).
# Supplementary: hyperparameters, DeLong (discovery + held-out), fairness
#       (all models), top-20 features (all models),
#       no-smell sensitivity (if 07 has run). Plus an INDEX.csv.
# =============================================================================

rm(list = ls()); gc()
while (sink.number() > 0) sink()
suppressMessages({ library(data.table) })

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)

ML_DIR    <- "output/ml_ready"
MODEL_DIR <- "output/models/xgboost"
P6_REP    <- "output/model_reports"
EVAL_DIR  <- file.path(P6_REP, "eval")
TAB_DIR   <- file.path(P6_REP, "fig_tables")
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
sink(file.path(P6_REP, "build_tables_log.txt"), split = TRUE)
on.exit({ while (sink.number() > 0) sink() }, add = TRUE)

cat("=== build publication tables ===\n")
OUTCOME <- "taste_2w_strict"
fmt_mean_sd <- function(x) { x <- x[!is.na(x)]; if (!length(x)) NA_character_ else sprintf("%.1f (%.1f)", mean(x), sd(x)) }
fmt_pct <- function(x, lvl = 1) { x <- x[!is.na(x)]; if (!length(x)) NA_character_ else sprintf("%d (%.1f%%)", sum(x == lvl), 100 * mean(x == lvl)) }

# ---- Table 1: cohort characteristics ----------------------------------------
cat("[Table 1] cohort characteristics\n")
tab1_rows <- list()
for (g in 1:3) {
  rds <- file.path(ML_DIR, sprintf("group%d_full.rds", g)); if (!file.exists(rds)) next
  d <- as.data.table(readRDS(rds))[!is.na(get(OUTCOME))]
  for (grp in c("all", "case", "control")) {
    sub <- switch(grp, all = d, case = d[get(OUTCOME) == 1], control = d[get(OUTCOME) == 0])
    tw <- if ("townsend" %in% names(sub)) fmt_mean_sd(sub$townsend) else
      if ("townsend_deprivation_index" %in% names(sub)) fmt_mean_sd(sub$townsend_deprivation_index) else NA_character_
    tab1_rows[[paste(g, grp, sep = "_")]] <- data.table(
      cohort = paste0("G", g), stratum = grp, n = nrow(sub),
      case_rate_pct = if (grp == "all") sprintf("%.2f", 100 * mean(sub[[OUTCOME]])) else NA_character_,
      age_mean_sd = if ("age" %in% names(sub)) fmt_mean_sd(sub$age) else NA_character_,
      sex_female_n_pct = if ("sex" %in% names(sub)) fmt_pct(sub$sex, 0) else NA_character_,
      townsend_mean_sd = tw,
      bmi_mean_sd = if ("bmi" %in% names(sub)) fmt_mean_sd(sub$bmi) else NA_character_,
      smell_any_n_pct = if ("smell_any" %in% names(sub)) fmt_pct(sub$smell_any, 1) else NA_character_)
  }
}
fwrite(rbindlist(tab1_rows, fill = TRUE), file.path(TAB_DIR, "table1_cohort_characteristics.csv"))

# ---- Table 2: model performance ---------------------------------------------
cat("[Table 2] model performance\n")
best    <- fread(file.path(P6_REP, "xgb_best_params.csv"))
eval_g1 <- fread(file.path(EVAL_DIR, "xgb_eval_summary.csv"))
ext_p   <- fread(file.path(EVAL_DIR, "heldout_val_summary.csv"))[variant == "platt"]

g1_rows <- eval_g1[, .(model_id, cohort = "G1_internal_OOF", n, n_cases,
  case_rate_pct = round(100 * case_rate, 2), AUC = round(auc, 4), AUC_lo = round(auc_lo, 4),
  AUC_hi = round(auc_hi, 4), PR_AUC = round(pr_auc, 4), Brier = round(brier, 4),
  LogLoss = round(logloss, 4))]
ext_rows <- ext_p[, .(model_id, cohort = toupper(cohort), n, n_cases,
  case_rate_pct = round(100 * n_cases / n, 2), AUC = round(auc, 4), AUC_lo = round(auc_lo, 4),
  AUC_hi = round(auc_hi, 4), PR_AUC = round(pr_auc, 4), Brier = round(brier, 4),
  LogLoss = round(logloss, 4))]
table2 <- rbind(g1_rows, ext_rows, fill = TRUE); setorder(table2, model_id, cohort)
fwrite(table2, file.path(TAB_DIR, "table2_model_performance.csv"))

# ---- Table 3: top-20 features (M1_TierA) ------------------------------------
cat("[Table 3] top-20 features (M1_TierA)\n")
imp1 <- fread(file.path(MODEL_DIR, "M1_TierA_importance.csv"))[order(-Gain)][1:20]
imp1[, Rank := seq_len(.N)]; setcolorder(imp1, c("Rank", "Feature", "Gain", "Cover", "Frequency"))
fwrite(imp1, file.path(TAB_DIR, "table3_top20_features_M1_TierA.csv"))

# ---- Table 4: fairness subgroup (M1_TierA) ----------------------------------
cat("[Table 4] fairness subgroup (M1_TierA)\n")
fair <- fread(file.path(EVAL_DIR, "fairness_subgroup_metrics.csv"))
fair1 <- fair[model_id == "M1_TierA", .(subgroup_var, subgroup_level, n, n_cases,
  case_rate_pct = round(100 * case_rate, 2), AUC = round(auc, 3), AUC_lo = round(auc_lo, 3),
  AUC_hi = round(auc_hi, 3), Brier = round(brier, 4),
  Sens_at_Spec90 = round(sens_at_spec90, 3), PPV_top10 = round(ppv_top10, 3))]
fwrite(fair1, file.path(TAB_DIR, "table4_fairness_subgroup_M1_TierA.csv"))

# ---- supplementary tables ---------------------------------------------------
fwrite(best[, .(model_id, sample_subset, n_train, n_cases, n_features, max_depth, eta, subsample,
  colsample_bytree, min_child_weight, nrounds, cv_mean_AUC = round(cv_mean_auc, 4),
  cv_sd_AUC = round(cv_sd_auc, 4), cv_mean_Brier = round(cv_mean_brier, 4),
  cv_mean_LogLoss = round(cv_mean_logloss, 4))],
  file.path(TAB_DIR, "suppl_table1_xgb_best_hyperparameters.csv"))

st2 <- fread(file.path(EVAL_DIR, "xgb_delong_pairwise.csv"))
fwrite(st2[, .(model_a, model_b, AUC_a = round(auc_a, 4), AUC_b = round(auc_b, 4),
  delta_AUC = round(delta_auc, 4), z = round(z, 3), p_value = signif(p_value, 3), n_common)],
  file.path(TAB_DIR, "suppl_table2_delong_internal_g1.csv"))

st3 <- fread(file.path(EVAL_DIR, "heldout_val_delong_vs_g1.csv"))
fwrite(st3[, .(model_id, cohort = toupper(cohort), AUC_G1 = round(auc_g1, 4), AUC_ext = round(auc_ext, 4),
  delta_AUC = round(delta_auc, 4), z = round(z, 3), p_value = signif(p_value, 3))],
  file.path(TAB_DIR, "suppl_table3_delong_g1_vs_g2g3.csv"))

st5 <- fair[, .(model_id, subgroup_var, subgroup_level, n, n_cases, case_rate_pct = round(100 * case_rate, 2),
  AUC = round(auc, 3), AUC_lo = round(auc_lo, 3), AUC_hi = round(auc_hi, 3), Brier = round(brier, 4),
  Sens_at_Spec90 = round(sens_at_spec90, 3), PPV_top10 = round(ppv_top10, 3))]
setorder(st5, model_id, subgroup_var, subgroup_level)
fwrite(st5, file.path(TAB_DIR, "suppl_table5_fairness_all_models.csv"))
fwrite(fread(file.path(EVAL_DIR, "fairness_gaps_summary.csv")), file.path(TAB_DIR, "suppl_table5b_fairness_gaps.csv"))

imp_all <- list()
for (mid in best$model_id) {
  f <- file.path(MODEL_DIR, sprintf("%s_importance.csv", mid)); if (!file.exists(f)) next
  x <- fread(f)[order(-Gain)][1:20]; x[, Rank := seq_len(.N)]; x[, model_id := mid]; imp_all[[mid]] <- x
}
st7 <- rbindlist(imp_all, fill = TRUE); setcolorder(st7, c("model_id", "Rank", "Feature", "Gain", "Cover", "Frequency"))
fwrite(st7, file.path(TAB_DIR, "suppl_table7_top20_features_all_models.csv"))

# ---- supplementary: no-smell sensitivity (if 07 has run) --------------------
sens_sum_f <- file.path(EVAL_DIR, "sensitivity_no_smell_summary.csv")
sens_dl_f  <- file.path(EVAL_DIR, "sensitivity_no_smell_delong_vs_original.csv")
sens_feat_f <- file.path(EVAL_DIR, "sensitivity_no_smell_feature_counts.csv")
if (file.exists(sens_sum_f) && file.exists(sens_dl_f)) {
  sens_sum <- fread(sens_sum_f); sens_dl <- fread(sens_dl_f)
  fwrite(sens_sum[, .(variant_id, original_id, cohort, n, n_cases, AUC = round(auc, 4),
    AUC_lo = round(auc_lo, 4), AUC_hi = round(auc_hi, 4), PR_AUC = round(pr_auc, 4),
    Brier = round(brier, 4), LogLoss = round(logloss, 4))][order(variant_id, cohort)],
    file.path(TAB_DIR, "suppl_table8a_sensitivity_no_smell_performance.csv"))
  fwrite(sens_dl[, .(variant_id, original_id, cohort, AUC_no_smell = round(AUC_no_smell, 4),
    AUC_original = round(AUC_original, 4), delta_AUC = round(delta_AUC, 4), z = round(z, 2),
    p_value = signif(p_value, 3), n_common, paired)][order(variant_id, cohort)],
    file.path(TAB_DIR, "suppl_table8b_sensitivity_no_smell_delong.csv"))
  if (file.exists(sens_feat_f)) fwrite(fread(sens_feat_f), file.path(TAB_DIR, "suppl_table8c_sensitivity_no_smell_features_removed.csv"))
} else {
  cat("  no-smell sensitivity CSVs not found - skipping Suppl 8\n")
}

# ---- INDEX ------------------------------------------------------------------
files <- setdiff(list.files(TAB_DIR, pattern = "\\.csv$"), "INDEX.csv")
fwrite(data.table(file = files), file.path(TAB_DIR, "INDEX.csv"))
cat(sprintf("\n%d table files written to %s\n", length(files), TAB_DIR))
cat("=== 08 complete ===\n")
