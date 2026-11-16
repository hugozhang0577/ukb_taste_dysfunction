# Reporting tables for the classification models: cohort characteristics, discrimination by
# model and cohort, hyperparameters, DeLong comparisons, top features, subgroup metrics and
# the no-smell sensitivity analysis.
library(data.table)
setwd(Sys.getenv("PROJECT_DIR", "."))

MODEL_DIR <- "output/models/xgboost"
REP_DIR   <- "output/model_reports"
EVAL_DIR  <- file.path(REP_DIR, "eval")
TAB_DIR   <- file.path(REP_DIR, "fig_tables")
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
OUTCOME <- "taste_2w_strict"

fmt_mean_sd <- function(x) { x <- x[!is.na(x)]; if (!length(x)) NA_character_ else sprintf("%.1f (%.1f)", mean(x), sd(x)) }
fmt_pct <- function(x, lvl = 1) { x <- x[!is.na(x)]; if (!length(x)) NA_character_ else sprintf("%d (%.1f%%)", sum(x == lvl), 100 * mean(x == lvl)) }

# Cohort characteristics, all / cases / controls
tab1 <- list()
for (g in 1:3) {
  d <- readRDS(file.path("output/ml_ready", sprintf("group%d_full.rds", g)))[!is.na(get(OUTCOME))]
  for (grp in c("all", "case", "control")) {
    sub <- switch(grp, all = d, case = d[get(OUTCOME) == 1], control = d[get(OUTCOME) == 0])
    tab1[[paste(g, grp)]] <- data.table(
      cohort = paste0("G", g), stratum = grp, n = nrow(sub),
      case_rate_pct = if (grp == "all") sprintf("%.2f", 100 * mean(sub[[OUTCOME]])) else NA_character_,
      age_mean_sd = fmt_mean_sd(sub$age), sex_female_n_pct = fmt_pct(sub$sex, 0),
      townsend_mean_sd = fmt_mean_sd(sub$townsend),
      bmi_mean_sd = fmt_mean_sd(sub$BMI),
      smell_any_n_pct = fmt_pct(sub$smell_any, 1))
  }
}
fwrite(rbindlist(tab1, fill = TRUE), file.path(TAB_DIR, "cohort_characteristics.csv"))

# Discrimination by model: group1 out-of-fold and held-out group2/group3 (Platt-calibrated)
best    <- fread(file.path(REP_DIR, "xgb_best_params.csv"))
eval_g1 <- fread(file.path(EVAL_DIR, "xgb_eval_summary.csv"))
heldout <- fread(file.path(EVAL_DIR, "heldout_val_summary.csv"))[variant == "platt"]
perf_cols <- function(d) d[, .(AUC = round(auc, 4), AUC_lo = round(auc_lo, 4), AUC_hi = round(auc_hi, 4),
                               PR_AUC = round(pr_auc, 4), Brier = round(brier, 4), LogLoss = round(logloss, 4))]
discrimination <- rbind(
  cbind(eval_g1[, .(model_id, cohort = "G1_internal_OOF", n, n_cases, case_rate_pct = round(100 * case_rate, 2))], perf_cols(eval_g1)),
  cbind(heldout[, .(model_id, cohort = toupper(cohort), n, n_cases, case_rate_pct = round(100 * n_cases / n, 2))], perf_cols(heldout)),
  fill = TRUE)
setorder(discrimination, model_id, cohort)
fwrite(discrimination, file.path(TAB_DIR, "discrimination_by_model.csv"))

# Top 20 features by gain
top20 <- function(mid) {
  x <- fread(file.path(MODEL_DIR, sprintf("%s_importance.csv", mid)))[order(-Gain)][1:20]
  x[, Rank := seq_len(.N)][]
}
fwrite(setcolorder(top20("model1"), c("Rank", "Feature", "Gain", "Cover", "Frequency")),
       file.path(TAB_DIR, "top_features_model1.csv"))
top_all <- rbindlist(lapply(best$model_id, function(mid) top20(mid)[, model_id := mid]), fill = TRUE)
fwrite(setcolorder(top_all, c("model_id", "Rank", "Feature", "Gain", "Cover", "Frequency")),
       file.path(TAB_DIR, "top_features_all_models.csv"))

# Subgroup (fairness) metrics
fair <- fread(file.path(EVAL_DIR, "fairness_subgroup_metrics.csv"))
fair_cols <- function(d) d[, .(subgroup_var, subgroup_level, n, n_cases, case_rate_pct = round(100 * case_rate, 2),
  AUC = round(auc, 3), AUC_lo = round(auc_lo, 3), AUC_hi = round(auc_hi, 3), Brier = round(brier, 4),
  Sens_at_Spec90 = round(sens_at_spec90, 3), PPV_top10 = round(ppv_top10, 3))]
fwrite(fair_cols(fair[model_id == "model1"]), file.path(TAB_DIR, "subgroup_model1.csv"))
fair_all <- cbind(fair[, .(model_id)], fair_cols(fair))
setorder(fair_all, model_id, subgroup_var, subgroup_level)
fwrite(fair_all, file.path(TAB_DIR, "subgroup_all_models.csv"))
fwrite(fread(file.path(EVAL_DIR, "fairness_gaps_summary.csv")), file.path(TAB_DIR, "subgroup_gaps.csv"))

# Hyperparameters and cross-validated performance of the selected configuration
fwrite(best[, .(model_id, sample_subset, n_train, n_cases, n_features, max_depth, eta, subsample,
  colsample_bytree, min_child_weight, nrounds, cv_mean_AUC = round(cv_mean_auc, 4),
  cv_sd_AUC = round(cv_sd_auc, 4), cv_mean_Brier = round(cv_mean_brier, 4),
  cv_mean_LogLoss = round(cv_mean_logloss, 4))],
  file.path(TAB_DIR, "hyperparameters_by_model.csv"))

# DeLong tests: between models in group1, and group1 vs held-out cohorts
dl <- fread(file.path(EVAL_DIR, "xgb_delong_pairwise.csv"))
fwrite(dl[, .(model_a, model_b, AUC_a = round(auc_a, 4), AUC_b = round(auc_b, 4),
  delta_AUC = round(delta_auc, 4), z = round(z, 3), p_value = signif(p_value, 3), n_common)],
  file.path(TAB_DIR, "delong_pairwise_group1.csv"))
dh <- fread(file.path(EVAL_DIR, "heldout_val_delong_vs_g1.csv"))
fwrite(dh[, .(model_id, cohort = toupper(cohort), AUC_G1 = round(auc_g1, 4), AUC_heldout = round(auc_heldout, 4),
  delta_AUC = round(delta_auc, 4), z = round(z, 3), p_value = signif(p_value, 3))],
  file.path(TAB_DIR, "delong_group1_vs_heldout.csv"))

# No-smell sensitivity analysis
sens_sum <- fread(file.path(EVAL_DIR, "sensitivity_no_smell_summary.csv"))
fwrite(sens_sum[, .(variant_id, original_id, cohort, n, n_cases, AUC = round(auc, 4),
  AUC_lo = round(auc_lo, 4), AUC_hi = round(auc_hi, 4), PR_AUC = round(pr_auc, 4),
  Brier = round(brier, 4), LogLoss = round(logloss, 4))][order(variant_id, cohort)],
  file.path(TAB_DIR, "no_smell_discrimination.csv"))
sens_dl <- fread(file.path(EVAL_DIR, "sensitivity_no_smell_delong_vs_original.csv"))
fwrite(sens_dl[, .(variant_id, original_id, cohort, AUC_no_smell = round(AUC_no_smell, 4),
  AUC_original = round(AUC_original, 4), delta_AUC = round(delta_AUC, 4), z = round(z, 2),
  p_value = signif(p_value, 3), n_common, paired)][order(variant_id, cohort)],
  file.path(TAB_DIR, "no_smell_delong.csv"))
fwrite(fread(file.path(EVAL_DIR, "sensitivity_no_smell_feature_counts.csv")),
       file.path(TAB_DIR, "no_smell_features_removed.csv"))
