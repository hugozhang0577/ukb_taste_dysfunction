#!/usr/bin/env Rscript
# =============================================================================
# 07_heldout_validation.R  (held-out cross-population validation)
#
# Apply the G1-trained boosters (06b) + G1-fitted calibrators (06h) to the two
# held-out ancestry-stratified cohorts (G2 Other White; G3 Non-White). These are
# UKB internal held-out partitions, not external samples. STRICTLY held out: no
# training, no tuning, no calibrator refitting (a transportability test).
#
# Per model x cohort: AUC (95% DeLong CI), PR-AUC, Brier, LogLoss, calibration
# intercept/slope of the calibrated predictions, sensitivity @ spec {0.80,0.90,
# 0.95}, top-k% capture; plus a DeLong comparison of G1 OOF vs each cohort.
#
# Input:  models/xgboost/{model_id}_booster.rds + calibrators/ + _oof_calibrated.rds
#         data/ml_ready/group{2,3}_{full,olink,nmr,olink_nmr}.rds
# Output: models/xgboost/predictions/{model_id}_{g2,g3}_preds.rds
#         model_reports/eval/heldout_val_{summary,thresholded,delong_vs_g1}.csv + plots
# =============================================================================

rm(list = ls()); gc()
while (sink.number() > 0) sink()
suppressMessages({ library(data.table); library(xgboost); library(pROC); library(PRROC) })

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)

ML_DIR    <- "output/ml_ready"
MODEL_DIR <- "output/models/xgboost"
CAL_DIR   <- file.path(MODEL_DIR, "calibrators")
PRED_DIR  <- file.path(MODEL_DIR, "predictions")
P6_REP    <- "output/model_reports"
EVAL_DIR  <- file.path(P6_REP, "eval")
PLOT_DIR  <- file.path(EVAL_DIR, "plots")
for (d in c(PRED_DIR, EVAL_DIR, PLOT_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
sink(file.path(P6_REP, "heldout_validation_log.txt"), split = TRUE)
on.exit({ while (sink.number() > 0) sink() }, add = TRUE)

cat("=== held-out cross-population validation (G2 & G3) ===\n")
OUTCOME <- "taste_2w_strict"; EPS <- 1e-15; COHORTS <- c("g2", "g3")
best_params <- fread(file.path(P6_REP, "xgb_best_params.csv"))
models <- best_params$model_id

# ---- helpers ----------------------------------------------------------------
clip <- function(p) pmax(pmin(p, 1 - EPS), EPS)
safe_logit <- function(p) qlogis(clip(p)); safe_log <- function(p) log(clip(p))
apply_platt <- function(cal, p_new) as.numeric(predict(cal$model, newdata = data.frame(lp = safe_logit(p_new)), type = "response"))
apply_isotonic <- function(cal, p_new) clip(cal$apply_fun(p_new))
eval_metrics <- function(y, p) {
  if (length(unique(y)) < 2 || length(p) < 10) return(NULL)
  roc_obj <- pROC::roc(y, p, quiet = TRUE, direction = "<")
  auc_ci  <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
  pr_obj  <- PRROC::pr.curve(scores.class0 = p[y == 1], scores.class1 = p[y == 0], curve = FALSE)
  lp <- safe_logit(p)
  list(roc = roc_obj, row = data.table(
    auc = auc_ci[2], auc_lo = auc_ci[1], auc_hi = auc_ci[3], pr_auc = pr_obj$auc.integral,
    brier = mean((p - y)^2), logloss = -mean(y * safe_log(p) + (1 - y) * safe_log(1 - p)),
    cal_intercept = unname(coef(glm(y ~ offset(lp), family = binomial()))[1]),
    cal_slope = unname(coef(glm(y ~ lp, family = binomial()))[2])))
}
sens_at_spec <- function(roc_obj, spec_targets = c(0.80, 0.90, 0.95))
  as.data.table(pROC::coords(roc_obj, x = spec_targets, input = "specificity",
                             ret = c("specificity", "sensitivity", "threshold"), transpose = FALSE))
top_k_capture <- function(y, p, ks = c(0.05, 0.10, 0.20)) {
  out <- list()
  for (k in ks) {
    thr <- quantile(p, 1 - k, na.rm = TRUE); sel <- p >= thr
    out[[as.character(k)]] <- data.table(top_k = k, threshold = as.numeric(thr),
      n_flagged = sum(sel), cases_caught = sum(y == 1 & sel),
      sensitivity = sum(y == 1 & sel) / sum(y == 1), ppv = sum(y == 1 & sel) / max(sum(sel), 1))
  }
  rbindlist(out)
}
build_X <- function(dt, feat_names) {
  missing <- setdiff(feat_names, names(dt))
  if (length(missing) > 0) { for (m in missing) dt[, (m) := NA_real_] }
  X <- as.matrix(dt[, feat_names, with = FALSE]); mode(X) <- "numeric"; X
}

# ---- main loop --------------------------------------------------------------
summary_rows <- list(); thr_rows <- list(); delong_rows <- list()
for (mid in models) {
  cat(strrep("-", 60), "\nModel:", mid, "\n", sep = "")
  boost_pkg <- readRDS(file.path(MODEL_DIR, sprintf("%s_booster.rds", mid)))
  bst <- xgb.load.raw(boost_pkg$booster); feats <- boost_pkg$feature_names
  subset <- best_params[model_id == mid, sample_subset]
  platt_cal <- readRDS(file.path(CAL_DIR, sprintf("%s_platt.rds", mid)))
  iso_path <- file.path(CAL_DIR, sprintf("%s_isotonic.rds", mid))
  iso_cal <- if (file.exists(iso_path)) readRDS(iso_path) else NULL

  g1_oof <- as.data.table(readRDS(file.path(MODEL_DIR, sprintf("%s_oof_calibrated.rds", mid))))
  g1_oof <- g1_oof[!is.na(pred_raw) & !is.na(y)]
  roc_cache <- list(g1 = pROC::roc(g1_oof$y, g1_oof$pred_platt, quiet = TRUE, direction = "<"))

  for (ch in COHORTS) {
    rds <- file.path(ML_DIR, sprintf("group%s_%s.rds", sub("g", "", ch), subset))
    if (!file.exists(rds)) { cat(sprintf("  [%s] MISSING %s - skip\n", ch, rds)); next }
    dt <- as.data.table(readRDS(rds))
    if ("years_baseline_to_taste" %in% feats && !"years_baseline_to_taste" %in% names(dt))
      dt[, years_baseline_to_taste := age - age_baseline]
    dt <- dt[!is.na(get(OUTCOME))]
    X <- build_X(dt, feats); y <- as.integer(dt[[OUTCOME]])
    pred_raw <- predict(bst, xgb.DMatrix(X, missing = NA))
    pred_platt <- apply_platt(platt_cal, pred_raw)
    pred_iso <- if (!is.null(iso_cal)) apply_isotonic(iso_cal, pred_raw) else rep(NA_real_, length(pred_raw))

    eid_col <- intersect(c("eid", "eid_code", "ID"), names(dt))[1]
    preds <- data.table(eid = if (!is.na(eid_col)) dt[[eid_col]] else NA_integer_,
                        y = y, pred_raw = pred_raw, pred_platt = pred_platt, pred_iso = pred_iso)
    saveRDS(preds, file.path(PRED_DIR, sprintf("%s_%s_preds.rds", mid, ch)))
    cat(sprintf("  [%s]  N=%d  cases=%d (%.2f%%)\n", ch, length(y), sum(y), 100 * mean(y)))

    ev <- eval_metrics(y, pred_platt)
    if (is.null(ev)) { cat("    cohort too small / degenerate; skip metrics\n"); next }
    roc_cache[[ch]] <- ev$roc
    row <- ev$row; row[, `:=`(model_id = mid, cohort = ch, variant = "platt", n = length(y), n_cases = sum(y))]
    setcolorder(row, c("model_id", "cohort", "variant", "n", "n_cases"))
    summary_rows[[paste(mid, ch, "platt", sep = "_")]] <- row
    cat(sprintf("    AUC=%.4f [%.4f, %.4f]  PR-AUC=%.4f  CalInt=%.3f  CalSlope=%.3f\n",
                ev$row$auc, ev$row$auc_lo, ev$row$auc_hi, ev$row$pr_auc, ev$row$cal_intercept, ev$row$cal_slope))

    for (v in c("raw", "iso")) {
      p_v <- switch(v, raw = pred_raw, iso = pred_iso)
      if (all(is.na(p_v))) next
      ev_v <- eval_metrics(y, p_v); if (is.null(ev_v)) next
      r_v <- ev_v$row; r_v[, `:=`(model_id = mid, cohort = ch, variant = v, n = length(y), n_cases = sum(y))]
      setcolorder(r_v, c("model_id", "cohort", "variant", "n", "n_cases"))
      summary_rows[[paste(mid, ch, v, sep = "_")]] <- r_v
    }
    ss <- sens_at_spec(ev$roc); ss[, `:=`(model_id = mid, cohort = ch, metric = "sens_at_spec")]
    tk <- top_k_capture(y, pred_platt); tk[, `:=`(model_id = mid, cohort = ch, metric = "top_k_capture")]
    thr_rows[[paste(mid, ch, "ss", sep = "_")]] <- ss; thr_rows[[paste(mid, ch, "tk", sep = "_")]] <- tk

    dt_test <- pROC::roc.test(roc_cache$g1, ev$roc, method = "delong", paired = FALSE)
    delong_rows[[paste(mid, ch, sep = "_")]] <- data.table(model_id = mid, cohort = ch,
      auc_g1 = as.numeric(pROC::auc(roc_cache$g1)), auc_ext = as.numeric(pROC::auc(ev$roc)),
      delta_auc = as.numeric(pROC::auc(ev$roc)) - as.numeric(pROC::auc(roc_cache$g1)),
      z = unname(dt_test$statistic), p_value = dt_test$p.value)
  }

  if (length(roc_cache) >= 2) {
    png(file.path(PLOT_DIR, sprintf("%s_heldout_roc.png", mid)), 1000, 900, res = 150)
    cols <- c(g1 = "#1f77b4", g2 = "#ff7f0e", g3 = "#2ca02c"); lbls <- names(roc_cache)
    plot(roc_cache[[lbls[1]]], col = cols[lbls[1]], lwd = 2, legacy.axes = TRUE,
         main = sprintf("%s - held-out cross-population validation", mid))
    for (l in lbls[-1]) lines(roc_cache[[l]], col = cols[l], lwd = 2)
    legend("bottomright", legend = sprintf("%s (AUC=%.3f, N=%d)", toupper(lbls),
           sapply(roc_cache, function(r) as.numeric(pROC::auc(r))),
           sapply(roc_cache, function(r) length(r$cases) + length(r$controls))),
           col = cols[lbls], lwd = 2, bty = "n", cex = 0.85)
    dev.off()
  }
}

fwrite(rbindlist(summary_rows, fill = TRUE), file.path(EVAL_DIR, "heldout_val_summary.csv"))
fwrite(rbindlist(thr_rows, fill = TRUE), file.path(EVAL_DIR, "heldout_val_thresholded.csv"))
fwrite(rbindlist(delong_rows, fill = TRUE), file.path(EVAL_DIR, "heldout_val_delong_vs_g1.csv"))
cat("\n=== 05 complete ===\nNext: 06_fairness_audit.R\n")
