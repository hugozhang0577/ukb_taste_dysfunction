# Apply the group1 models and calibrators unchanged (no retraining) to the held-out
# cohorts group2 and group3; compare AUROC with the group1 out-of-fold AUROC (unpaired DeLong).
library(data.table)
library(xgboost)
library(pROC)
library(PRROC)
setwd(Sys.getenv("PROJECT_DIR", "."))

ML_DIR    <- "output/ml_ready"
MODEL_DIR <- "output/models/xgboost"
CAL_DIR   <- file.path(MODEL_DIR, "calibrators")
PRED_DIR  <- file.path(MODEL_DIR, "predictions")
EVAL_DIR  <- "output/model_reports/eval"
for (d in c(PRED_DIR, EVAL_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

OUTCOME <- "taste_2w_strict"; EPS <- 1e-15; COHORTS <- c("g2", "g3")
best_params <- fread("output/model_reports/xgb_best_params.csv")

clip       <- function(p) pmax(pmin(p, 1 - EPS), EPS)
safe_logit <- function(p) qlogis(clip(p)); safe_log <- function(p) log(clip(p))
eval_metrics <- function(y, p) {
  if (length(unique(y)) < 2 || length(p) < 10) return(NULL)
  roc_obj <- pROC::roc(y, p, quiet = TRUE, direction = "<")
  auc_ci  <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
  lp <- safe_logit(p)
  list(roc = roc_obj, row = data.table(
    auc = auc_ci[2], auc_lo = auc_ci[1], auc_hi = auc_ci[3],
    pr_auc = PRROC::pr.curve(scores.class0 = p[y == 1], scores.class1 = p[y == 0])$auc.integral,
    brier = mean((p - y)^2), logloss = -mean(y * safe_log(p) + (1 - y) * safe_log(1 - p)),
    cal_intercept = unname(coef(glm(y ~ offset(lp), family = binomial()))[1]),
    cal_slope = unname(coef(glm(y ~ lp, family = binomial()))[2])))
}
sens_at_spec <- function(roc_obj, spec = c(0.80, 0.90, 0.95))
  as.data.table(pROC::coords(roc_obj, x = spec, input = "specificity",
                             ret = c("specificity", "sensitivity", "threshold"), transpose = FALSE))
top_k_capture <- function(y, p, ks = c(0.05, 0.10, 0.20))
  rbindlist(lapply(ks, function(k) {
    thr <- quantile(p, 1 - k, na.rm = TRUE); sel <- p >= thr
    data.table(top_k = k, threshold = as.numeric(thr), n_flagged = sum(sel),
               cases_caught = sum(y == 1 & sel), sensitivity = sum(y == 1 & sel) / sum(y == 1),
               ppv = sum(y == 1 & sel) / max(sum(sel), 1))
  }))
# Features absent from a held-out cohort enter as missing
build_X <- function(dt, feats) {
  for (m in setdiff(feats, names(dt))) dt[, (m) := NA_real_]
  X <- as.matrix(dt[, feats, with = FALSE]); mode(X) <- "numeric"; X
}

summary_rows <- list(); thr_rows <- list(); delong_rows <- list()
for (mid in best_params$model_id) {
  pkg   <- readRDS(file.path(MODEL_DIR, sprintf("%s_booster.rds", mid)))
  bst   <- xgb.load.raw(pkg$booster); feats <- pkg$feature_names
  platt <- readRDS(file.path(CAL_DIR, sprintf("%s_platt.rds", mid)))
  iso_path <- file.path(CAL_DIR, sprintf("%s_isotonic.rds", mid))
  iso   <- if (file.exists(iso_path)) readRDS(iso_path) else NULL   # fitted only when N >= 20,000
  g1    <- readRDS(file.path(MODEL_DIR, sprintf("%s_oof_calibrated.rds", mid)))[!is.na(pred_raw) & !is.na(y)]
  roc_g1 <- pROC::roc(g1$y, g1$pred_platt, quiet = TRUE, direction = "<")
  subset <- best_params[model_id == mid, sample_subset]

  for (ch in COHORTS) {
    dt <- readRDS(file.path(ML_DIR, sprintf("group%s_%s.rds", sub("g", "", ch), subset)))
    dt <- dt[!is.na(get(OUTCOME))]
    y  <- as.integer(dt[[OUTCOME]])
    pred_raw   <- predict(bst, xgb.DMatrix(build_X(dt, feats), missing = NA))
    pred_platt <- as.numeric(predict(platt$model, newdata = data.frame(lp = safe_logit(pred_raw)), type = "response"))
    pred_iso   <- if (!is.null(iso)) clip(iso$apply_fun(pred_raw)) else rep(NA_real_, length(pred_raw))
    saveRDS(data.table(eid = dt$eid, y = y, pred_raw = pred_raw, pred_platt = pred_platt, pred_iso = pred_iso),
            file.path(PRED_DIR, sprintf("%s_%s_preds.rds", mid, ch)))

    ev <- eval_metrics(y, pred_platt)
    if (is.null(ev)) next
    for (v in c("platt", "raw", "iso")) {
      p_v  <- switch(v, platt = pred_platt, raw = pred_raw, iso = pred_iso)
      ev_v <- if (v == "platt") ev else if (all(is.na(p_v))) NULL else eval_metrics(y, p_v)
      if (is.null(ev_v)) next
      summary_rows[[paste(mid, ch, v)]] <- cbind(
        data.table(model_id = mid, cohort = ch, variant = v, n = length(y), n_cases = sum(y)), ev_v$row)
    }
    thr_rows[[paste(mid, ch, "ss")]] <- sens_at_spec(ev$roc)[, `:=`(model_id = mid, cohort = ch, metric = "sens_at_spec")]
    thr_rows[[paste(mid, ch, "tk")]] <- top_k_capture(y, pred_platt)[, `:=`(model_id = mid, cohort = ch, metric = "top_k_capture")]

    tt <- pROC::roc.test(roc_g1, ev$roc, method = "delong", paired = FALSE)
    delong_rows[[paste(mid, ch)]] <- data.table(model_id = mid, cohort = ch,
      auc_g1 = as.numeric(pROC::auc(roc_g1)), auc_heldout = as.numeric(pROC::auc(ev$roc)),
      delta_auc = as.numeric(pROC::auc(ev$roc)) - as.numeric(pROC::auc(roc_g1)),
      z = unname(tt$statistic), p_value = tt$p.value)
  }
}

fwrite(rbindlist(summary_rows, fill = TRUE), file.path(EVAL_DIR, "heldout_val_summary.csv"))
fwrite(rbindlist(thr_rows, fill = TRUE), file.path(EVAL_DIR, "heldout_val_thresholded.csv"))
fwrite(rbindlist(delong_rows, fill = TRUE), file.path(EVAL_DIR, "heldout_val_delong_vs_g1.csv"))
