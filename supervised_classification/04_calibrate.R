#!/usr/bin/env Rscript
# =============================================================================
# 06_calibrate.R  (probability recalibration)
#
# 06b uses scale_pos_weight to handle class imbalance, which preserves ranking
# (AUC) but distorts raw probabilities. Refit two calibrators on the G1 OOF
# predictions and save them for reuse on G2/G3 in 06i:
#   Platt scaling   - glm(y ~ qlogis(p_raw), binomial); parametric, stable
#   Isotonic        - nonparametric monotone fit (only if N >= 20K)
#
# Input:  models/xgboost/{model_id}_oof_preds.rds
# Output: models/xgboost/calibrators/{model_id}_{platt,isotonic}.rds
#         models/xgboost/{model_id}_oof_calibrated.rds
#         model_reports/eval/xgb_calibration_comparison.csv + plots + calibrate_log.txt
# =============================================================================

rm(list = ls()); gc()
while (sink.number() > 0) sink()
suppressMessages({ library(data.table); library(pROC) })

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)

MODEL_DIR <- "output/models/xgboost"
CAL_DIR   <- file.path(MODEL_DIR, "calibrators")
P6_REP    <- "output/model_reports"
EVAL_DIR  <- file.path(P6_REP, "eval")
PLOT_DIR  <- file.path(EVAL_DIR, "plots")
dir.create(CAL_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)
sink(file.path(P6_REP, "calibrate_log.txt"), split = TRUE)
on.exit({ while (sink.number() > 0) sink() }, add = TRUE)

cat("=== recalibration (Platt + Isotonic) ===\n")
ISO_MIN_N <- 20000L; EPS <- 1e-15
models <- fread(file.path(P6_REP, "xgb_best_params.csv"))$model_id
cat("models:", paste(models, collapse = ", "), "\n")

# ---- helpers ----------------------------------------------------------------
clip <- function(p) pmax(pmin(p, 1 - EPS), EPS)
safe_logit <- function(p) qlogis(clip(p))
safe_log   <- function(p) log(clip(p))
fit_platt <- function(p_raw, y) glm(y ~ lp, family = binomial(), data = data.frame(y = y, lp = safe_logit(p_raw)))
apply_platt <- function(fit, p_new) as.numeric(predict(fit, newdata = data.frame(lp = safe_logit(p_new)), type = "response"))
fit_isotonic <- function(p_raw, y) {
  ord <- order(p_raw); fit <- isoreg(p_raw[ord], y[ord])
  approxfun(p_raw[ord], fit$yf, method = "linear", rule = 2, ties = "ordered")
}
apply_isotonic <- function(iso_fun, p_new) clip(iso_fun(p_new))
eval_metrics <- function(y, p, label) {
  auc <- as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE, direction = "<")))
  lp <- safe_logit(p)
  data.table(variant = label, auc = auc, brier = mean((p - y)^2),
             logloss = -mean(y * safe_log(p) + (1 - y) * safe_log(1 - p)),
             cal_intercept = unname(coef(glm(y ~ offset(lp), family = binomial()))[1]),
             cal_slope = unname(coef(glm(y ~ lp, family = binomial()))[2]))
}
reliability <- function(p, y) {
  q <- unique(quantile(p, probs = seq(0, 1, 0.1), na.rm = TRUE))
  if (length(q) < 3) return(data.table(mean_pred = p, obs_rate = y))
  g <- cut(p, q, include.lowest = TRUE, labels = FALSE)
  data.table(p = p, y = y, g = g)[, .(mean_pred = mean(p), obs_rate = mean(y)), by = g][order(g)]
}

# ---- main loop --------------------------------------------------------------
comparison_rows <- list()
for (mid in models) {
  cat(strrep("-", 60), "\nCalibrating:", mid, "\n", sep = "")
  oof <- as.data.table(readRDS(file.path(MODEL_DIR, sprintf("%s_oof_preds.rds", mid))))
  oof <- oof[!is.na(pred) & !is.na(y)]; N <- nrow(oof)
  cat(sprintf("  N=%d  cases=%d (%.2f%%)\n", N, sum(oof$y), 100 * mean(oof$y)))

  platt_fit <- fit_platt(oof$pred, oof$y)
  oof[, pred_platt := apply_platt(platt_fit, pred)]
  saveRDS(list(method = "platt", model = platt_fit, fit_n = N, fit_cases = sum(oof$y), model_id = mid),
          file.path(CAL_DIR, sprintf("%s_platt.rds", mid)))

  use_iso <- N >= ISO_MIN_N
  if (use_iso) {
    iso_fun <- fit_isotonic(oof$pred, oof$y)
    oof[, pred_iso := apply_isotonic(iso_fun, pred)]
    saveRDS(list(method = "isotonic", apply_fun = iso_fun, fit_n = N, fit_cases = sum(oof$y), model_id = mid),
            file.path(CAL_DIR, sprintf("%s_isotonic.rds", mid)))
  } else {
    oof[, pred_iso := NA_real_]
    cat(sprintf("  [isotonic] skipped (N=%d < %d)\n", N, ISO_MIN_N))
  }

  rows <- rbind(eval_metrics(oof$y, oof$pred, "raw"), eval_metrics(oof$y, oof$pred_platt, "platt"))
  if (use_iso) rows <- rbind(rows, eval_metrics(oof$y, oof$pred_iso, "isotonic"))
  rows[, `:=`(model_id = mid, n = N)]; setcolorder(rows, c("model_id", "n", "variant"))
  comparison_rows[[mid]] <- rows; print(rows)

  setcolorder(oof, c("eid", "fold", "y", "pred", "pred_platt", "pred_iso"))
  setnames(oof, "pred", "pred_raw")
  saveRDS(oof, file.path(MODEL_DIR, sprintf("%s_oof_calibrated.rds", mid)))

  png(file.path(PLOT_DIR, sprintf("%s_calibration_comparison.png", mid)), 1000, 900, res = 150)
  rel_raw <- reliability(oof$pred_raw, oof$y); rel_platt <- reliability(oof$pred_platt, oof$y)
  lim <- c(0, max(rel_raw$mean_pred, rel_platt$mean_pred, rel_raw$obs_rate, rel_platt$obs_rate, 0.5, na.rm = TRUE))
  plot(rel_raw$mean_pred, rel_raw$obs_rate, type = "b", pch = 19, col = "#d62728", xlim = lim, ylim = lim,
       xlab = "Mean predicted risk (decile)", ylab = "Observed event rate",
       main = sprintf("%s calibration (raw vs Platt%s)", mid, if (use_iso) " vs Isotonic" else ""))
  abline(0, 1, lty = 2, col = "grey40")
  lines(rel_platt$mean_pred, rel_platt$obs_rate, type = "b", pch = 19, col = "#1f77b4")
  if (use_iso) {
    rel_iso <- reliability(oof$pred_iso, oof$y)
    lines(rel_iso$mean_pred, rel_iso$obs_rate, type = "b", pch = 19, col = "#2ca02c")
    legend("topleft", c("Raw", "Platt", "Isotonic", "Ideal"), col = c("#d62728", "#1f77b4", "#2ca02c", "grey40"),
           lty = c(1, 1, 1, 2), pch = c(19, 19, 19, NA), bty = "n")
  } else {
    legend("topleft", c("Raw", "Platt", "Ideal"), col = c("#d62728", "#1f77b4", "grey40"),
           lty = c(1, 1, 2), pch = c(19, 19, NA), bty = "n")
  }
  dev.off()
}

fwrite(rbindlist(comparison_rows), file.path(EVAL_DIR, "xgb_calibration_comparison.csv"))
cat("\n=== 04 complete ===\nNext: 05_heldout_validation.R\n")
