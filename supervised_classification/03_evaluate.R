#!/usr/bin/env Rscript
# =============================================================================
# 05_evaluate.R  (evaluation)
#
# Consume the OOF predictions from 06b and compute, per tier-model:
#   ROC-AUC (95% DeLong CI), PR-AUC, Brier, log-loss, calibration intercept/slope,
#   Hosmer-Lemeshow p (10 deciles), sensitivity @ specificity {0.80,0.90,0.95},
#   top-k% risk capture, pairwise DeLong tests, and diagnostic ROC/PR/calibration plots.
#
# Input:  models/xgboost/{model_id}_oof_preds.rds
#         model_reports/xgb_best_params.csv (model list)
# Output: model_reports/eval/{xgb_eval_summary, xgb_eval_thresholded, xgb_delong_pairwise}.csv
#         model_reports/eval/plots/*.png + model_reports/evaluate_log.txt
# =============================================================================

rm(list = ls()); gc()
while (sink.number() > 0) sink()
suppressMessages({ library(data.table); library(pROC); library(PRROC) })

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)

MODEL_DIR <- "output/models/xgboost"
P6_REP    <- "output/model_reports"
EVAL_DIR  <- file.path(P6_REP, "eval")
PLOT_DIR  <- file.path(EVAL_DIR, "plots")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)
sink(file.path(P6_REP, "evaluate_log.txt"), split = TRUE)
on.exit({ while (sink.number() > 0) sink() }, add = TRUE)

cat("=== XGBoost evaluation (OOF) ===\n")
models <- fread(file.path(P6_REP, "xgb_best_params.csv"))$model_id
cat("models:", paste(models, collapse = ", "), "\n")

# ---- helpers ----------------------------------------------------------------
safe_log <- function(p) log(pmax(pmin(p, 1 - 1e-15), 1e-15))

eval_one <- function(oof) {
  y <- oof$y; p <- oof$pred
  roc_obj <- pROC::roc(y, p, quiet = TRUE, direction = "<")
  auc_ci  <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
  pr_obj  <- PRROC::pr.curve(scores.class0 = p[y == 1], scores.class1 = p[y == 0], curve = TRUE)
  brier   <- mean((p - y)^2)
  logloss <- -mean(y * safe_log(p) + (1 - y) * safe_log(1 - p))
  lp      <- qlogis(pmax(pmin(p, 1 - 1e-15), 1e-15))
  cal_int <- coef(glm(y ~ offset(lp), family = binomial()))[1]
  cal_slp <- coef(glm(y ~ lp, family = binomial()))[2]
  grp     <- cut(p, quantile(p, probs = seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)
  hl_tab  <- data.table(grp = grp, y = y, p = p)[, .(obs = sum(y), exp = sum(p), n = .N), by = grp][order(grp)]
  hl_tab[, stat := (obs - exp)^2 / (exp * (1 - exp / n))]
  hl_chi  <- sum(hl_tab$stat, na.rm = TRUE); hl_p <- pchisq(hl_chi, df = 8, lower.tail = FALSE)
  list(auc = auc_ci[2], auc_lo = auc_ci[1], auc_hi = auc_ci[3], pr_auc = pr_obj$auc.integral,
       brier = brier, logloss = logloss, cal_intercept = unname(cal_int), cal_slope = unname(cal_slp),
       hl_chisq = hl_chi, hl_p = hl_p, roc_obj = roc_obj, pr_obj = pr_obj)
}

sens_at_spec <- function(roc_obj, spec_targets = c(0.80, 0.90, 0.95)) {
  as.data.table(pROC::coords(roc_obj, x = spec_targets, input = "specificity",
                             ret = c("specificity", "sensitivity", "threshold"), transpose = FALSE))
}

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

# ---- per-model evaluation ---------------------------------------------------
summary_rows <- list(); thr_rows <- list(); roc_objs <- list()
for (mid in models) {
  cat(strrep("-", 60), "\nEvaluating:", mid, "\n", sep = "")
  oof <- as.data.table(readRDS(file.path(MODEL_DIR, sprintf("%s_oof_preds.rds", mid))))
  oof <- oof[!is.na(pred) & !is.na(y)]
  e <- eval_one(oof); roc_objs[[mid]] <- e$roc_obj
  cat(sprintf("  N=%d  cases=%d  AUC=%.4f [%.4f, %.4f]  PR-AUC=%.4f  Brier=%.4f  HL p=%.3g\n",
              nrow(oof), sum(oof$y), e$auc, e$auc_lo, e$auc_hi, e$pr_auc, e$brier, e$hl_p))
  summary_rows[[mid]] <- data.table(model_id = mid, n = nrow(oof), n_cases = sum(oof$y),
    case_rate = mean(oof$y), auc = e$auc, auc_lo = e$auc_lo, auc_hi = e$auc_hi, pr_auc = e$pr_auc,
    brier = e$brier, logloss = e$logloss, cal_intercept = e$cal_intercept, cal_slope = e$cal_slope,
    hl_chisq = e$hl_chisq, hl_p = e$hl_p)
  ss <- sens_at_spec(e$roc_obj); ss[, `:=`(model_id = mid, metric = "sens_at_spec")]
  tk <- top_k_capture(oof$y, oof$pred); tk[, `:=`(model_id = mid, metric = "top_k_capture")]
  thr_rows[[paste0(mid, "_ss")]] <- ss; thr_rows[[paste0(mid, "_tk")]] <- tk

  png(file.path(PLOT_DIR, sprintf("%s_roc.png", mid)), 900, 900, res = 150)
  plot(e$roc_obj, main = sprintf("%s  AUC=%.3f [%.3f,%.3f]", mid, e$auc, e$auc_lo, e$auc_hi), legacy.axes = TRUE)
  dev.off()
  png(file.path(PLOT_DIR, sprintf("%s_pr.png", mid)), 900, 900, res = 150)
  plot(e$pr_obj$curve[, 1], e$pr_obj$curve[, 2], type = "l", xlab = "Recall", ylab = "Precision",
       main = sprintf("%s  PR-AUC=%.3f (baseline=%.3f)", mid, e$pr_auc, mean(oof$y)))
  abline(h = mean(oof$y), lty = 2, col = "grey50"); dev.off()
  cg <- cut(oof$pred, quantile(oof$pred, probs = seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)
  cal <- oof[, .(mean_pred = mean(pred), obs_rate = mean(y), n = .N), by = cg][order(cg)]
  png(file.path(PLOT_DIR, sprintf("%s_calibration.png", mid)), 900, 900, res = 150)
  plot(cal$mean_pred, cal$obs_rate, pch = 19, xlim = c(0, max(cal$mean_pred, cal$obs_rate)),
       ylim = c(0, max(cal$mean_pred, cal$obs_rate)), xlab = "Mean predicted risk (decile)",
       ylab = "Observed event rate",
       main = sprintf("%s calibration  int=%.2f slope=%.2f", mid, e$cal_intercept, e$cal_slope))
  abline(0, 1, lty = 2, col = "grey40"); dev.off()
}

# ---- pairwise DeLong --------------------------------------------------------
cat("\npairwise DeLong tests\n")
pairs <- CJ(a = models, b = models)[a < b]; delong_rows <- list()
for (i in seq_len(nrow(pairs))) {
  a <- pairs$a[i]; b <- pairs$b[i]
  oa <- as.data.table(readRDS(file.path(MODEL_DIR, sprintf("%s_oof_preds.rds", a))))
  ob <- as.data.table(readRDS(file.path(MODEL_DIR, sprintf("%s_oof_preds.rds", b))))
  common <- merge(oa[, .(eid, y_a = y, p_a = pred)], ob[, .(eid, y_b = y, p_b = pred)], by = "eid")
  stopifnot(all(common$y_a == common$y_b))
  ra <- pROC::roc(common$y_a, common$p_a, quiet = TRUE, direction = "<")
  rb <- pROC::roc(common$y_b, common$p_b, quiet = TRUE, direction = "<")
  dt <- pROC::roc.test(ra, rb, method = "delong")
  delong_rows[[i]] <- data.table(model_a = a, model_b = b,
    auc_a = as.numeric(pROC::auc(ra)), auc_b = as.numeric(pROC::auc(rb)),
    delta_auc = as.numeric(pROC::auc(rb)) - as.numeric(pROC::auc(ra)),
    z = unname(dt$statistic), p_value = dt$p.value, n_common = nrow(common))
}
delong_dt <- rbindlist(delong_rows)

# ---- combined ROC -----------------------------------------------------------
png(file.path(PLOT_DIR, "combined_roc.png"), 1200, 1000, res = 150)
cols <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b")
plot(roc_objs[[models[1]]], col = cols[1], lwd = 2, legacy.axes = TRUE, main = "XGBoost OOF ROC - all tier models")
for (i in seq_along(models)[-1]) lines(roc_objs[[models[i]]], col = cols[i], lwd = 2)
legend("bottomright", legend = sprintf("%s (AUC=%.3f)", models,
       sapply(roc_objs, function(r) as.numeric(pROC::auc(r)))),
       col = cols[seq_along(models)], lwd = 2, bty = "n", cex = 0.8)
dev.off()

fwrite(rbindlist(summary_rows), file.path(EVAL_DIR, "xgb_eval_summary.csv"))
fwrite(rbindlist(thr_rows, fill = TRUE), file.path(EVAL_DIR, "xgb_eval_thresholded.csv"))
fwrite(delong_dt, file.path(EVAL_DIR, "xgb_delong_pairwise.csv"))
cat("\n=== 03 complete ===\n")
