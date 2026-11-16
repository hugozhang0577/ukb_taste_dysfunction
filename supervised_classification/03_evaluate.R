# Out-of-fold performance in group1 for each model (AUROC with DeLong CI, AUPRC, Brier,
# log loss, calibration intercept and slope, Hosmer-Lemeshow), thresholded metrics, and
# paired DeLong tests between models on their common participants.
library(data.table)
library(pROC)
library(PRROC)
setwd(Sys.getenv("PROJECT_DIR", "."))

MODEL_DIR <- "output/models/xgboost"
REP_DIR   <- "output/model_reports"
EVAL_DIR  <- file.path(REP_DIR, "eval")
dir.create(EVAL_DIR, showWarnings = FALSE, recursive = TRUE)

models <- fread(file.path(REP_DIR, "xgb_best_params.csv"))$model_id
read_oof <- function(mid) readRDS(file.path(MODEL_DIR, sprintf("%s_oof_preds.rds", mid)))
safe_log <- function(p) log(pmax(pmin(p, 1 - 1e-15), 1e-15))

eval_one <- function(y, p) {
  roc_obj <- pROC::roc(y, p, quiet = TRUE, direction = "<")
  auc_ci  <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
  lp      <- qlogis(pmax(pmin(p, 1 - 1e-15), 1e-15))
  # Hosmer-Lemeshow over deciles of predicted risk
  grp <- cut(p, quantile(p, probs = seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)
  hl  <- data.table(grp = grp, y = y, p = p)[, .(obs = sum(y), exp = sum(p), n = .N), by = grp][order(grp)]
  hl_chi <- sum((hl$obs - hl$exp)^2 / (hl$exp * (1 - hl$exp / hl$n)), na.rm = TRUE)
  list(roc = roc_obj, row = data.table(
    auc = auc_ci[2], auc_lo = auc_ci[1], auc_hi = auc_ci[3],
    pr_auc = PRROC::pr.curve(scores.class0 = p[y == 1], scores.class1 = p[y == 0])$auc.integral,
    brier = mean((p - y)^2), logloss = -mean(y * safe_log(p) + (1 - y) * safe_log(1 - p)),
    cal_intercept = unname(coef(glm(y ~ offset(lp), family = binomial()))[1]),
    cal_slope = unname(coef(glm(y ~ lp, family = binomial()))[2]),
    hl_chisq = hl_chi, hl_p = pchisq(hl_chi, df = 8, lower.tail = FALSE)))
}

sens_at_spec <- function(roc_obj, spec = c(0.80, 0.90, 0.95))
  as.data.table(pROC::coords(roc_obj, x = spec, input = "specificity",
                             ret = c("specificity", "sensitivity", "threshold"), transpose = FALSE))

# Treat the top 5%, 10% and 20% of predicted risk as screen-positive
top_k_capture <- function(y, p, ks = c(0.05, 0.10, 0.20))
  rbindlist(lapply(ks, function(k) {
    thr <- quantile(p, 1 - k, na.rm = TRUE); sel <- p >= thr
    data.table(top_k = k, threshold = as.numeric(thr), n_flagged = sum(sel),
               cases_caught = sum(y == 1 & sel), sensitivity = sum(y == 1 & sel) / sum(y == 1),
               ppv = sum(y == 1 & sel) / max(sum(sel), 1))
  }))

summary_rows <- list(); thr_rows <- list()
for (mid in models) {
  oof <- read_oof(mid)[!is.na(pred) & !is.na(y)]
  e <- eval_one(oof$y, oof$pred)
  summary_rows[[mid]] <- data.table(model_id = mid, n = nrow(oof), n_cases = sum(oof$y),
                                    case_rate = mean(oof$y), e$row)
  thr_rows[[length(thr_rows) + 1]] <- sens_at_spec(e$roc)[, `:=`(model_id = mid, metric = "sens_at_spec")]
  thr_rows[[length(thr_rows) + 1]] <- top_k_capture(oof$y, oof$pred)[, `:=`(model_id = mid, metric = "top_k_capture")]
}

pairs <- CJ(a = models, b = models)[a < b]
delong <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  a <- pairs$a[i]; b <- pairs$b[i]
  common <- merge(read_oof(a)[, .(eid, y, p_a = pred)], read_oof(b)[, .(eid, p_b = pred)], by = "eid")
  ra <- pROC::roc(common$y, common$p_a, quiet = TRUE, direction = "<")
  rb <- pROC::roc(common$y, common$p_b, quiet = TRUE, direction = "<")
  tt <- pROC::roc.test(ra, rb, method = "delong")
  data.table(model_a = a, model_b = b, auc_a = as.numeric(pROC::auc(ra)), auc_b = as.numeric(pROC::auc(rb)),
             delta_auc = as.numeric(pROC::auc(rb)) - as.numeric(pROC::auc(ra)),
             z = unname(tt$statistic), p_value = tt$p.value, n_common = nrow(common))
}))

fwrite(rbindlist(summary_rows), file.path(EVAL_DIR, "xgb_eval_summary.csv"))
fwrite(rbindlist(thr_rows, fill = TRUE), file.path(EVAL_DIR, "xgb_eval_thresholded.csv"))
fwrite(delong, file.path(EVAL_DIR, "xgb_delong_pairwise.csv"))
