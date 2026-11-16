# Recalibrate the out-of-fold predictions of each model in group1: Platt scaling, and
# isotonic regression where the analysis sample reaches 20,000. The fitted calibrators are
# reused unchanged on the held-out cohorts.
library(data.table)
library(pROC)
setwd(Sys.getenv("PROJECT_DIR", "."))

MODEL_DIR <- "output/models/xgboost"
CAL_DIR   <- file.path(MODEL_DIR, "calibrators")
REP_DIR   <- "output/model_reports"
EVAL_DIR  <- file.path(REP_DIR, "eval")
dir.create(CAL_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(EVAL_DIR, showWarnings = FALSE, recursive = TRUE)

ISO_MIN_N <- 20000L; EPS <- 1e-15
models <- fread(file.path(REP_DIR, "xgb_best_params.csv"))$model_id

clip       <- function(p) pmax(pmin(p, 1 - EPS), EPS)
safe_logit <- function(p) qlogis(clip(p))
safe_log   <- function(p) log(clip(p))
fit_platt  <- function(p, y) glm(y ~ lp, family = binomial(), data = data.frame(y = y, lp = safe_logit(p)))
apply_platt <- function(fit, p) as.numeric(predict(fit, newdata = data.frame(lp = safe_logit(p)), type = "response"))
fit_isotonic <- function(p, y) {
  ord <- order(p); fit <- isoreg(p[ord], y[ord])
  approxfun(p[ord], fit$yf, method = "linear", rule = 2, ties = "ordered")
}
eval_metrics <- function(y, p, label) {
  lp <- safe_logit(p)
  data.table(variant = label, auc = as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE, direction = "<"))),
             brier = mean((p - y)^2), logloss = -mean(y * safe_log(p) + (1 - y) * safe_log(1 - p)),
             cal_intercept = unname(coef(glm(y ~ offset(lp), family = binomial()))[1]),
             cal_slope = unname(coef(glm(y ~ lp, family = binomial()))[2]))
}

comparison <- list()
for (mid in models) {
  oof <- readRDS(file.path(MODEL_DIR, sprintf("%s_oof_preds.rds", mid)))[!is.na(pred) & !is.na(y)]
  N <- nrow(oof)
  platt_fit <- fit_platt(oof$pred, oof$y)
  oof[, pred_platt := apply_platt(platt_fit, pred)]
  saveRDS(list(method = "platt", model = platt_fit, fit_n = N, fit_cases = sum(oof$y), model_id = mid),
          file.path(CAL_DIR, sprintf("%s_platt.rds", mid)))
  rows <- rbind(eval_metrics(oof$y, oof$pred, "raw"), eval_metrics(oof$y, oof$pred_platt, "platt"))

  if (N >= ISO_MIN_N) {
    iso_fun <- fit_isotonic(oof$pred, oof$y)
    oof[, pred_iso := clip(iso_fun(pred))]
    saveRDS(list(method = "isotonic", apply_fun = iso_fun, fit_n = N, fit_cases = sum(oof$y), model_id = mid),
            file.path(CAL_DIR, sprintf("%s_isotonic.rds", mid)))
    rows <- rbind(rows, eval_metrics(oof$y, oof$pred_iso, "isotonic"))
  } else {
    oof[, pred_iso := NA_real_]
  }
  comparison[[mid]] <- cbind(data.table(model_id = mid, n = N), rows)

  setnames(oof, "pred", "pred_raw")
  setcolorder(oof, c("eid", "fold", "y", "pred_raw", "pred_platt", "pred_iso"))
  saveRDS(oof, file.path(MODEL_DIR, sprintf("%s_oof_calibrated.rds", mid)))
}
fwrite(rbindlist(comparison), file.path(EVAL_DIR, "xgb_calibration_comparison.csv"))
