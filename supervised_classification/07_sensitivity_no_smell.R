# Sensitivity analysis: retrain each of the six models without the four smell items
# (smell_any, smell_time, smell_extent, smell_2w_strict_OR) using the same CV design, apply
# the variants to group2/group3, and compare AUROC with the original model (paired DeLong).
library(data.table)
library(xgboost)
library(pROC)
library(PRROC)
set.seed(42)
setwd(Sys.getenv("PROJECT_DIR", "."))

ML_DIR   <- "output/ml_ready"
XGB_DIR  <- "output/models/xgboost"
OUT_DIR  <- file.path(XGB_DIR, "sensitivity_no_smell")
PRED_DIR <- file.path(OUT_DIR, "predictions")
EVAL_DIR <- "output/model_reports/eval"
for (d in c(PRED_DIR, EVAL_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

OUTCOME <- "taste_2w_strict"; N_FOLDS <- 5; NTHREAD <- 4
MAX_ROUNDS <- 2000; EARLY_STOP <- 50; EPS <- 1e-15
# Analysis samples below 20,000 use a reduced grid
SMALL_N_THR <- 20000L
GRID_FULL  <- CJ(max_depth = c(3, 5, 7), eta = c(0.03, 0.05, 0.1), subsample = 0.8, colsample_bytree = 0.8, min_child_weight = 5)
GRID_SMALL <- CJ(max_depth = c(3, 5), eta = c(0.03, 0.05), subsample = 0.8, colsample_bytree = 0.8, min_child_weight = 5)

clip <- function(p) pmax(pmin(p, 1 - EPS), EPS)
safe_logit <- function(p) qlogis(clip(p)); safe_log <- function(p) log(clip(p))
make_folds <- function(y, k, seed = 42) {
  set.seed(seed)
  fold <- integer(length(y))
  fold[y == 1] <- sample(rep_len(seq_len(k), sum(y == 1)))
  fold[y == 0] <- sample(rep_len(seq_len(k), sum(y == 0)))
  fold
}
# Same CV design as 02_train_xgboost.R
cv_one <- function(X, y, folds, params) {
  oof <- rep(NA_real_, length(y)); metrics <- list()
  for (k in seq_len(max(folds))) {
    tr <- which(folds != k); va <- which(folds == k)
    es <- tr[make_folds(y[tr], N_FOLDS, seed = 42 + k) == 1]; fi <- setdiff(tr, es)
    spw <- sum(y[fi] == 0) / max(sum(y[fi] == 1), 1)
    dtr <- xgb.DMatrix(X[fi, , drop = FALSE], label = y[fi], missing = NA)
    des <- xgb.DMatrix(X[es, , drop = FALSE], label = y[es], missing = NA,
                       weight = ifelse(y[es] == 1, spw, 1))
    dva <- xgb.DMatrix(X[va, , drop = FALSE], label = y[va], missing = NA)
    prm <- c(params, list(objective = "binary:logistic", eval_metric = "auc", eval_metric = "logloss",
                          scale_pos_weight = spw, nthread = NTHREAD))
    fit <- xgb.train(params = prm, data = dtr, nrounds = MAX_ROUNDS, watchlist = list(es = des),
                     early_stopping_rounds = EARLY_STOP, verbose = 0)
    pr <- predict(fit, dva, iterationrange = c(1, fit$best_iteration))
    oof[va] <- pr
    metrics[[k]] <- data.table(auc = as.numeric(pROC::auc(y[va], pr, quiet = TRUE)),
                               best_iter = fit$best_iteration,
                               es_auc = fit$evaluation_log[iter == fit$best_iteration, es_auc])
  }
  metrics <- rbindlist(metrics)
  list(oof = oof, metrics = metrics, mean_best_iter = max(1L, round(mean(metrics$best_iter))))
}
build_X <- function(dt, feats) {
  for (m in setdiff(feats, names(dt))) dt[, (m) := NA_real_]
  X <- as.matrix(dt[, feats, with = FALSE]); mode(X) <- "numeric"; X
}
eval_metrics <- function(y, p) {
  if (length(unique(y)) < 2 || sum(y) < 10) return(NULL)
  roc_obj <- pROC::roc(y, p, quiet = TRUE, direction = "<")
  auc_ci  <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
  lp <- safe_logit(p)
  data.table(auc = auc_ci[2], auc_lo = auc_ci[1], auc_hi = auc_ci[3],
    pr_auc = PRROC::pr.curve(scores.class0 = p[y == 1], scores.class1 = p[y == 0])$auc.integral,
    brier = mean((p - y)^2), logloss = -mean(y * safe_log(p) + (1 - y) * safe_log(1 - p)),
    cal_intercept = unname(coef(glm(y ~ offset(lp), family = binomial()))[1]),
    cal_slope = unname(coef(glm(y ~ lp, family = binomial()))[2]))
}
# Paired DeLong on the participants scored by both the variant and the original model
delong_vs_original <- function(v, o, var_id, orig_id, ch) {
  m <- merge(v[, .(eid, y_a = y, p_a = pred_platt)], o[, .(eid, y_b = y, p_b = pred_platt)], by = "eid")
  if (nrow(m) < 50 || sum(m$y_a) < 10) return(NULL)
  ra <- pROC::roc(m$y_a, m$p_a, quiet = TRUE, direction = "<")
  rb <- pROC::roc(m$y_b, m$p_b, quiet = TRUE, direction = "<")
  tt <- pROC::roc.test(ra, rb, method = "delong", paired = TRUE)
  data.table(variant_id = var_id, original_id = orig_id, cohort = ch,
             AUC_no_smell = as.numeric(pROC::auc(ra)), AUC_original = as.numeric(pROC::auc(rb)),
             delta_AUC = as.numeric(pROC::auc(ra)) - as.numeric(pROC::auc(rb)),
             z = unname(tt$statistic), p_value = tt$p.value, n_common = nrow(m), paired = TRUE)
}

recipe <- fread("output/feature_manifest/model_feature_recipes.csv")
best_params_orig <- fread("output/model_reports/xgb_best_params.csv")
summary_rows <- list(); delong_rows <- list(); feat_count_rows <- list()

for (orig_id in paste0("model", 1:6)) {
  var_id <- paste0(orig_id, "_no_smell")
  feats_full  <- recipe[model_id == orig_id, feature_id]
  smell_feats <- grep("^smell", feats_full, ignore.case = TRUE, value = TRUE)
  feats <- setdiff(feats_full, smell_feats)
  feat_count_rows[[var_id]] <- data.table(original_id = orig_id, variant_id = var_id,
    n_feat_original = length(feats_full), n_feat_after_smell_strip = length(feats),
    n_smell_removed = length(smell_feats), smell_feats_removed = paste(smell_feats, collapse = ";"))

  subset <- best_params_orig[model_id == orig_id, sample_subset]
  dt <- readRDS(file.path(ML_DIR, sprintf("group1_%s.rds", subset)))[!is.na(get(OUTCOME))]
  X <- build_X(dt, feats); y <- as.integer(dt[[OUTCOME]]); N <- length(y); n_pos <- sum(y)
  grid  <- if (N < SMALL_N_THR) GRID_SMALL else GRID_FULL
  folds <- make_folds(y, N_FOLDS)

  grid_dt <- rbindlist(lapply(seq_len(nrow(grid)), function(g) {
    cv <- cv_one(X, y, folds, as.list(grid[g]))
    data.table(grid[g], mean_auc = mean(cv$metrics$auc), mean_es_auc = mean(cv$metrics$es_auc))
  }))
  best <- grid_dt[which.max(mean_es_auc)]   # chosen on the inner split, not the outer folds
  best_params <- as.list(best[, .(max_depth, eta, subsample, colsample_bytree, min_child_weight)])
  final_cv <- cv_one(X, y, folds, best_params)
  oof_dt <- data.table(eid = dt$eid, fold = folds, y = y, pred_raw = final_cv$oof)

  final_params <- c(best_params, list(objective = "binary:logistic", eval_metric = "auc",
                                      scale_pos_weight = (N - n_pos) / max(n_pos, 1), nthread = NTHREAD))
  final_fit <- xgb.train(params = final_params, data = xgb.DMatrix(X, label = y, missing = NA),
                         nrounds = final_cv$mean_best_iter, verbose = 0)
  saveRDS(list(booster = xgb.save.raw(final_fit), feature_names = feats, params = final_params,
               nrounds = final_cv$mean_best_iter, n_train = N, n_cases = n_pos,
               variant_id = var_id, original_id = orig_id),
          file.path(OUT_DIR, sprintf("%s_booster.rds", var_id)))
  fwrite(xgb.importance(feature_names = feats, model = final_fit),
         file.path(OUT_DIR, sprintf("%s_importance.csv", var_id)))

  platt_fit <- glm(y ~ lp_oof, family = binomial(),
                   data = data.frame(y = oof_dt$y, lp_oof = safe_logit(oof_dt$pred_raw)))
  apply_platt <- function(p) as.numeric(predict(platt_fit, newdata = data.frame(lp_oof = safe_logit(p)), type = "response"))
  oof_dt[, pred_platt := apply_platt(pred_raw)]
  saveRDS(oof_dt, file.path(OUT_DIR, sprintf("%s_oof_calibrated.rds", var_id)))

  ev <- eval_metrics(oof_dt$y, oof_dt$pred_platt)
  if (!is.null(ev)) summary_rows[[paste(var_id, "G1")]] <- cbind(
    data.table(variant_id = var_id, original_id = orig_id, cohort = "G1", n = N, n_cases = n_pos), ev)
  delong_rows[[paste(var_id, "G1")]] <- delong_vs_original(
    oof_dt, readRDS(file.path(XGB_DIR, sprintf("%s_oof_calibrated.rds", orig_id))), var_id, orig_id, "G1")

  for (gnum in 2:3) {
    ch <- paste0("G", gnum)
    d2 <- readRDS(file.path(ML_DIR, sprintf("group%d_%s.rds", gnum, subset)))[!is.na(get(OUTCOME))]
    if (nrow(d2) < 100) next
    y2 <- as.integer(d2[[OUTCOME]])
    raw2 <- predict(final_fit, xgb.DMatrix(build_X(d2, feats), missing = NA))
    preds <- data.table(eid = d2$eid, y = y2, pred_raw = raw2, pred_platt = apply_platt(raw2))
    saveRDS(preds, file.path(PRED_DIR, sprintf("%s_%s_preds.rds", var_id, tolower(ch))))
    ev <- eval_metrics(y2, preds$pred_platt)
    if (!is.null(ev)) summary_rows[[paste(var_id, ch)]] <- cbind(
      data.table(variant_id = var_id, original_id = orig_id, cohort = ch, n = length(y2), n_cases = sum(y2)), ev)
    delong_rows[[paste(var_id, ch)]] <- delong_vs_original(
      preds, readRDS(file.path(XGB_DIR, "predictions", sprintf("%s_%s_preds.rds", orig_id, tolower(ch)))),
      var_id, orig_id, ch)
  }
}

fwrite(rbindlist(summary_rows, fill = TRUE), file.path(EVAL_DIR, "sensitivity_no_smell_summary.csv"))
fwrite(rbindlist(delong_rows, fill = TRUE), file.path(EVAL_DIR, "sensitivity_no_smell_delong_vs_original.csv"))
fwrite(rbindlist(feat_count_rows, fill = TRUE), file.path(EVAL_DIR, "sensitivity_no_smell_feature_counts.csv"))
