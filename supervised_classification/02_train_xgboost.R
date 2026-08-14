#!/usr/bin/env Rscript
# =============================================================================
# 04_train_xgboost.R  (primary algorithm)
#
# Train one XGBoost classifier per tier-model on Group 1 only. Save the booster,
# out-of-fold predictions (with eid), CV metrics, best hyperparameters, and
# gain-based feature importance.
#
# Input:  feature_manifest/model_feature_recipes.csv (from 06a)
#         data/ml_ready/group1_{full,olink,nmr,olink_nmr}.rds (immutable)
# Output: models/xgboost/{model_id}_{booster,oof_preds}.rds + _importance.csv
#         model_reports/xgb_cv_metrics.csv + xgb_best_params.csv + train_xgboost_log.txt
# =============================================================================

rm(list = ls()); gc()
while (sink.number() > 0) sink()
suppressMessages({ library(data.table); library(xgboost); library(pROC) })
set.seed(42)

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)

MF_DIR     <- "output/feature_manifest"
ML_DIR     <- "output/ml_ready"
MODEL_DIR  <- "output/models/xgboost"
P6_REP_DIR <- "output/model_reports"
dir.create(MODEL_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(P6_REP_DIR, showWarnings = FALSE, recursive = TRUE)

sink(file.path(P6_REP_DIR, "train_xgboost_log.txt"), split = TRUE)
on.exit({ while (sink.number() > 0) sink() }, add = TRUE)

cat("=== XGBoost training ===\n")

# ---- configuration ----------------------------------------------------------
OUTCOME    <- "taste_2w_strict"
N_FOLDS    <- 5
NTHREAD    <- 4              # fixed for cross-platform reproducibility
MAX_ROUNDS <- 2000
EARLY_STOP <- 50

# Initial sweep on max_depth x eta (other hyperparameters at small-sample /
# imbalanced defaults).
GRID <- CJ(max_depth = c(3, 5, 7), eta = c(0.03, 0.05, 0.1),
           subsample = 0.8, colsample_bytree = 0.8, min_child_weight = 5)
cat("hyperparameter grid:", nrow(GRID), "combinations\n")

recipe <- fread(file.path(MF_DIR, "model_feature_recipes.csv"))
models <- recipe[, unique(model_id)]
cat("models to train:", paste(models, collapse = ", "), "\n\n")

# ---- helpers ----------------------------------------------------------------
make_folds <- function(y, k, seed = 42) {
  set.seed(seed)
  idx_pos <- which(y == 1); idx_neg <- which(y == 0)
  fold <- integer(length(y))
  fold[idx_pos] <- sample(rep_len(seq_len(k), length(idx_pos)))
  fold[idx_neg] <- sample(rep_len(seq_len(k), length(idx_neg)))
  fold
}

cv_one <- function(X, y, folds, params, nrounds = MAX_ROUNDS, early_stop = EARLY_STOP) {
  oof_pred <- rep(NA_real_, length(y)); fold_metrics <- list()
  best_iters <- integer(max(folds))
  for (k in seq_len(max(folds))) {
    tr <- which(folds != k); va <- which(folds == k)
    dtr <- xgb.DMatrix(X[tr, , drop = FALSE], label = y[tr], missing = NA)
    dva <- xgb.DMatrix(X[va, , drop = FALSE], label = y[va], missing = NA)
    n_pos <- sum(y[tr] == 1); n_neg <- sum(y[tr] == 0)
    params_k <- c(params, list(objective = "binary:logistic", eval_metric = "auc",
                               scale_pos_weight = n_neg / max(n_pos, 1), nthread = NTHREAD))
    fit <- xgb.train(params = params_k, data = dtr, nrounds = nrounds,
                     watchlist = list(val = dva), early_stopping_rounds = early_stop, verbose = 0)
    pr <- predict(fit, dva, iterationrange = c(1, fit$best_iteration))
    oof_pred[va] <- pr; best_iters[k] <- fit$best_iteration
    fold_metrics[[k]] <- data.table(
      fold = k, auc = as.numeric(pROC::auc(y[va], pr, quiet = TRUE)),
      logloss = -mean(y[va]*log(pmax(pr, 1e-15)) + (1-y[va])*log(pmax(1-pr, 1e-15))),
      brier = mean((pr - y[va])^2), best_iter = fit$best_iteration)
  }
  list(oof = oof_pred, metrics = rbindlist(fold_metrics), mean_best_iter = round(mean(best_iters)))
}

# ---- main loop --------------------------------------------------------------
all_cv_metrics <- list(); all_best_params <- list()
for (mid in models) {
  cat(strrep("-", 60), "\nTraining:", mid, "\n", strrep("-", 60), "\n", sep = "")
  subset <- recipe[model_id == mid, unique(sample_subset)]
  feats  <- recipe[model_id == mid, feature_id]
  dt <- as.data.table(readRDS(file.path(ML_DIR, sprintf("group1_%s.rds", subset))))

  if ("years_baseline_to_taste" %in% feats && !"years_baseline_to_taste" %in% names(dt))
    dt[, years_baseline_to_taste := age - age_baseline]
  dt <- dt[!is.na(get(OUTCOME))]

  missing_feats <- setdiff(feats, names(dt))
  if (length(missing_feats) > 0) {
    warning(sprintf("[%s] %d recipe features not in rds: %s", mid, length(missing_feats),
                    paste(head(missing_feats, 10), collapse = ", ")))
    feats <- intersect(feats, names(dt))
  }
  X <- as.matrix(dt[, feats, with = FALSE]); mode(X) <- "numeric"
  y <- as.integer(dt[[OUTCOME]])
  cat(sprintf("  subset=%s  N=%d  cases=%d (%.2f%%)  features=%d\n",
              subset, length(y), sum(y), 100 * mean(y), length(feats)))
  folds <- make_folds(y, N_FOLDS)

  grid_results <- list()
  for (g in seq_len(nrow(GRID))) {
    params_g <- as.list(GRID[g]); cv <- cv_one(X, y, folds, params_g)
    grid_results[[g]] <- data.table(model_id = mid, grid_idx = g, as.data.table(params_g),
      mean_auc = mean(cv$metrics$auc), sd_auc = sd(cv$metrics$auc),
      mean_logloss = mean(cv$metrics$logloss), mean_brier = mean(cv$metrics$brier),
      mean_best_iter = cv$mean_best_iter)
    cat(sprintf("  grid %2d/%d  depth=%d eta=%.2f  AUC=%.4f +/- %.4f  best_iter=%d\n",
                g, nrow(GRID), params_g$max_depth, params_g$eta,
                mean(cv$metrics$auc), sd(cv$metrics$auc), cv$mean_best_iter))
  }
  grid_dt <- rbindlist(grid_results); best <- grid_dt[which.max(mean_auc)]
  cat(sprintf("  >>> BEST: depth=%d eta=%.2f  AUC=%.4f  nrounds=%d\n",
              best$max_depth, best$eta, best$mean_auc, best$mean_best_iter))

  best_params <- list(max_depth = best$max_depth, eta = best$eta, subsample = best$subsample,
                      colsample_bytree = best$colsample_bytree, min_child_weight = best$min_child_weight)
  final_cv <- cv_one(X, y, folds, best_params)

  fm <- copy(final_cv$metrics); fm[, model_id := mid]; fm[, sample_subset := subset]
  setcolorder(fm, c("model_id", "sample_subset", "fold")); all_cv_metrics[[mid]] <- fm

  eid_col <- intersect(c("eid", "eid_code", "ID"), names(dt))[1]
  oof_dt <- data.table(eid = if (!is.na(eid_col)) dt[[eid_col]] else NA_integer_,
                       fold = folds, y = y, pred = final_cv$oof)
  saveRDS(oof_dt, file.path(MODEL_DIR, sprintf("%s_oof_preds.rds", mid)))

  n_pos <- sum(y == 1); n_neg <- sum(y == 0)
  dall <- xgb.DMatrix(X, label = y, missing = NA)
  final_params <- c(best_params, list(objective = "binary:logistic", eval_metric = "auc",
                                      scale_pos_weight = n_neg / max(n_pos, 1), nthread = NTHREAD))
  final_fit <- xgb.train(params = final_params, data = dall, nrounds = final_cv$mean_best_iter, verbose = 0)
  saveRDS(list(booster = xgb.save.raw(final_fit), feature_names = feats, params = final_params,
               nrounds = final_cv$mean_best_iter, n_train = length(y), n_cases = sum(y)),
          file.path(MODEL_DIR, sprintf("%s_booster.rds", mid)))
  fwrite(xgb.importance(feature_names = feats, model = final_fit),
         file.path(MODEL_DIR, sprintf("%s_importance.csv", mid)))

  all_best_params[[mid]] <- data.table(model_id = mid, sample_subset = subset,
    n_train = length(y), n_cases = sum(y), n_features = length(feats), as.data.table(best_params),
    nrounds = final_cv$mean_best_iter, cv_mean_auc = best$mean_auc, cv_sd_auc = best$sd_auc,
    cv_mean_brier = best$mean_brier, cv_mean_logloss = best$mean_logloss)
  cat(sprintf("  saved %s_booster.rds / _oof_preds.rds / _importance.csv\n\n", mid))
}

fwrite(rbindlist(all_cv_metrics), file.path(P6_REP_DIR, "xgb_cv_metrics.csv"))
fwrite(rbindlist(all_best_params), file.path(P6_REP_DIR, "xgb_best_params.csv"))
cat("=== 02 complete ===\nNext: 03_evaluate.R\n")
