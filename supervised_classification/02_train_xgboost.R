# XGBoost for each nested model in group1: stratified 5-fold CV over a 3 x 3 grid of tree
# depth and learning rate, boosting rounds set by early stopping on an inner split of each
# training fold; out-of-fold predictions are kept and a final model is fitted on all of group1.
library(data.table)
library(xgboost)
library(pROC)
set.seed(42)
setwd(Sys.getenv("PROJECT_DIR", "."))

ML_DIR    <- "output/ml_ready"
MODEL_DIR <- "output/models/xgboost"
REP_DIR   <- "output/model_reports"
dir.create(MODEL_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(REP_DIR,   showWarnings = FALSE, recursive = TRUE)

OUTCOME    <- "taste_2w_strict"
N_FOLDS    <- 5
NTHREAD    <- 4
MAX_ROUNDS <- 2000
EARLY_STOP <- 50
GRID <- CJ(max_depth = c(3, 5, 7), eta = c(0.03, 0.05, 0.1),
           subsample = 0.8, colsample_bytree = 0.8, min_child_weight = 5)

recipe <- fread("output/feature_manifest/model_feature_recipes.csv")

make_folds <- function(y, k, seed = 42) {
  set.seed(seed)
  fold <- integer(length(y))
  fold[y == 1] <- sample(rep_len(seq_len(k), sum(y == 1)))
  fold[y == 0] <- sample(rep_len(seq_len(k), sum(y == 0)))
  fold
}

# Within each training fold, a stratified fifth is held back for early stopping and grid
# choice; the outer fold is only scored. Classes are weighted by scale_pos_weight.
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
    prm <- c(params, list(objective = "binary:logistic", eval_metric = "auc",
                          eval_metric = "logloss",   # the last metric drives early stopping
                          scale_pos_weight = spw, nthread = NTHREAD))
    fit <- xgb.train(params = prm, data = dtr, nrounds = MAX_ROUNDS, watchlist = list(es = des),
                     early_stopping_rounds = EARLY_STOP, verbose = 0)
    pr <- predict(fit, dva, iterationrange = c(1, fit$best_iteration))
    oof[va] <- pr
    metrics[[k]] <- data.table(
      auc = as.numeric(pROC::auc(y[va], pr, quiet = TRUE)),
      logloss = -mean(y[va] * log(pmax(pr, 1e-15)) + (1 - y[va]) * log(pmax(1 - pr, 1e-15))),
      brier = mean((pr - y[va])^2), best_iter = fit$best_iteration,
      es_auc = fit$evaluation_log[iter == fit$best_iteration, es_auc])
  }
  metrics <- rbindlist(metrics)
  list(oof = oof, metrics = metrics, mean_best_iter = round(mean(metrics$best_iter)))
}

best_rows <- list()
for (mid in unique(recipe$model_id)) {
  subset <- recipe[model_id == mid, unique(sample_subset)]
  feats  <- recipe[model_id == mid, feature_id]
  dt <- readRDS(file.path(ML_DIR, sprintf("group1_%s.rds", subset)))
  dt <- dt[!is.na(get(OUTCOME))]
  X <- as.matrix(dt[, feats, with = FALSE]); mode(X) <- "numeric"
  y <- as.integer(dt[[OUTCOME]])
  folds <- make_folds(y, N_FOLDS)

  grid_dt <- rbindlist(lapply(seq_len(nrow(GRID)), function(g) {
    cv <- cv_one(X, y, folds, as.list(GRID[g]))
    data.table(GRID[g], mean_auc = mean(cv$metrics$auc), sd_auc = sd(cv$metrics$auc),
               mean_logloss = mean(cv$metrics$logloss), mean_brier = mean(cv$metrics$brier),
               mean_best_iter = cv$mean_best_iter, mean_es_auc = mean(cv$metrics$es_auc))
  }))
  best <- grid_dt[which.max(mean_es_auc)]   # chosen on the inner split, not the outer folds
  best_params <- as.list(best[, .(max_depth, eta, subsample, colsample_bytree, min_child_weight)])
  final_cv <- cv_one(X, y, folds, best_params)
  saveRDS(data.table(eid = dt$eid, fold = folds, y = y, pred = final_cv$oof),
          file.path(MODEL_DIR, sprintf("%s_oof_preds.rds", mid)))

  final_params <- c(best_params, list(objective = "binary:logistic", eval_metric = "auc",
                                      scale_pos_weight = sum(y == 0) / max(sum(y == 1), 1),
                                      nthread = NTHREAD))
  final_fit <- xgb.train(params = final_params, data = xgb.DMatrix(X, label = y, missing = NA),
                         nrounds = final_cv$mean_best_iter, verbose = 0)
  saveRDS(list(booster = xgb.save.raw(final_fit), feature_names = feats, params = final_params,
               nrounds = final_cv$mean_best_iter, n_train = length(y), n_cases = sum(y)),
          file.path(MODEL_DIR, sprintf("%s_booster.rds", mid)))
  fwrite(xgb.importance(feature_names = feats, model = final_fit),
         file.path(MODEL_DIR, sprintf("%s_importance.csv", mid)))

  best_rows[[mid]] <- data.table(model_id = mid, sample_subset = subset, n_train = length(y),
    n_cases = sum(y), n_features = length(feats), as.data.table(best_params),
    nrounds = final_cv$mean_best_iter, cv_mean_auc = best$mean_auc, cv_sd_auc = best$sd_auc,
    cv_mean_brier = best$mean_brier, cv_mean_logloss = best$mean_logloss)
}
fwrite(rbindlist(best_rows), file.path(REP_DIR, "xgb_best_params.csv"))
