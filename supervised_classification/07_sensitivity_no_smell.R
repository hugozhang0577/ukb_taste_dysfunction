#!/usr/bin/env Rscript
# =============================================================================
# 10_sensitivity_no_smell.R  (sensitivity analysis)
#
# Retrain all six tier-models (M1..M6) with every smell-related feature removed
# (variants M1a..M6a), to quantify how much of each tier's AUC depends on
# self-reported smell status and whether omics/EHR/genotype features show
# incremental value once smell dominance is removed.
#
# Per variant: strip smell-* from the tree pool, 5-fold CV grid search (reduced
# grid for small-N tiers), best-params OOF, full-G1 refit + gain importance,
# Platt calibration, held-out G2/G3 validation, paired DeLong vs the original.
#
# Input:  feature_manifest/model_feature_recipes.csv + ml_ready/ + models/xgboost/ (originals)
# Output: models/xgboost/sensitivity_no_smell/* +
#         model_reports/eval/sensitivity_no_smell_{summary,delong_vs_original,grid_search,feature_counts}.csv + plots
# =============================================================================

rm(list = ls()); gc()
while (sink.number() > 0) sink()
suppressMessages({ library(data.table); library(xgboost); library(pROC); library(PRROC) })
set.seed(42)

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)

MF_DIR   <- "output/feature_manifest"
ML_DIR   <- "output/ml_ready"
XGB_DIR  <- "output/models/xgboost"
OUT_DIR  <- file.path(XGB_DIR, "sensitivity_no_smell")
CAL_DIR  <- file.path(OUT_DIR, "calibrators")
PRED_DIR <- file.path(OUT_DIR, "predictions")
P6_REP   <- "output/model_reports"
EVAL_DIR <- file.path(P6_REP, "eval")
PLOT_DIR <- file.path(EVAL_DIR, "plots")
for (d in c(OUT_DIR, CAL_DIR, PRED_DIR, EVAL_DIR, PLOT_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
sink(file.path(P6_REP, "sensitivity_no_smell_log.txt"), split = TRUE)
on.exit({ while (sink.number() > 0) sink() }, add = TRUE)

cat("=== sensitivity - all 6 tiers without smell features ===\n")
OUTCOME <- "taste_2w_strict"; N_FOLDS <- 5; NTHREAD <- 4
MAX_ROUNDS <- 2000; EARLY_STOP <- 50; EPS <- 1e-15; SMALL_N_THR <- 20000L
GRID_FULL  <- CJ(max_depth = c(3, 5, 7), eta = c(0.03, 0.05, 0.1), subsample = 0.8, colsample_bytree = 0.8, min_child_weight = 5)
GRID_SMALL <- CJ(max_depth = c(3, 5), eta = c(0.03, 0.05), subsample = 0.8, colsample_bytree = 0.8, min_child_weight = 5)
VARIANTS <- data.table(
  original   = c("M1_TierA", "M2_TierAB", "M3_TierABC", "M4_TierD_Olink", "M5_TierD_NMR", "M6_TierD_Full"),
  variant_id = c("M1a_NoSmell", "M2a_NoSmell", "M3a_NoSmell", "M4a_NoSmell", "M5a_NoSmell", "M6a_NoSmell"))

# ---- helpers ----------------------------------------------------------------
clip <- function(p) pmax(pmin(p, 1 - EPS), EPS); safe_logit <- function(p) qlogis(clip(p)); safe_log <- function(p) log(clip(p))
make_folds <- function(y, k, seed = 42) {
  set.seed(seed); idx_pos <- which(y == 1); idx_neg <- which(y == 0); f <- integer(length(y))
  f[idx_pos] <- sample(rep_len(seq_len(k), length(idx_pos))); f[idx_neg] <- sample(rep_len(seq_len(k), length(idx_neg))); f
}
cv_one <- function(X, y, folds, params) {
  oof_pred <- rep(NA_real_, length(y)); fold_metrics <- list(); best_iters <- integer(max(folds))
  for (k in seq_len(max(folds))) {
    tr <- which(folds != k); va <- which(folds == k)
    dtr <- xgb.DMatrix(X[tr, , drop = FALSE], label = y[tr], missing = NA)
    dva <- xgb.DMatrix(X[va, , drop = FALSE], label = y[va], missing = NA)
    n_pos <- sum(y[tr] == 1); n_neg <- sum(y[tr] == 0)
    prm <- c(params, list(objective = "binary:logistic", eval_metric = "auc",
                          scale_pos_weight = n_neg / max(n_pos, 1), nthread = NTHREAD))
    fit <- xgb.train(params = prm, data = dtr, nrounds = MAX_ROUNDS, watchlist = list(val = dva),
                     early_stopping_rounds = EARLY_STOP, verbose = 0)
    pr <- predict(fit, dva, iterationrange = c(1, fit$best_iteration))
    oof_pred[va] <- pr; best_iters[k] <- fit$best_iteration
    fold_metrics[[k]] <- data.table(fold = k, auc = as.numeric(pROC::auc(y[va], pr, quiet = TRUE)),
      logloss = -mean(y[va]*safe_log(pr) + (1-y[va])*safe_log(1-pr)), brier = mean((pr - y[va])^2),
      best_iter = fit$best_iteration)
  }
  list(oof = oof_pred, metrics = rbindlist(fold_metrics), mean_best_iter = max(1L, round(mean(best_iters))))
}
build_X <- function(dt, feat_names) { for (m in setdiff(feat_names, names(dt))) dt[, (m) := NA_real_]
  X <- as.matrix(dt[, feat_names, with = FALSE]); mode(X) <- "numeric"; X }
eval_metrics <- function(y, p) {
  if (length(unique(y)) < 2 || sum(y) < 10) return(NULL)
  roc_obj <- pROC::roc(y, p, quiet = TRUE, direction = "<")
  auc_ci <- tryCatch(as.numeric(pROC::ci.auc(roc_obj, method = "delong")), error = function(e) rep(NA_real_, 3))
  pr_obj <- tryCatch(PRROC::pr.curve(scores.class0 = p[y == 1], scores.class1 = p[y == 0], curve = FALSE),
                     error = function(e) list(auc.integral = NA_real_))
  lp <- safe_logit(p)
  list(roc = roc_obj, row = data.table(auc = auc_ci[2], auc_lo = auc_ci[1], auc_hi = auc_ci[3],
    pr_auc = pr_obj$auc.integral, brier = mean((p - y)^2),
    logloss = -mean(y*safe_log(p) + (1-y)*safe_log(1-p)),
    cal_intercept = tryCatch(unname(coef(glm(y ~ offset(lp), family = binomial()))[1]), error = function(e) NA_real_),
    cal_slope = tryCatch(unname(coef(glm(y ~ lp, family = binomial()))[2]), error = function(e) NA_real_)))
}

recipe <- fread(file.path(MF_DIR, "model_feature_recipes.csv"))
best_params_orig <- fread(file.path(P6_REP, "xgb_best_params.csv"))
summary_rows <- list(); grid_rows <- list(); delong_rows <- list(); feat_count_rows <- list()

# ---- main loop --------------------------------------------------------------
for (v_i in seq_len(nrow(VARIANTS))) {
  orig_id <- VARIANTS$original[v_i]; var_id <- VARIANTS$variant_id[v_i]
  cat(strrep("=", 70), "\n[", v_i, "/", nrow(VARIANTS), "] ", var_id, " (original ", orig_id, ")\n", sep = "")
  feats_full <- recipe[model_id == orig_id & pool == "tree", feature_id]
  smell_feats <- grep("^smell", feats_full, ignore.case = TRUE, value = TRUE)
  feats <- setdiff(feats_full, smell_feats)
  cat(sprintf("  features %d -> %d (removed %d smell)\n", length(feats_full), length(feats), length(smell_feats)))
  feat_count_rows[[var_id]] <- data.table(original_id = orig_id, variant_id = var_id,
    n_feat_original = length(feats_full), n_feat_after_smell_strip = length(feats),
    n_smell_removed = length(smell_feats), smell_feats_removed = paste(smell_feats, collapse = ";"))

  subset <- best_params_orig[model_id == orig_id, sample_subset]
  rds_g1 <- file.path(ML_DIR, sprintf("group1_%s.rds", subset))
  if (!file.exists(rds_g1)) { cat("  MISSING", rds_g1, "- skip\n"); next }
  dt <- as.data.table(readRDS(rds_g1))
  if ("years_baseline_to_taste" %in% feats && !"years_baseline_to_taste" %in% names(dt))
    dt[, years_baseline_to_taste := age - age_baseline]
  dt <- dt[!is.na(get(OUTCOME))]
  X <- build_X(dt, feats); y <- as.integer(dt[[OUTCOME]]); N <- length(y); n_pos <- sum(y); n_neg <- N - n_pos
  cat(sprintf("  G1: N=%d cases=%d subset=%s\n", N, n_pos, subset))
  grid <- if (N < SMALL_N_THR) GRID_SMALL else GRID_FULL
  folds <- make_folds(y, N_FOLDS)

  grid_res <- list()
  for (g in seq_len(nrow(grid))) {
    p_g <- as.list(grid[g]); cv <- cv_one(X, y, folds, p_g)
    grid_res[[g]] <- data.table(variant_id = var_id, grid_idx = g, as.data.table(p_g),
      mean_auc = mean(cv$metrics$auc), sd_auc = sd(cv$metrics$auc), mean_best_iter = cv$mean_best_iter)
  }
  grid_dt <- rbindlist(grid_res); grid_rows[[var_id]] <- grid_dt; best <- grid_dt[which.max(mean_auc)]
  cat(sprintf("  >>> BEST depth=%d eta=%.2f AUC=%.4f\n", best$max_depth, best$eta, best$mean_auc))

  best_params <- list(max_depth = best$max_depth, eta = best$eta, subsample = best$subsample,
                      colsample_bytree = best$colsample_bytree, min_child_weight = best$min_child_weight)
  final_cv <- cv_one(X, y, folds, best_params)
  eid_col <- intersect(c("eid", "eid_code", "ID"), names(dt))[1]
  oof_dt <- data.table(eid = if (!is.na(eid_col)) dt[[eid_col]] else NA_integer_, fold = folds, y = y, pred_raw = final_cv$oof)

  dall <- xgb.DMatrix(X, label = y, missing = NA)
  final_params <- c(best_params, list(objective = "binary:logistic", eval_metric = "auc",
                                      scale_pos_weight = n_neg / max(n_pos, 1), nthread = NTHREAD))
  final_fit <- xgb.train(params = final_params, data = dall, nrounds = final_cv$mean_best_iter, verbose = 0)
  saveRDS(list(booster = xgb.save.raw(final_fit), feature_names = feats, params = final_params,
               nrounds = final_cv$mean_best_iter, n_train = N, n_cases = n_pos, variant_id = var_id, original_id = orig_id),
          file.path(OUT_DIR, sprintf("%s_booster.rds", var_id)))
  fwrite(xgb.importance(feature_names = feats, model = final_fit), file.path(OUT_DIR, sprintf("%s_importance.csv", var_id)))
  saveRDS(oof_dt, file.path(OUT_DIR, sprintf("%s_oof_preds.rds", var_id)))

  platt_fit <- tryCatch(glm(oof_dt$y ~ lp_oof, family = binomial(),
                            data = data.frame(y = oof_dt$y, lp_oof = safe_logit(oof_dt$pred_raw))),
                        error = function(e) NULL)
  apply_platt <- function(fit, p) if (is.null(fit)) p else
    as.numeric(predict(fit, newdata = data.frame(lp_oof = safe_logit(p)), type = "response"))
  oof_dt[, pred_platt := apply_platt(platt_fit, pred_raw)]
  saveRDS(oof_dt, file.path(OUT_DIR, sprintf("%s_oof_calibrated.rds", var_id)))
  if (!is.null(platt_fit)) saveRDS(list(method = "platt", model = platt_fit, fit_n = N, fit_cases = n_pos, model_id = var_id),
                                   file.path(CAL_DIR, sprintf("%s_platt.rds", var_id)))

  ev_g1 <- eval_metrics(oof_dt$y, oof_dt$pred_platt)
  if (!is.null(ev_g1)) summary_rows[[paste(var_id, "G1", sep = "_")]] <- cbind(
    data.table(variant_id = var_id, original_id = orig_id, cohort = "G1", n = N, n_cases = n_pos), ev_g1$row)

  for (gnum in 2:3) {
    ch <- paste0("G", gnum); rds <- file.path(ML_DIR, sprintf("group%d_%s.rds", gnum, subset))
    if (!file.exists(rds)) next
    d2 <- as.data.table(readRDS(rds))
    if ("years_baseline_to_taste" %in% feats && !"years_baseline_to_taste" %in% names(d2))
      d2[, years_baseline_to_taste := age - age_baseline]
    d2 <- d2[!is.na(get(OUTCOME))]; if (nrow(d2) < 100) next
    X2 <- build_X(d2, feats); y2 <- as.integer(d2[[OUTCOME]])
    raw2 <- predict(final_fit, xgb.DMatrix(X2, missing = NA)); platt2 <- apply_platt(platt_fit, raw2)
    saveRDS(data.table(eid = if ("eid" %in% names(d2)) d2$eid else NA_integer_, y = y2, pred_raw = raw2, pred_platt = platt2),
            file.path(PRED_DIR, sprintf("%s_%s_preds.rds", var_id, tolower(ch))))
    ev <- eval_metrics(y2, platt2); if (is.null(ev)) next
    summary_rows[[paste(var_id, ch, sep = "_")]] <- cbind(
      data.table(variant_id = var_id, original_id = orig_id, cohort = ch, n = length(y2), n_cases = sum(y2)), ev$row)
  }

  for (ch in c("G1", "G2", "G3")) {
    f_var <- if (ch == "G1") file.path(OUT_DIR, sprintf("%s_oof_calibrated.rds", var_id))
             else file.path(PRED_DIR, sprintf("%s_%s_preds.rds", var_id, tolower(ch)))
    f_orig <- if (ch == "G1") file.path(XGB_DIR, sprintf("%s_oof_calibrated.rds", orig_id))
              else file.path(XGB_DIR, "predictions", sprintf("%s_%s_preds.rds", orig_id, tolower(ch)))
    if (!file.exists(f_var) || !file.exists(f_orig)) next
    a <- as.data.table(readRDS(f_var)); b <- as.data.table(readRDS(f_orig))
    if (!("eid" %in% names(a)) || !("eid" %in% names(b))) next
    m <- merge(a[, .(eid, y_a = y, p_a = pred_platt)], b[, .(eid, y_b = y, p_b = pred_platt)], by = "eid")
    if (nrow(m) < 50 || sum(m$y_a) < 10) next
    ra <- pROC::roc(m$y_a, m$p_a, quiet = TRUE, direction = "<")
    rb <- pROC::roc(m$y_b, m$p_b, quiet = TRUE, direction = "<")
    tt <- tryCatch(pROC::roc.test(ra, rb, method = "delong", paired = TRUE), error = function(e) NULL)
    if (is.null(tt)) next
    delong_rows[[paste(var_id, ch, sep = "_")]] <- data.table(variant_id = var_id, original_id = orig_id, cohort = ch,
      AUC_no_smell = as.numeric(pROC::auc(ra)), AUC_original = as.numeric(pROC::auc(rb)),
      delta_AUC = as.numeric(pROC::auc(ra)) - as.numeric(pROC::auc(rb)),
      z = unname(tt$statistic), p_value = tt$p.value, n_common = nrow(m), paired = TRUE)
  }
  gc()
}

summary_dt <- rbindlist(summary_rows, fill = TRUE); delong_dt <- rbindlist(delong_rows, fill = TRUE)
fwrite(summary_dt, file.path(EVAL_DIR, "sensitivity_no_smell_summary.csv"))
fwrite(delong_dt, file.path(EVAL_DIR, "sensitivity_no_smell_delong_vs_original.csv"))
fwrite(rbindlist(grid_rows, fill = TRUE), file.path(EVAL_DIR, "sensitivity_no_smell_grid_search.csv"))
fwrite(rbindlist(feat_count_rows, fill = TRUE), file.path(EVAL_DIR, "sensitivity_no_smell_feature_counts.csv"))

# ---- delta-AUC forest plot --------------------------------------------------
if (nrow(delong_dt) > 0) {
  png(file.path(PLOT_DIR, "sensitivity_no_smell_delta_auc.png"), 1300, 900, res = 150)
  par(mar = c(4, 14, 3, 2)); dd <- copy(delong_dt)
  dd[, label := paste0(variant_id, " / ", cohort)]
  dd[, cohort_col := fifelse(cohort == "G1", "#1f77b4", fifelse(cohort == "G2", "#ff7f0e", "#2ca02c"))]
  setorder(dd, variant_id, cohort); y_pos <- seq_len(nrow(dd))
  plot(dd$delta_AUC, y_pos, pch = 19, col = dd$cohort_col, xlim = range(c(dd$delta_AUC, 0)) + c(-0.02, 0.02),
       ylim = c(0.5, nrow(dd) + 0.5), xlab = "delta-AUC (no-smell - original); negative = smell was informative",
       ylab = "", yaxt = "n", main = "Sensitivity: AUC change from removing smell features")
  axis(2, at = y_pos, labels = dd$label, las = 1, cex.axis = 0.7); abline(v = 0, lty = 2, col = "grey40")
  legend("bottomleft", c("G1", "G2", "G3"), col = c("#1f77b4", "#ff7f0e", "#2ca02c"), pch = 19, bty = "n")
  dev.off()
}
cat("\n=== 07 complete ===\n")
print(summary_dt[, .(variant_id, cohort, n, n_cases, AUC = round(auc, 4), CalSlope = round(cal_slope, 3))])
