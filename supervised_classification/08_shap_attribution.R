# TreeSHAP attributions for model 1 and its no-smell variant in group1: mean absolute SHAP
# value per feature on a case-stratified sample of 30,000, with the direction of effect.
library(data.table)
library(xgboost)
set.seed(42)
setwd(Sys.getenv("PROJECT_DIR", "."))

XGB_DIR  <- "output/models/xgboost"
EVAL_DIR <- "output/model_reports/eval"
dir.create(EVAL_DIR, recursive = TRUE, showWarnings = FALSE)
OUTCOME <- "taste_2w_strict"
N_SHAP  <- 30000

MODELS <- c(model1          = file.path(XGB_DIR, "model1_booster.rds"),
            model1_no_smell = file.path(XGB_DIR, "sensitivity_no_smell", "model1_no_smell_booster.rds"))

# Stratified by outcome so the sample keeps the case fraction of the full cohort
strat_sample <- function(y, n, seed = 42) {
  if (length(y) <= n) return(seq_along(y))
  set.seed(seed); pos <- which(y == 1); neg <- which(y == 0)
  n_pos <- round(n * length(pos) / length(y))
  c(sample(pos, min(n_pos, length(pos))), sample(neg, min(n - n_pos, length(neg))))
}

dt <- readRDS("output/ml_ready/group1_full.rds")[!is.na(get(OUTCOME))]
for (id in names(MODELS)) {
  obj   <- readRDS(MODELS[[id]])
  bst   <- xgb.load.raw(obj$booster); feats <- obj$feature_names
  idx   <- strat_sample(as.integer(dt[[OUTCOME]]), N_SHAP)
  Xs    <- as.matrix(dt[idx, feats, with = FALSE]); mode(Xs) <- "numeric"
  pc    <- predict(bst, xgb.DMatrix(Xs, missing = NA), predcontrib = TRUE)
  S     <- pc[, -ncol(pc), drop = FALSE]; colnames(S) <- feats   # last column is the bias term

  # Direction = rank correlation between a feature and its own attribution (describes the model)
  direction <- vapply(feats, function(f) {
    r <- suppressWarnings(cor(Xs[, f], S[, f], method = "spearman", use = "pairwise.complete.obs"))
    if (is.na(r)) "not monotone" else if (r >= 0) "higher value, higher risk" else "higher value, lower risk"
  }, character(1))

  imp <- data.table(feature = feats, mean_abs_shap = as.numeric(colMeans(abs(S), na.rm = TRUE)[feats]),
                    direction = direction)[order(-mean_abs_shap)]
  fwrite(imp, file.path(EVAL_DIR, sprintf("shap_%s_mean_abs.csv", id)))
}
