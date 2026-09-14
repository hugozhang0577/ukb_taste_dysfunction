#> in   output/ml_ready/group1_full.rds
suppressMessages({ library(data.table); library(xgboost) })
set.seed(42)
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)

ML_DIR   <- "output/ml_ready"
XGB_DIR  <- "output/models/xgboost"
EVAL_DIR <- "output/model_reports/eval"
dir.create(EVAL_DIR, recursive = TRUE, showWarnings = FALSE)
OUTCOME <- "taste_2w_strict"
# Attributions are computed on a stratified sample: the outcome is rare, so a
# simple random sample leaves too few cases for a stable per-feature mean.
N_SHAP  <- 30000

MODELS <- list(
  list(id = "M1_TierA",    booster = file.path(XGB_DIR, "M1_TierA_booster.rds")),
  list(id = "M1a_NoSmell", booster = file.path(XGB_DIR, "sensitivity_no_smell", "M1a_NoSmell_booster.rds")))

strat_sample <- function(y, n, seed = 42) {
  if (length(y) <= n) return(seq_along(y))
  set.seed(seed); pos <- which(y == 1); neg <- which(y == 0)
  n_pos <- round(n * length(pos) / length(y))
  c(sample(pos, min(n_pos, length(pos))), sample(neg, min(n - n_pos, length(neg))))
}

for (m in MODELS) {
  cat("---", m$id, "---\n")
  obj <- readRDS(m$booster); bst <- xgb.load.raw(obj$booster); feats0 <- obj$feature_names
  dt <- readRDS(file.path(ML_DIR, "group1_full.rds")); setDT(dt)
  if ("years_baseline_to_taste" %in% feats0 && !"years_baseline_to_taste" %in% names(dt))
    dt[, years_baseline_to_taste := age - age_baseline]
  dt <- dt[!is.na(get(OUTCOME))]
  feats <- intersect(feats0, names(dt))
  y <- as.integer(dt[[OUTCOME]]); idx <- strat_sample(y, N_SHAP)
  Xs <- as.matrix(dt[idx, feats, with = FALSE]); mode(Xs) <- "numeric"
  pc <- predict(bst, xgb.DMatrix(Xs, missing = NA), predcontrib = TRUE)
  S <- pc[, -ncol(pc), drop = FALSE]; colnames(S) <- feats
  mean_abs <- colMeans(abs(S), na.rm = TRUE)

  # Direction is the rank correlation between a feature and its own
  # attribution: it describes the fitted model, it is not a causal claim.
  direction <- vapply(feats, function(f) {
    r <- suppressWarnings(cor(Xs[, f], S[, f], method = "spearman",
                              use = "pairwise.complete.obs"))
    if (is.na(r)) "not monotone"
    else if (r >= 0) "higher value, higher risk" else "higher value, lower risk"
  }, character(1))

  imp <- data.table(feature = feats, mean_abs_shap = as.numeric(mean_abs[feats]),
                    direction = direction)[order(-mean_abs_shap)]
  fwrite(imp, file.path(EVAL_DIR, sprintf("shap_%s_mean_abs.csv", m$id)))
  # Export and dx upload to RAP  (the SHAP importance tables)

  cat(sprintf("  %d features | top: %s\n", length(feats), imp$feature[1]))
}
cat("DONE\n")
