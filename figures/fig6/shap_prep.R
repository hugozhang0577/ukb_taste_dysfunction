# =============================================================================
# Fig 6 panels E and F — SHAP feature importance
# =============================================================================
#
# TreeSHAP attributions are recomputed here rather than read from disk: the
# training step stores boosters and out-of-fold predictions, not per-sample
# attributions, which are large and cheap to regenerate.
#
# Two models are summarised, the universal-access tier with and without the
# smell features, so panels E and F show what the model uses when the strongest
# predictor is available and what it falls back on when it is not.
#
# Attributions are computed on a stratified sample rather than the full cohort:
# the outcome is rare, so a simple random sample would leave too few cases for a
# stable per-feature mean, and TreeSHAP over the whole cohort buys no extra
# precision for a summary statistic.
#
# Outputs:
#   output/model_reports/eval/shap_{model}_mean_abs.csv   feeds panels E and F
#   output/figures/fig6/shap_cache_{model}.rds            per-sample attributions
# =============================================================================
suppressMessages({ library(data.table); library(xgboost) })
set.seed(42)
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)

ML_DIR  <- "output/ml_ready"
XGB_DIR <- "output/models/xgboost"
OUT      <- "output/figures/fig6"
EVAL_DIR <- "output/model_reports/eval"
for (d in c(OUT, EVAL_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
OUTCOME <- "taste_2w_strict"
N_SHAP  <- 30000     # participants sampled for the attribution computation
N_PLOT  <- 2500      # per-sample points retained in the cache
TOPN    <- 15        # features retained in the cache

MODELS <- list(
  list(id = "M1_TierA",    booster = file.path(XGB_DIR, "M1_TierA_booster.rds")),
  list(id = "M1a_NoSmell", booster = file.path(XGB_DIR, "sensitivity_no_smell", "M1a_NoSmell_booster.rds")))

strat_sample <- function(y, n, seed = 42) {
  if (length(y) <= n) return(seq_along(y))
  set.seed(seed); pos <- which(y == 1); neg <- which(y == 0)
  n_pos <- round(n * length(pos) / length(y))
  c(sample(pos, min(n_pos, length(pos))), sample(neg, min(n - n_pos, length(neg))))
}
rank01 <- function(v) { r <- rank(v, na.last = "keep", ties.method = "average")
  (r - min(r, na.rm = TRUE)) / (max(r, na.rm = TRUE) - min(r, na.rm = TRUE)) }

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
  top <- names(sort(mean_abs, decreasing = TRUE))[1:min(TOPN, length(feats))]
  set.seed(7); pidx <- if (nrow(S) > N_PLOT) sample(nrow(S), N_PLOT) else seq_len(nrow(S))
  shap  <- S[pidx, top, drop = FALSE]
  fval01 <- apply(Xs[pidx, top, drop = FALSE], 2, rank01)
  saveRDS(list(feats = top, mean_abs = mean_abs[top], shap = shap, fval01 = fval01),
          file.path(OUT, sprintf("shap_cache_%s.rds", m$id)))

  # Importance table for the bar panels. Direction is the rank correlation
  # between a feature's value and its own attribution: positive means a higher
  # value pushes the prediction towards case. It is a description of what the
  # fitted model does, not a causal claim, so the wording stays on "risk".
  direction <- vapply(feats, function(f) {
    r <- suppressWarnings(cor(Xs[, f], S[, f], method = "spearman",
                              use = "pairwise.complete.obs"))
    if (is.na(r)) "not monotone"
    else if (r >= 0) "higher value, higher risk" else "higher value, lower risk"
  }, character(1))

  imp <- data.table(feature = feats, mean_abs_shap = as.numeric(mean_abs[feats]),
                    direction = direction)[order(-mean_abs_shap)]
  fwrite(imp, file.path(EVAL_DIR, sprintf("shap_%s_mean_abs.csv", m$id)))
  # Export and dx upload to RAP  (SHAP importance tables feed Fig 6 panels E and F)

  cat(sprintf("  cached %d feats x %d pts | top: %s\n", length(top), length(pidx), top[1]))
}
cat("DONE\n")
