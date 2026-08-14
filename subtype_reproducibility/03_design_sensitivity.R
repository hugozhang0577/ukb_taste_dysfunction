#!/usr/bin/env Rscript
# =============================================================================
# 03_design_sensitivity.R  (design-grid re-fit)
# =============================================================================
#
# Re-fits the subtyping pipeline under a different MOFA design choice and
# re-clusters at the same k, so the resulting partition can be compared with the
# main solution by ARI (04_design_ari.R).
#
# View weighting is the factor this script varies. MOFA's default lets each view
# contribute in proportion to its raw variance, so the two large clinical views
# dominate the fit and the small assay views contribute little; balanced
# weighting equalises per-view variance before fitting. Everything else — views,
# features, likelihoods, factor count, training options, seed, and the
# clustering algorithm and its parameters — is identical to the main run, so any
# difference in the partition is attributable to the weighting alone.
#
# k is fixed rather than re-selected. The point is to compare partitions at the
# same granularity; the natural k of a different design is a separate question
# and is reported in the sweep table this script writes.
#
# Usage: Rscript 03_design_sensitivity.R <M|F> [k]     (k defaults to 4)
# Input:  output/subtyping/inputs/clustering_input_g1_{m,f}.rds   (from 01)
# Output: output/subtyping/mofa/mofa_factors_g1_{m,f}_balanced.rds
#         output/subtyping/clusters/cluster_assignments_g1_{m,f}_balanced_k{k}.rds
#         output/subtyping/reports/design_balanced_{variance,sweep,log}_{m,f}.*
# =============================================================================

rm(list = ls()); gc()
suppressPackageStartupMessages({
  library(data.table)
  library(MOFA2)
  library(cluster)     # silhouette
})

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR); set.seed(20260413)

ARGS <- commandArgs(trailingOnly = TRUE)
SEX  <- if (length(ARGS) >= 1 && nzchar(ARGS[1])) toupper(ARGS[1]) else ""
stopifnot(SEX %in% c("M", "F"))
K    <- if (length(ARGS) >= 2 && nzchar(ARGS[2])) as.integer(ARGS[2]) else 4L
SUF  <- paste0("_", tolower(SEX))
DESIGN <- "balanced"

P7_DIR  <- "output/subtyping"
RPT_DIR <- file.path(P7_DIR, "reports")
for (d in c(file.path(P7_DIR, "mofa"), file.path(P7_DIR, "clusters"), RPT_DIR))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
while (sink.number() > 0) sink()
sink(file.path(RPT_DIR, sprintf("design_%s_log%s.txt", DESIGN, SUF)), split = TRUE)
on.exit({ while (sink.number() > 0) sink() }, add = TRUE)

cat("=== design sensitivity:", DESIGN, "view weighting ===\n")
cat("sex:", SEX, "| k fixed at:", K, "\n")

# ---- [1] the same input the main fit uses, not regenerated -------------------
in_path <- file.path(P7_DIR, "inputs", sprintf("clustering_input_g1%s.rds", SUF))
if (!file.exists(in_path))
  stop("clustering input not found: ", in_path,
       "\n  run unsupervised_subtyping/01_prepare_clustering_features.R first")
data_list <- readRDS(in_path)
n_samples <- unique(vapply(data_list, ncol, integer(1)))
stopifnot(length(n_samples) == 1)
cat(sprintf("[1] views=%d  samples=%d\n", length(data_list), n_samples))

# ---- [2] MOFA, identical to the main fit except the weighting ----------------
mofa <- create_mofa(data_list)

data_opts <- get_default_data_options(mofa)
data_opts$scale_views <- TRUE          # the single difference from the main fit
cat("[2] data_opts$scale_views =", data_opts$scale_views, "(main fit: FALSE)\n")

model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- 15
lik_map <- c(olink = "gaussian", nmr = "gaussian", clinical_cont = "gaussian",
             clinical_bin = "gaussian", phecode = "bernoulli")
model_opts$likelihoods <- lik_map[names(data_list)]

train_opts <- get_default_training_options(mofa)
train_opts$convergence_mode      <- "slow"
train_opts$maxiter               <- 1000
train_opts$seed                  <- 20260413
train_opts$drop_factor_threshold <- 0
train_opts$verbose               <- TRUE

mofa <- prepare_mofa(mofa, data_options = data_opts,
                     model_options = model_opts, training_options = train_opts)

cat("[3] training (30-90 min)\n")
hdf5_path <- file.path(P7_DIR, "mofa", sprintf("mofa_model_g1%s_%s.hdf5", SUF, DESIGN))
if (file.exists(hdf5_path)) file.remove(hdf5_path)
t0 <- Sys.time()
mofa <- run_mofa(mofa, outfile = hdf5_path, use_basilisk = TRUE)
cat(sprintf("    done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

Z <- get_factors(mofa, factors = "all", as.data.frame = FALSE)[[1]]
saveRDS(Z, file.path(P7_DIR, "mofa", sprintf("mofa_factors_g1%s_%s.rds", SUF, DESIGN)))

ve_mat   <- get_variance_explained(mofa)$r2_per_factor[[1]]
activity <- data.table(factor = rownames(ve_mat),
                       total_var = rowSums(ve_mat), n_views_active = rowSums(ve_mat > 1))
fwrite(activity, file.path(RPT_DIR, sprintf("design_%s_variance%s.csv", DESIGN, SUF)))

# ---- [4] cluster, identical parameters to the main clustering step -----------
VAR_THRESH <- 6                                   # active-factor rule, as in the main run
factors_keep <- activity[total_var >= VAR_THRESH, factor]
Z_act <- Z[, factors_keep, drop = FALSE]
stopifnot(all(is.finite(Z_act)))
cat(sprintf("[4] active factors (total_var >= %g%%): %d of %d\n",
            VAR_THRESH, length(factors_keep), ncol(Z)))

D <- dist(Z_act)
sweep_rows <- list()
for (k in 2:8) {
  km <- kmeans(Z_act, centers = k, nstart = 25, iter.max = 100)
  sweep_rows[[length(sweep_rows) + 1]] <- data.table(
    k = k, silhouette = mean(silhouette(km$cluster, D)[, "sil_width"]),
    min_cluster_pct = min(table(km$cluster)) / nrow(Z_act))
}
sweep <- rbindlist(sweep_rows)
fwrite(sweep, file.path(RPT_DIR, sprintf("design_%s_sweep%s.csv", DESIGN, SUF)))
cat("    natural k by silhouette:", sweep[which.max(silhouette), k],
    "| k used for the comparison:", K, "\n")
print(sweep)

km_final <- kmeans(Z_act, centers = K, nstart = 25, iter.max = 100)
final_lab <- km_final$cluster
cat("    cluster sizes:", paste(table(final_lab), collapse = "/"), "\n")

out <- file.path(P7_DIR, "clusters",
                 sprintf("cluster_assignments_g1%s_%s_k%d.rds", SUF, DESIGN, K))
saveRDS(list(design = DESIGN, factors_used = factors_keep, Z_active = Z_act,
             sweep = sweep, final_k = K, final_labels = final_lab,
             silhouette = mean(silhouette(final_lab, D)[, "sil_width"])), out)
# Export and dx upload to RAP  (the design-cell labels 04_design_ari.R compares)
cat("[5] written:", out, "\n")
