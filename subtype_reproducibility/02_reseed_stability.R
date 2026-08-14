#!/usr/bin/env Rscript
# =============================================================================
# 02_reseed_stability.R
#
# Reseed reproducibility check for the k=4 subtype solution. Re-fits MOFA on the
# SAME full feature set with identical settings, changing ONLY the random seed,
# then re-clusters (k-means, k=4) and computes the Adjusted Rand Index against
# the main solution. Run several seeds to gauge the distribution of ARI.
#
# Everything is identical to the main fit (unsupervised_subtyping/02, /03)
# except train_opts$seed.
#
# Usage: Rscript 02_reseed_stability.R <M|F> <SEED>
# Output (appended): subtyping/reports/reseed_stability_summary.csv
# =============================================================================

rm(list = ls()); gc()
suppressPackageStartupMessages({ library(data.table); library(MOFA2); library(cluster) })

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)

ARGS <- commandArgs(trailingOnly = TRUE)
SEX  <- if (length(ARGS) >= 1 && nzchar(ARGS[1])) toupper(ARGS[1]) else "M"
SEED <- if (length(ARGS) >= 2 && nzchar(ARGS[2])) as.integer(ARGS[2]) else stop("provide a SEED, e.g. 99999")
stopifnot(SEX %in% c("M", "F"))
SEX_SUF <- paste0("_", tolower(SEX))
SUF     <- sprintf("%s_reseed%d", SEX_SUF, SEED)

P7_DIR  <- "output/subtyping"
RPT_DIR <- file.path(P7_DIR, "reports")
VAR_THRESH <- 6   # active-factor rule, same as the main clustering step
K <- 4

adj_rand_index <- function(a, b) {
  tab <- table(a, b); n <- sum(tab); comb2 <- function(x) sum(x * (x - 1) / 2)
  sij <- comb2(as.vector(tab)); sa <- comb2(rowSums(tab)); sb <- comb2(colSums(tab))
  expected <- sa * sb / comb2(n); maxidx <- (sa + sb) / 2
  if ((maxidx - expected) == 0) return(1)
  (sij - expected) / (maxidx - expected)
}

cat(sprintf("=== Reseed reproducibility: SEX=%s SEED=%d ===\n", SEX, SEED))

# ---- [1] load the same clustering input as the main fit ------------------------------------
in_path <- file.path(P7_DIR, "inputs", sprintf("clustering_input_g1%s.rds", SEX_SUF))
if (!file.exists(in_path)) stop("input not found: ", in_path)
data_list <- readRDS(in_path); eid_vec <- attr(data_list, "eid")
cat(sprintf("[1] views=%d, samples=%d\n", length(data_list), unique(sapply(data_list, ncol))))

# ---- [2] MOFA - identical to the main fit, only the seed changes ---------------------
set.seed(SEED)
mofa <- create_mofa(data_list)
data_opts <- get_default_data_options(mofa); data_opts$scale_views <- FALSE
model_opts <- get_default_model_options(mofa); model_opts$num_factors <- 15
model_opts$likelihoods <- c(olink = "gaussian", nmr = "gaussian", clinical_cont = "gaussian",
                            clinical_bin = "gaussian", phecode = "bernoulli")[names(data_list)]
train_opts <- get_default_training_options(mofa)
train_opts$convergence_mode <- "slow"; train_opts$maxiter <- 1000
train_opts$seed <- SEED; train_opts$drop_factor_threshold <- 0; train_opts$verbose <- FALSE
mofa <- prepare_mofa(mofa, data_options = data_opts, model_options = model_opts, training_options = train_opts)
cat(sprintf("[2] training MOFA (seed=%d)\n", SEED))
hdf5 <- file.path(P7_DIR, "mofa", sprintf("mofa_model_g1%s.hdf5", SUF))
if (file.exists(hdf5)) file.remove(hdf5)
mofa <- run_mofa(mofa, outfile = hdf5, use_basilisk = TRUE)

# ---- [3] active factors + k-means k=4 ---------------------------------------
Z <- get_factors(mofa, factors = "all", as.data.frame = FALSE)[[1]]; rownames(Z) <- eid_vec
ve <- get_variance_explained(mofa)$r2_per_factor[[1]]; ve[is.nan(ve)] <- 0
keep <- names(rowSums(ve))[rowSums(ve) >= VAR_THRESH]
cat(sprintf("[3] active factors (total_var>=%g): %d\n", VAR_THRESH, length(keep)))
Z_act <- Z[, keep, drop = FALSE]
set.seed(SEED); km <- kmeans(Z_act, centers = K, nstart = 25, iter.max = 100)
reseed_lab <- km$cluster; names(reseed_lab) <- rownames(Z_act)
sil <- mean(silhouette(km$cluster, dist(Z_act))[, "sil_width"])
cat(sprintf("    silhouette=%.3f  sizes: %s\n", sil, paste(table(reseed_lab), collapse = "/")))
saveRDS(list(final_labels = reseed_lab, factors_used = keep, silhouette = sil),
        file.path(P7_DIR, "clusters", sprintf("cluster_assignments_g1%s_k%d.rds", SUF, K)))

# ---- [4] ARI vs main solution -----------------------------------------------
main_path <- file.path(P7_DIR, "clusters", sprintf("cluster_assignments_g1%s_k4.rds", SEX_SUF))
if (!file.exists(main_path)) stop("main labels not found: ", main_path)
main_lab <- readRDS(main_path)$final_labels
common <- intersect(names(main_lab), names(reseed_lab))
ari <- adj_rand_index(as.integer(main_lab[common]), as.integer(reseed_lab[common]))
cat(sprintf("\n[4] ARI(main vs reseed) = %.4f (n_common=%d)\n", ari, length(common)))

# ---- [5] append to summary --------------------------------------------------
row <- data.table(sex = SEX, seed = SEED, n_factors = length(keep), silhouette = round(sil, 3),
                  ARI_vs_main = round(ari, 4), cluster_sizes = paste(table(reseed_lab), collapse = "/"))
sum_f <- file.path(RPT_DIR, "reseed_stability_summary.csv")
fwrite(row, sum_f, append = file.exists(sum_f))
cat("appended to:", sum_f, "\n")
