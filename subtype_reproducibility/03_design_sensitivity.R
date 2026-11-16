# Design sensitivity (variance-balanced view weighting): MOFA+ refitted with scale_views = TRUE, otherwise
# identical to the primary fit, then k-means on the active factors (k = 2-8 sweep; k = 4 for comparison).
library(data.table)
library(MOFA2)
library(cluster)
source("unsupervised_subtyping/_subtype_map.R")
setwd(Sys.getenv("PROJECT_DIR", "."))

SEX  <- tolower(commandArgs(trailingOnly = TRUE)[1])   # Usage: Rscript subtype_reproducibility/03_design_sensitivity.R M|F
K    <- 4L
SEED <- PRIMARY_SEED[[SEX]]
set.seed(SEED)

data_list <- readRDS(sprintf("output/subtyping/inputs/clustering_input_g1_%s.rds", SEX))
mofa <- create_mofa(data_list)
data_opts <- get_default_data_options(mofa)
data_opts$scale_views <- TRUE   # the only difference from the primary fit
model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- 15
model_opts$likelihoods <- c(olink = "gaussian", nmr = "gaussian", clinical_cont = "gaussian",
                            clinical_bin = "gaussian", prior_dx = "bernoulli")[names(data_list)]
train_opts <- get_default_training_options(mofa)
train_opts$convergence_mode <- "slow"
train_opts$maxiter <- 1000
train_opts$seed <- SEED
train_opts$drop_factor_threshold <- 0
mofa <- prepare_mofa(mofa, data_options = data_opts, model_options = model_opts, training_options = train_opts)
mofa <- run_mofa(mofa, outfile = sprintf("output/subtyping/mofa/mofa_model_g1_%s_balanced.hdf5", SEX), use_basilisk = TRUE)

Z <- get_factors(mofa, factors = "all", as.data.frame = FALSE)[[1]]
rownames(Z) <- attr(data_list, "eid")
saveRDS(Z, sprintf("output/subtyping/mofa/mofa_factors_g1_%s_balanced.rds", SEX))
ve_mat <- get_variance_explained(mofa)$r2_per_factor[[1]]
factors_keep <- rownames(ve_mat)[which(rowSums(ve_mat) >= 6)]   # active factors, as in the primary fit
Z_act <- Z[, factors_keep, drop = FALSE]
D <- dist(Z_act)

sweep <- rbindlist(lapply(2:8, function(k) {
  km <- kmeans(Z_act, centers = k, nstart = 25, iter.max = 100)
  data.table(k = k, silhouette = mean(silhouette(km$cluster, D)[, "sil_width"]),
             min_cluster_pct = min(table(km$cluster)) / nrow(Z_act))
}))
fwrite(sweep, sprintf("output/subtyping/reports/design_balanced_sweep_%s.csv", SEX))
final_lab <- kmeans(Z_act, centers = K, nstart = 25, iter.max = 100)$cluster
saveRDS(list(design = "balanced", factors_used = factors_keep, Z_active = Z_act, sweep = sweep, final_k = K,
             final_labels = final_lab, silhouette = mean(silhouette(final_lab, D)[, "sil_width"])),
        sprintf("output/subtyping/clusters/cluster_assignments_g1_%s_balanced_k%d.rds", SEX, K))
