#!/usr/bin/env Rscript
# =============================================================================
# 02_mofa_fit.R  (subtyping core)
#
# Fit MOFA+ on the G1 case-only multi-view data from 01_prepare_clustering_features. Produces factor scores
# (samples x K), per-view loadings, and a variance decomposition for downstream
# clustering (03) and the discriminator step (04).
#
# Model:
#   views       : olink, nmr, clinical_cont, clinical_bin, phecode
#   likelihoods : gaussian (first 4), bernoulli (phecode)
#   num_factors : 15 (pruned post-hoc by variance threshold in 03_cluster_on_factors)
#   training    : slow convergence; seed = 20260413; native NA handling
#
# Optional CLI arg "M"/"F" for the sex-stratified fit.
#
# Input:  subtyping/inputs/clustering_input_g1{,_m,_f}.rds
# Output: subtyping/mofa/mofa_{model.hdf5,object,factors,loadings}_g1*  +  subtyping/reports/factor_activity + variance_decomposition + plots
# =============================================================================

rm(list = ls()); gc()
suppressPackageStartupMessages({ library(data.table); library(MOFA2); library(ggplot2) })

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR); set.seed(20260413)

ARGS <- commandArgs(trailingOnly = TRUE)
SEX  <- if (length(ARGS) >= 1 && nzchar(ARGS[1])) toupper(ARGS[1]) else ""
stopifnot(SEX %in% c("", "M", "F"))
SUF  <- if (nzchar(SEX)) paste0("_", tolower(SEX)) else ""

P7_DIR  <- "output/subtyping"
RPT_DIR <- file.path(P7_DIR, "reports")
FIG_DIR <- file.path(RPT_DIR, "fig")
dir.create(file.path(P7_DIR, "mofa"), showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
while (sink.number() > 0) sink()
sink(file.path(RPT_DIR, sprintf("mofa_training_log%s.txt", SUF)), split = TRUE)
on.exit(sink(), add = TRUE)

cat("=== MOFA+ fit (G1 cases) ===\n")
cat(sprintf("sex stratification: %s\n", if (nzchar(SEX)) SEX else "none"))

# ---- [1] load input ---------------------------------------------------------
data_list <- readRDS(file.path(P7_DIR, "inputs", sprintf("clustering_input_g1%s.rds", SUF)))
n_samples <- unique(sapply(data_list, ncol)); stopifnot(length(n_samples) == 1)
cat(sprintf("  views=%d, samples=%d\n", length(data_list), n_samples))
for (v in names(data_list))
  cat(sprintf("    %-14s : %4d features, %.1f%% NA\n", v, nrow(data_list[[v]]), 100 * mean(is.na(data_list[[v]]))))

# ---- [2] construct MOFA object ----------------------------------------------
cat("[2] building MOFA object\n")
mofa <- create_mofa(data_list)
data_opts <- get_default_data_options(mofa)
data_opts$scale_views <- FALSE  # Gaussian views already z-scored in 01

model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- 15
# clinical_bin uses the gaussian likelihood: with ~24.5% missingness the Bernoulli
# ELBO is numerically unstable; gaussian on 0/1 data is the recommended workaround
# and retains factor interpretation. phecode remains bernoulli (low NA, stable).
lik_map <- c(olink = "gaussian", nmr = "gaussian", clinical_cont = "gaussian",
             clinical_bin = "gaussian", phecode = "bernoulli")
model_opts$likelihoods <- lik_map[names(data_list)]
cat("  likelihoods:\n"); print(model_opts$likelihoods)

train_opts <- get_default_training_options(mofa)
train_opts$convergence_mode <- "slow"
train_opts$maxiter <- 1000
train_opts$seed <- 20260413
train_opts$drop_factor_threshold <- 0   # prune post-hoc in 03_cluster_on_factors
train_opts$verbose <- TRUE
mofa <- prepare_mofa(mofa, data_options = data_opts, model_options = model_opts, training_options = train_opts)

# ---- [3] train --------------------------------------------------------------
cat("\n[3] training MOFA+ (30-90 min)\n")
hdf5_path <- file.path(P7_DIR, "mofa", sprintf("mofa_model_g1%s.hdf5", SUF))
if (file.exists(hdf5_path)) file.remove(hdf5_path)
t0 <- Sys.time()
mofa <- run_mofa(mofa, outfile = hdf5_path, use_basilisk = TRUE)
cat(sprintf("  training wall time: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
saveRDS(mofa, file.path(P7_DIR, "mofa", sprintf("mofa_object_g1%s.rds", SUF)))

# ---- [4] extract factors / loadings / variance ------------------------------
cat("\n[4] extracting outputs\n")
Z <- get_factors(mofa, factors = "all", as.data.frame = FALSE)[[1]]   # samples x K
rownames(Z) <- attr(data_list, "eid")
cat(sprintf("  factor matrix: %d samples x %d factors\n", nrow(Z), ncol(Z)))
saveRDS(Z, file.path(P7_DIR, "mofa", sprintf("mofa_factors_g1%s.rds", SUF)))
saveRDS(get_weights(mofa, views = "all", factors = "all", as.data.frame = FALSE),
        file.path(P7_DIR, "mofa", sprintf("mofa_loadings_g1%s.rds", SUF)))

ve_mat <- get_variance_explained(mofa)$r2_per_factor[[1]]  # factor x view, percent
# NaN entries (Bernoulli views with near-zero predictor-scale variance) -> 0.
# Does not affect clustering, which uses Z not ve_mat.
n_nan <- sum(is.nan(ve_mat)); if (n_nan > 0) { cat(sprintf("  %d NaN -> 0\n", n_nan)); ve_mat[is.nan(ve_mat)] <- 0 }
fwrite(as.data.table(ve_mat, keep.rownames = "factor"),
       file.path(RPT_DIR, sprintf("variance_decomposition%s.csv", SUF)))
cat("  variance decomposition (factor x view, %):\n"); print(round(ve_mat, 2))
if (max(ve_mat) == 0) warning("all variance explained == 0; check factor SDs")

activity <- data.table(factor = rownames(ve_mat), max_var_any_view = apply(ve_mat, 1, max),
                       total_var = rowSums(ve_mat), n_views_active = rowSums(ve_mat > 1))
fwrite(activity, file.path(RPT_DIR, sprintf("factor_activity%s.csv", SUF)))

# ---- [5] figures ------------------------------------------------------------
cat("\n[5] figures\n")
png(file.path(FIG_DIR, sprintf("variance_heatmap%s.png", SUF)), width = 1400, height = 1000, res = 180)
print(plot_variance_explained(mofa, x = "view", y = "factor", max_r2 = max(ve_mat, na.rm = TRUE)))
dev.off()
p_bar <- ggplot(activity, aes(x = reorder(factor, -total_var), y = total_var)) +
  geom_col(fill = "steelblue") +
  labs(x = "Factor", y = "Total variance explained (% sum across views)", title = "MOFA+ factor activity (G1 cases)") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(FIG_DIR, sprintf("factor_variance_bar%s.png", SUF)), p_bar, width = 8, height = 5, dpi = 180)
cat("\n=== done ===\n")
