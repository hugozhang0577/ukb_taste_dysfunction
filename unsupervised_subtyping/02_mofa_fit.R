# MOFA+ on group-1 cases of one sex (primary-fit initialization): Gaussian likelihood for the Olink,
# NMR and clinical views, Bernoulli for the prior-diagnosis view; 15 factors.
library(data.table)
library(MOFA2)
source("unsupervised_subtyping/_subtype_map.R")
setwd(Sys.getenv("PROJECT_DIR", "."))

SEX  <- tolower(commandArgs(trailingOnly = TRUE)[1])   # Usage: Rscript unsupervised_subtyping/02_mofa_fit.R M|F
SUF  <- paste0("_", SEX)
SEED <- PRIMARY_SEED[[SEX]]
set.seed(SEED)
for (d in c("output/subtyping/mofa", "output/subtyping/reports")) dir.create(d, showWarnings = FALSE, recursive = TRUE)

data_list <- readRDS(sprintf("output/subtyping/inputs/clustering_input_g1%s.rds", SUF))
mofa <- create_mofa(data_list)
data_opts <- get_default_data_options(mofa)
data_opts$scale_views <- FALSE   # Gaussian views z-scored in 01
model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- 15
model_opts$likelihoods <- c(olink = "gaussian", nmr = "gaussian", clinical_cont = "gaussian",
                            clinical_bin = "gaussian", prior_dx = "bernoulli")[names(data_list)]
train_opts <- get_default_training_options(mofa)
train_opts$convergence_mode <- "slow"
train_opts$maxiter <- 1000
train_opts$seed <- SEED
train_opts$drop_factor_threshold <- 0   # inactive factors are removed in 03
mofa <- prepare_mofa(mofa, data_options = data_opts, model_options = model_opts, training_options = train_opts)
mofa <- run_mofa(mofa, outfile = sprintf("output/subtyping/mofa/mofa_model_g1%s.hdf5", SUF), use_basilisk = TRUE)

Z <- get_factors(mofa, factors = "all", as.data.frame = FALSE)[[1]]   # samples x factors
rownames(Z) <- attr(data_list, "eid")
saveRDS(Z, sprintf("output/subtyping/mofa/mofa_factors_g1%s.rds", SUF))
saveRDS(get_weights(mofa, views = "all", factors = "all", as.data.frame = FALSE),
        sprintf("output/subtyping/mofa/mofa_loadings_g1%s.rds", SUF))

ve_mat <- get_variance_explained(mofa)$r2_per_factor[[1]]   # factor x view, percent
ve_mat[is.nan(ve_mat)] <- 0
fwrite(as.data.table(ve_mat, keep.rownames = "factor"),
       sprintf("output/subtyping/reports/variance_decomposition%s.csv", SUF))
fwrite(data.table(factor = rownames(ve_mat), max_var_any_view = apply(ve_mat, 1, max),
                  total_var = rowSums(ve_mat), n_views_active = rowSums(ve_mat > 1)),
       sprintf("output/subtyping/reports/factor_activity%s.csv", SUF))
