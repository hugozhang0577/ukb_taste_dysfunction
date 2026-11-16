# Seed re-runs: MOFA+ and k-means (k = 4) refitted from each of the 20 non-primary initializations per sex;
# the summary compares all 21 fits by ARI and scores recovery of the psychosomatic subtype (B) in both sexes.
# Usage: Rscript subtype_reproducibility/02_reseed_stability.R M|F <seed>  (each of other_seeds()), then ... summary
library(data.table)
library(MOFA2)
library(cluster)
source("unsupervised_subtyping/_subtype_map.R")
setwd(Sys.getenv("PROJECT_DIR", "."))

ARGS    <- tolower(commandArgs(trailingOnly = TRUE))
CL_DIR  <- "output/subtyping/clusters"
RPT_DIR <- "output/subtyping/reports"
K <- 4
COSINE_THR <- 0.80   # pre-specified recovery criterion for the psychosomatic subtype

adj_rand_index <- function(a, b) {
  tab <- table(a, b); n <- sum(tab); comb2 <- function(x) sum(x * (x - 1) / 2)
  sij <- comb2(as.vector(tab)); sa <- comb2(rowSums(tab)); sb <- comb2(colSums(tab))
  expected <- sa * sb / comb2(n); maxidx <- (sa + sb) / 2
  if ((maxidx - expected) == 0) return(1)
  (sij - expected) / (maxidx - expected)
}

# ==== one re-initialization: identical to the primary fit except the seed
if (ARGS[1] != "summary") {
  SEX <- ARGS[1]; SEED <- as.integer(ARGS[2])
  data_list <- readRDS(sprintf("output/subtyping/inputs/clustering_input_g1_%s.rds", SEX))
  set.seed(SEED)
  mofa <- create_mofa(data_list)
  data_opts <- get_default_data_options(mofa); data_opts$scale_views <- FALSE
  model_opts <- get_default_model_options(mofa); model_opts$num_factors <- 15
  model_opts$likelihoods <- c(olink = "gaussian", nmr = "gaussian", clinical_cont = "gaussian",
                              clinical_bin = "gaussian", prior_dx = "bernoulli")[names(data_list)]
  train_opts <- get_default_training_options(mofa)
  train_opts$convergence_mode <- "slow"; train_opts$maxiter <- 1000
  train_opts$seed <- SEED; train_opts$drop_factor_threshold <- 0
  mofa <- prepare_mofa(mofa, data_options = data_opts, model_options = model_opts, training_options = train_opts)
  mofa <- run_mofa(mofa, outfile = sprintf("output/subtyping/mofa/mofa_model_g1_%s_reseed%d.hdf5", SEX, SEED),
                   use_basilisk = TRUE)

  Z <- get_factors(mofa, factors = "all", as.data.frame = FALSE)[[1]]; rownames(Z) <- attr(data_list, "eid")
  ve <- get_variance_explained(mofa)$r2_per_factor[[1]]; ve[is.nan(ve)] <- 0
  keep <- names(rowSums(ve))[rowSums(ve) >= 6]
  Z_act <- Z[, keep, drop = FALSE]
  set.seed(SEED); km <- kmeans(Z_act, centers = K, nstart = 25, iter.max = 100)
  saveRDS(list(final_labels = setNames(km$cluster, rownames(Z_act)), factors_used = keep,
               silhouette = mean(silhouette(km$cluster, dist(Z_act))[, "sil_width"])),
          file.path(CL_DIR, sprintf("cluster_assignments_g1_%s_reseed%d_k%d.rds", SEX, SEED, K)))
  quit(save = "no")
}

# ==== summary: the 20 re-initializations vs the primary fit, and median ARI of each fit to the rest
lab_of <- function(p) { o <- readRDS(p); setNames(as.integer(o$final_labels), names(o$final_labels)) }
primary_path <- function(sx) file.path(CL_DIR, sprintf("cluster_assignments_g1_%s_k4.rds", sx))
reseed_path  <- function(sx, s) file.path(CL_DIR, sprintf("cluster_assignments_g1_%s_reseed%d_k%d.rds", sx, s, K))
rows <- list(); pw <- list()
for (sx in c("m", "f")) {
  prim <- lab_of(primary_path(sx))
  L <- list(); L[[as.character(PRIMARY_SEED[[sx]])]] <- prim
  for (s in other_seeds(sx)) {
    o <- readRDS(reseed_path(sx, s)); lb <- as.integer(o$final_labels[names(prim)])
    L[[as.character(s)]] <- lb
    rows[[length(rows) + 1]] <- data.table(sex = toupper(sx), seed = s, n_factors = length(o$factors_used),
      silhouette = round(o$silhouette, 3), ARI_vs_main = round(adj_rand_index(prim, lb), 4),
      cluster_sizes = paste(table(o$final_labels), collapse = "/"))
  }
  nm <- names(L); A <- matrix(1, length(nm), length(nm), dimnames = list(nm, nm))
  for (i in seq_along(nm)) for (j in seq_along(nm)) if (i < j) A[i, j] <- A[j, i] <- adj_rand_index(L[[i]], L[[j]])
  mdn <- vapply(nm, function(n) median(A[n, setdiff(nm, n)]), numeric(1))
  mn  <- vapply(nm, function(n) mean(A[n, setdiff(nm, n)]), numeric(1))
  a <- rbindlist(rows)[sex == toupper(sx), ARI_vs_main]; P <- A[upper.tri(A)]
  pw[[sx]] <- data.table(sex = toupper(sx), primary_seed = PRIMARY_SEED[[sx]], n_fits = length(L),
    median_ARI_others_vs_primary = round(median(a), 3), n_others_ge06 = sum(a >= 0.6),
    median_pairwise_ARI = round(median(P), 3),
    iqr_pairwise_ARI = paste(round(quantile(P, c(.25, .75)), 3), collapse = "-"),
    highest_median_ARI_seed = as.integer(nm[order(-mdn, -mn)[1]]), primary_median_ARI = round(mdn[[1]], 3),
    highest_median_ARI = round(max(mdn), 3))
}
fwrite(setorder(rbindlist(rows), sex, seed), file.path(RPT_DIR, "reseed_stability_summary.csv"))
fwrite(rbindlist(pw), file.path(RPT_DIR, "reseed_pairwise_ari.csv"))

# ==== psychosomatic subtype (B) recovery, same rule in both sexes: the cluster closest to the female-primary
# B centroid has cosine >= 0.80 and is that fit's highest-depression cluster
cases <- as.data.table(readRDS("output/ml_ready/group1_full.rds"))[taste_2w_strict == 1]
depv  <- cases$phq2_score
feats <- intersect(unique(c(fread(file.path(RPT_DIR, "fig_tables", "table_cluster_topfeatures_m_k4.csv"))$feature_id,
                            fread(file.path(RPT_DIR, "fig_tables", "table_cluster_topfeatures_f_k4.csv"))$feature_id)),
                   names(cases))
M <- as.matrix(cases[, feats, with = FALSE]); storage.mode(M) <- "numeric"
eid <- as.character(cases$eid); Mz <- M * NA_real_
for (s in c(0, 1)) {   # z-score within sex
  r <- which(cases$sex == s); Ms <- M[r, , drop = FALSE]
  mu <- apply(Ms, 2, mean, na.rm = TRUE); sdv <- apply(Ms, 2, sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv == 0] <- 1; Mz[r, ] <- sweep(sweep(Ms, 2, mu, "-"), 2, sdv, "/")
}
Mz[!is.finite(Mz)] <- NA
cent <- function(lv) t(vapply(1:K, function(cl) {
  r <- match(names(lv)[lv == cl], eid); colMeans(Mz[r[!is.na(r)], , drop = FALSE], na.rm = TRUE)
}, numeric(ncol(Mz))))
clustdep <- function(lv) vapply(1:K, function(cl) {
  r <- match(names(lv)[lv == cl], eid); mean(depv[r[!is.na(r)]], na.rm = TRUE)
}, numeric(1))
cosine <- function(a, b) { ok <- is.finite(a) & is.finite(b); if (sum(ok) < 5) return(NA_real_)
  sum(a[ok] * b[ok]) / (sqrt(sum(a[ok]^2)) * sqrt(sum(b[ok]^2))) }

flab <- lab_of(primary_path("f"))
fBc  <- cent(flab)[which.max(clustdep(flab)), ]   # female B = highest-depression female cluster
brows <- list()
for (sx in c("m", "f")) {
  fits <- c(setNames(primary_path(sx), PRIMARY_SEED[[sx]]),
            setNames(vapply(other_seeds(sx), function(s) reseed_path(sx, s), character(1)), other_seeds(sx)))
  for (s in names(fits)) {
    ml <- lab_of(fits[[s]]); Cm <- cent(ml); dep <- clustdep(ml)
    cs <- vapply(1:K, function(i) cosine(Cm[i, ], fBc), numeric(1))
    best <- which.max(cs); is_topdep <- best == which.max(dep)
    brows[[length(brows) + 1]] <- data.table(sex = toupper(sx), seed = as.integer(s),
      primary = s == as.character(PRIMARY_SEED[[sx]]), best_cluster = best,
      cos_to_female_B = round(cs[best], 3), cluster_depression = round(dep[best], 2),
      is_highest_depression = is_topdep, B_recovered = cs[best] >= COSINE_THR && is_topdep)
  }
}
fwrite(rbindlist(brows), file.path(RPT_DIR, "b_recovery_frequency.csv"))
