#!/usr/bin/env Rscript
# =============================================================================
# 04_crosspopulation_projection.R  (held-out cross-population validation)
#
# Validate the G1 subtypes on the held-out G2 + G3 cohorts by projection:
#   1. load G1 MOFA loadings + G1 cluster centroids in factor space
#   2. for each cohort {g2, g3}:
#        a. load the (G1-z-scored) clustering input
#        b. project samples into G1 factor space via the joint-Gaussian-views
#           posterior  Z_new = (sum_v W_v^T W_v + ridge*I)^-1 sum_v W_v^T X_v
#           (NA-aware per sample; phecode/Bernoulli view skipped)
#        c. restrict to G1 active factors
#        d. assign each sample to the nearest G1 centroid -> projected labels
#        e. independently k-means cluster (k = k_G1), Hungarian-align to G1
#        f. concordance: ARI(projected, independent), silhouette, prevalence,
#           mean centroid distance
#
# CLI: arg1 = ""/"M"/"F" (matches the per-sex stratification of 02 and 03).
#
# Input:  subtyping/inputs/clustering_input_g{2,3}*.rds + subtyping/mofa/mofa_loadings_g1*.rds
#         + subtyping/clusters/cluster_assignments_g1*_k4.rds
# Output: subtyping/clusters/cluster_assignments_projected_g{2,3}*.rds + subtyping/reports/crosspop_*
# =============================================================================

rm(list = ls()); gc()
suppressPackageStartupMessages({ library(data.table); library(cluster); library(clue); library(mclust) })

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR); set.seed(20260413)

ARGS <- commandArgs(trailingOnly = TRUE)
SEX  <- if (length(ARGS) >= 1 && nzchar(ARGS[1])) toupper(ARGS[1]) else ""
stopifnot(SEX %in% c("", "M", "F"))
SEX_SUF <- if (nzchar(SEX)) paste0("_", tolower(SEX)) else ""
KSUF <- "_k4"; SUF <- paste0(SEX_SUF, KSUF)

P7_DIR  <- "output/subtyping"
RPT_DIR <- file.path(P7_DIR, "reports")
CL_DIR  <- file.path(P7_DIR, "clusters")
dir.create(RPT_DIR, showWarnings = FALSE, recursive = TRUE); dir.create(CL_DIR, showWarnings = FALSE, recursive = TRUE)
while (sink.number() > 0) sink()
sink(file.path(RPT_DIR, sprintf("crosspop_log%s.txt", SUF)), split = TRUE)
on.exit(sink(), add = TRUE)

cat("=== held-out cross-population validation (G2, G3) ===\n")
cat(sprintf("sex: %s | solution: %s\n", if (nzchar(SEX)) SEX else "none", KSUF))

# ---- [1] G1 loadings + centroids --------------------------------------------
W_list <- readRDS(file.path(P7_DIR, "mofa", sprintf("mofa_loadings_g1%s.rds", SEX_SUF)))
cl_path <- file.path(CL_DIR, sprintf("cluster_assignments_g1%s%s.rds", SEX_SUF, KSUF))
if (!file.exists(cl_path)) cl_path <- file.path(CL_DIR, sprintf("cluster_assignments_g1%s.rds", SEX_SUF))
g1_cl <- readRDS(cl_path)
factors_used <- g1_cl$factors_used; Z_g1_active <- g1_cl$Z_active
labels_g1 <- g1_cl$final_labels; K <- g1_cl$final_k
cat(sprintf("  G1 active factors: %d | k=%d (%s)\n", length(factors_used), K, g1_cl$final_method))
centroids_g1 <- matrix(NA_real_, K, length(factors_used), dimnames = list(seq_len(K), factors_used))
for (cl in seq_len(K)) centroids_g1[cl, ] <- colMeans(Z_g1_active[labels_g1 == cl, , drop = FALSE])

# ---- helpers ----------------------------------------------------------------
project_factors <- function(X_list, W_list, factors_used,
                            views_use = c("olink", "nmr", "clinical_cont", "clinical_bin"), ridge = 0.01) {
  views_use <- intersect(intersect(views_use, names(W_list)), names(X_list))
  if (length(views_use) == 0) stop("no usable Gaussian views overlap")
  n <- ncol(X_list[[views_use[1]]]); K_full <- ncol(W_list[[views_use[1]]])
  Z_new <- matrix(NA_real_, n, K_full, dimnames = list(NULL, colnames(W_list[[views_use[1]]])))
  for (i in seq_len(n)) {
    A <- ridge * diag(K_full); b <- rep(0, K_full); contrib <- 0L
    for (v in views_use) {
      x_i <- X_list[[v]][, i]; obs <- which(!is.na(x_i))
      if (length(obs) < 5) next
      W_obs <- W_list[[v]][obs, , drop = FALSE]
      A <- A + crossprod(W_obs); b <- b + as.numeric(crossprod(W_obs, x_i[obs])); contrib <- contrib + 1L
    }
    if (contrib == 0L) next
    Z_new[i, ] <- as.numeric(solve(A, b))
  }
  Z_new[, factors_used, drop = FALSE]
}
nearest_centroid <- function(Z, centroids) {
  D <- as.matrix(dist(rbind(centroids, Z)))[seq_len(nrow(centroids)), -seq_len(nrow(centroids)), drop = FALSE]
  apply(D, 2, which.min)
}
hungarian_align <- function(labs_new, centroids_new, centroids_g1) {
  K <- nrow(centroids_g1); cost <- matrix(NA_real_, K, K)
  for (a in seq_len(K)) for (b in seq_len(K)) cost[a, b] <- sqrt(sum((centroids_new[a, ] - centroids_g1[b, ])^2))
  perm <- as.integer(clue::solve_LSAP(cost)); perm[labs_new]
}

# ---- per-cohort projection --------------------------------------------------
prevalence_rows <- list(); concordance_rows <- list(); distance_rows <- list()
for (grp in c("g2", "g3")) {
  in_path <- file.path(P7_DIR, "inputs", sprintf("clustering_input_%s%s.rds", grp, SEX_SUF))
  if (!file.exists(in_path)) { cat("  SKIP (missing):", in_path, "\n"); next }
  X_list <- readRDS(in_path); eids <- attr(X_list, "eid")
  cat(sprintf("\n[%s] projecting %d cases\n", toupper(grp), attr(X_list, "n_cases")))
  Z_new <- project_factors(X_list, W_list, factors_used)
  ok <- complete.cases(Z_new)
  cat(sprintf("  complete projections: %d / %d\n", sum(ok), nrow(Z_new)))
  Z_new_ok <- Z_new[ok, , drop = FALSE]; eids_ok <- eids[ok]
  if (nrow(Z_new_ok) < K * 5) { cat("  too few projected samples; skipping\n"); next }

  proj_labs <- nearest_centroid(Z_new_ok, centroids_g1)
  centroid_dist <- vapply(seq_len(nrow(Z_new_ok)),
                          function(i) sqrt(sum((Z_new_ok[i, ] - centroids_g1[proj_labs[i], ])^2)), numeric(1))
  distance_rows[[grp]] <- data.table(cohort = grp, cluster = proj_labs, centroid_distance = centroid_dist)

  km_new <- kmeans(Z_new_ok, centers = K, nstart = 25, iter.max = 100)
  indep_labs <- hungarian_align(km_new$cluster, km_new$centers, centroids_g1)
  sil_proj <- mean(silhouette(proj_labs, dist(Z_new_ok))[, "sil_width"])
  ari <- mclust::adjustedRandIndex(proj_labs, indep_labs)

  saveRDS(list(cohort = grp, eid = eids_ok, Z_active_proj = Z_new_ok, centroids_g1 = centroids_g1,
               proj_labels = proj_labs, indep_labels = indep_labs, factors_used = factors_used,
               centroid_distance = centroid_dist, silhouette_proj = sil_proj, ari_proj_vs_indep = ari,
               final_k = K, source_g1_cluster_file = cl_path),
          file.path(CL_DIR, sprintf("cluster_assignments_projected_%s%s.rds", grp, SUF)))

  prevalence_rows[[grp]] <- as.data.table(table(proj_labs))[, .(cohort = grp,
    cluster = as.integer(as.character(proj_labs)), n = N, pct = round(100 * N / sum(N), 2))]
  concordance_rows[[grp]] <- data.table(cohort = grp, n = nrow(Z_new_ok), n_total_input = nrow(Z_new),
    pct_complete_projection = round(100 * mean(ok), 2), silhouette_projected = round(sil_proj, 3),
    ari_projected_vs_indep = round(ari, 3), centroid_dist_median = round(median(centroid_dist), 3),
    centroid_dist_p95 = round(quantile(centroid_dist, 0.95), 3))
}

# ---- G1 reference + write summaries -----------------------------------------
prevalence_rows[["g1"]] <- as.data.table(table(labels_g1))[, .(cohort = "g1",
  cluster = as.integer(as.character(labels_g1)), n = N, pct = round(100 * N / sum(N), 2))]
prev_dt <- rbindlist(prevalence_rows, use.names = TRUE)[order(cohort, cluster)]
fwrite(prev_dt, file.path(RPT_DIR, sprintf("crosspop_subtype_prevalence%s.csv", SUF)))
cat("\nsubtype prevalence by cohort:\n"); print(prev_dt)
if (length(concordance_rows) > 0) {
  conc_dt <- rbindlist(concordance_rows, use.names = TRUE)
  fwrite(conc_dt, file.path(RPT_DIR, sprintf("crosspop_concordance_summary%s.csv", SUF)))
  cat("\nconcordance summary:\n"); print(conc_dt)
}
if (length(distance_rows) > 0)
  fwrite(rbindlist(distance_rows, use.names = TRUE), file.path(RPT_DIR, sprintf("crosspop_centroid_distance%s.csv", SUF)))

# ---- acceptance heuristic ---------------------------------------------------
cat("\n=== acceptance heuristic ===\n")
if (exists("conc_dt")) for (i in seq_len(nrow(conc_dt))) {
  r <- conc_dt[i]
  cat(sprintf("  %s: ARI=%.3f (>=0.5? %s), sil=%.3f (>=0.05? %s)\n", r$cohort, r$ari_projected_vs_indep,
              ifelse(r$ari_projected_vs_indep >= 0.50, "OK", "WEAK"), r$silhouette_projected,
              ifelse(r$silhouette_projected >= 0.05, "OK", "WEAK")))
}
cat("\n=== done ===\n")
