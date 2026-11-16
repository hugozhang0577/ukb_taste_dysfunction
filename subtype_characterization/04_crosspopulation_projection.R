# Cross-population projection (descriptive): group-2 and group-3 cases of one sex projected onto the group-1 MOFA+
# factors (ridge least squares on the Gaussian-view loadings) and assigned to the nearest group-1 subtype centroid;
# subtype prevalence, centroid distances and cluster separation per cohort.
library(data.table)
library(cluster)
library(clue)
library(mclust)
setwd(Sys.getenv("PROJECT_DIR", "."))
set.seed(20260413)

SEX_SUF <- paste0("_", tolower(commandArgs(trailingOnly = TRUE)[1]))   # Usage: Rscript subtype_characterization/04_crosspopulation_projection.R M|F
SUF     <- paste0(SEX_SUF, "_k4")
RPT_DIR <- "output/subtyping/reports"; CL_DIR <- "output/subtyping/clusters"

W_list  <- readRDS(sprintf("output/subtyping/mofa/mofa_loadings_g1%s.rds", SEX_SUF))
cl_path <- file.path(CL_DIR, sprintf("cluster_assignments_g1%s.rds", SUF))
g1_cl   <- readRDS(cl_path)
factors_used <- g1_cl$factors_used; labels_g1 <- g1_cl$final_labels; K <- g1_cl$final_k
centroids_g1 <- matrix(NA_real_, K, length(factors_used), dimnames = list(seq_len(K), factors_used))
for (cl in seq_len(K)) centroids_g1[cl, ] <- colMeans(g1_cl$Z_active[labels_g1 == cl, , drop = FALSE])

# factor scores of new samples: ridge-regularized least squares on the observed features of each Gaussian view
project_factors <- function(X_list, W_list, factors_used,
                            views_use = c("olink", "nmr", "clinical_cont", "clinical_bin"), ridge = 0.01) {
  views_use <- intersect(intersect(views_use, names(W_list)), names(X_list))
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
    if (contrib > 0L) Z_new[i, ] <- as.numeric(solve(A, b))
  }
  Z_new[, factors_used, drop = FALSE]
}
nearest_centroid <- function(Z, centroids) {
  D <- as.matrix(dist(rbind(centroids, Z)))[seq_len(nrow(centroids)), -seq_len(nrow(centroids)), drop = FALSE]
  apply(D, 2, which.min)
}
hungarian_align <- function(labs_new, centroids_new, centroids_g1) {   # match de novo clusters to group-1 subtypes
  K <- nrow(centroids_g1); cost <- matrix(NA_real_, K, K)
  for (a in seq_len(K)) for (b in seq_len(K)) cost[a, b] <- sqrt(sum((centroids_new[a, ] - centroids_g1[b, ])^2))
  as.integer(solve_LSAP(cost))[labs_new]
}

prevalence_rows <- list(); concordance_rows <- list(); distance_rows <- list()
for (grp in c("g2", "g3")) {
  X_list <- readRDS(sprintf("output/subtyping/inputs/clustering_input_%s%s.rds", grp, SEX_SUF))
  Z_new  <- project_factors(X_list, W_list, factors_used)
  ok <- complete.cases(Z_new); Z_ok <- Z_new[ok, , drop = FALSE]

  proj_labs <- nearest_centroid(Z_ok, centroids_g1)
  centroid_dist <- vapply(seq_len(nrow(Z_ok)), function(i) sqrt(sum((Z_ok[i, ] - centroids_g1[proj_labs[i], ])^2)), numeric(1))
  km_new <- kmeans(Z_ok, centers = K, nstart = 25, iter.max = 100)   # de novo k-means in the projected space
  indep_labs <- hungarian_align(km_new$cluster, km_new$centers, centroids_g1)
  sil_proj <- mean(silhouette(proj_labs, dist(Z_ok))[, "sil_width"])
  ari <- adjustedRandIndex(proj_labs, indep_labs)
  saveRDS(list(cohort = grp, eid = attr(X_list, "eid")[ok], Z_active_proj = Z_ok, centroids_g1 = centroids_g1,
               proj_labels = proj_labs, indep_labels = indep_labs, factors_used = factors_used,
               centroid_distance = centroid_dist, silhouette_proj = sil_proj, ari_proj_vs_indep = ari,
               final_k = K, source_g1_cluster_file = cl_path),
          file.path(CL_DIR, sprintf("cluster_assignments_projected_%s%s.rds", grp, SUF)))

  distance_rows[[grp]] <- data.table(cohort = grp, cluster = proj_labs, centroid_distance = centroid_dist)
  prevalence_rows[[grp]] <- as.data.table(table(proj_labs))[, .(cohort = grp,
    cluster = as.integer(as.character(proj_labs)), n = N, pct = round(100 * N / sum(N), 2))]
  concordance_rows[[grp]] <- data.table(cohort = grp, n = nrow(Z_ok), n_total_input = nrow(Z_new),
    pct_complete_projection = round(100 * mean(ok), 2), silhouette_projected = round(sil_proj, 3),
    ari_projected_vs_indep = round(ari, 3), centroid_dist_median = round(median(centroid_dist), 3),
    centroid_dist_p95 = round(quantile(centroid_dist, 0.95), 3))
}
prevalence_rows[["g1"]] <- as.data.table(table(labels_g1))[, .(cohort = "g1",
  cluster = as.integer(as.character(labels_g1)), n = N, pct = round(100 * N / sum(N), 2))]
fwrite(rbindlist(prevalence_rows, use.names = TRUE)[order(cohort, cluster)],
       file.path(RPT_DIR, sprintf("crosspop_subtype_prevalence%s.csv", SUF)))
fwrite(rbindlist(concordance_rows, use.names = TRUE), file.path(RPT_DIR, sprintf("crosspop_concordance_summary%s.csv", SUF)))
fwrite(rbindlist(distance_rows, use.names = TRUE), file.path(RPT_DIR, sprintf("crosspop_centroid_distance%s.csv", SUF)))
