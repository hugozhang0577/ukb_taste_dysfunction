# Clustering of group-1 cases on the active MOFA+ factors (total variance explained >= 6%): k-means and PAM
# for k = 2-8 and HDBSCAN, with bootstrap Jaccard stability (100 x 80% subsamples); the reported solution
# is k-means with k = 4.
library(data.table)
library(cluster)
library(dbscan)
library(fpc)
source("unsupervised_subtyping/_subtype_map.R")
setwd(Sys.getenv("PROJECT_DIR", "."))

SEX <- tolower(commandArgs(trailingOnly = TRUE)[1])   # Usage: Rscript unsupervised_subtyping/03_cluster_on_factors.R M|F
K   <- 4L
SUF <- sprintf("_%s_k%d", SEX, K)
set.seed(PRIMARY_SEED[[SEX]])
dir.create("output/subtyping/clusters", showWarnings = FALSE, recursive = TRUE)

Z <- readRDS(sprintf("output/subtyping/mofa/mofa_factors_g1_%s.rds", SEX))
activity <- fread(sprintf("output/subtyping/reports/factor_activity_%s.csv", SEX))
factors_keep <- activity[total_var >= 6, factor]
Z_act <- Z[, factors_keep, drop = FALSE]
N <- nrow(Z_act); D <- dist(Z_act)

# ---- sweep: k-means / PAM (k = 2-8), HDBSCAN (minPts 30-200)
sweep_row <- function(method, k, lab, sil, X = Z_act)
  data.table(method = method, k = k, silhouette = sil, CH = calinhara(X, lab), min_cluster_pct = min(table(lab)) / N)
sweep_rows <- list(); labels_store <- list()
for (k in 2:8) {
  km <- kmeans(Z_act, centers = k, nstart = 25, iter.max = 100)
  sweep_rows[[length(sweep_rows) + 1]] <- sweep_row("kmeans", k, km$cluster, mean(silhouette(km$cluster, D)[, "sil_width"]))
  labels_store[[sprintf("kmeans_k%d", k)]] <- km$cluster
  pm <- pam(Z_act, k = k, metric = "euclidean", pamonce = 5)
  sweep_rows[[length(sweep_rows) + 1]] <- sweep_row("pam", k, pm$clustering, pm$silinfo$avg.width)
  labels_store[[sprintf("pam_k%d", k)]] <- pm$clustering
}
for (mp in c(30, 60, 100, 200)) {
  cl <- hdbscan(Z_act, minPts = mp)$cluster; ok <- cl > 0; n_cl <- length(unique(cl[ok]))
  sweep_rows[[length(sweep_rows) + 1]] <- if (n_cl >= 2)
    sweep_row(sprintf("hdbscan_mp%d", mp), n_cl, cl[ok],
              if (sum(ok) > 10) mean(silhouette(cl[ok], dist(Z_act[ok, ]))[, "sil_width"]) else NA_real_, Z_act[ok, ])
  else data.table(method = sprintf("hdbscan_mp%d", mp), k = n_cl, silhouette = NA_real_, CH = NA_real_, min_cluster_pct = NA_real_)
  labels_store[[sprintf("hdbscan_mp%d", mp)]] <- cl
}
sweep <- rbindlist(sweep_rows, fill = TRUE)

# ---- bootstrap Jaccard for k-means / PAM solutions with min cluster >= 5% and silhouette > 0
boot_jaccard <- function(method, k, B = 100, frac = 0.8) {
  base <- if (method == "kmeans") kmeans(Z_act, centers = k, nstart = 25, iter.max = 100)$cluster
          else pam(Z_act, k = k, pamonce = 5)$clustering
  mean(vapply(seq_len(B), function(b) {
    idx <- sample(N, size = floor(frac * N), replace = FALSE)
    sub_lab <- if (method == "kmeans") kmeans(Z_act[idx, ], centers = k, nstart = 10, iter.max = 100)$cluster
               else pam(Z_act[idx, ], k = k, pamonce = 5)$clustering
    base_sub <- base[idx]
    mean(sapply(seq_len(k), function(cl) { a <- base_sub == cl
      max(sapply(seq_len(k), function(cl2) { b2 <- sub_lab == cl2; sum(a & b2) / max(1, sum(a | b2)) })) }))
  }, numeric(1)))
}
cands <- sweep[method %in% c("kmeans", "pam") & min_cluster_pct >= 0.05 & silhouette > 0]
cands[, bootstrap_jaccard := NA_real_]
for (i in seq_len(nrow(cands))) set(cands, i, "bootstrap_jaccard", boot_jaccard(cands$method[i], cands$k[i]))
sweep <- merge(sweep, cands[, .(method, k, bootstrap_jaccard)], by = c("method", "k"), all.x = TRUE)
fwrite(sweep, sprintf("output/subtyping/reports/cluster_sweep%s.csv", SUF))

final <- sweep[method == "kmeans" & k == K]
saveRDS(list(factors_used = factors_keep, Z_active = Z_act, sweep = sweep, all_labels = labels_store,
             final_method = final$method, final_k = final$k, final_labels = labels_store[[sprintf("kmeans_k%d", K)]],
             final_silhouette = final$silhouette, final_jaccard = final$bootstrap_jaccard),
        sprintf("output/subtyping/clusters/cluster_assignments_g1%s.rds", SUF))
