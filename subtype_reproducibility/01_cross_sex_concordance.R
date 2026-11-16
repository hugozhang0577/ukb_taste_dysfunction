# Cross-sex concordance: cosine similarity between male and female k = 4 cluster centroids on the union of
# top discriminating features (z-scored within sex), with greedy one-to-one matching of the clusters.
library(data.table)
setwd(Sys.getenv("PROJECT_DIR", "."))
set.seed(20260413)
RPT_DIR <- "output/subtyping/reports"

cl_m  <- readRDS("output/subtyping/clusters/cluster_assignments_g1_m_k4.rds")
cl_f  <- readRDS("output/subtyping/clusters/cluster_assignments_g1_f_k4.rds")
cases <- as.data.table(readRDS("output/ml_ready/group1_full.rds"))[taste_2w_strict == 1]
feats <- intersect(unique(c(fread(file.path(RPT_DIR, "fig_tables/table_cluster_topfeatures_m_k4.csv"))$feature_id,
                            fread(file.path(RPT_DIR, "fig_tables/table_cluster_topfeatures_f_k4.csv"))$feature_id)),
                   names(cases))

# z-score within sex
M <- as.matrix(cases[, feats, with = FALSE]); storage.mode(M) <- "numeric"; Mz <- M * NA_real_
for (s in c(0L, 1L)) {
  rows <- cases$sex == s; Ms <- M[rows, , drop = FALSE]
  mu <- apply(Ms, 2, mean, na.rm = TRUE); sdv <- apply(Ms, 2, sd, na.rm = TRUE); sdv[!is.finite(sdv) | sdv == 0] <- 1
  Mz[rows, ] <- sweep(sweep(Ms, 2, mu, "-"), 2, sdv, "/")
}
Mz[!is.finite(Mz)] <- NA; rownames(Mz) <- cases$eid

centroid <- function(cl) {
  eids <- as.integer(rownames(cl$Z_active)); lab <- cl$final_labels
  t(vapply(seq_len(cl$final_k), function(k) {
    r <- match(eids[lab == k], as.integer(rownames(Mz))); colMeans(Mz[r[!is.na(r)], , drop = FALSE], na.rm = TRUE)
  }, numeric(ncol(Mz))))
}
cosine <- function(a, b) { ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 5) return(NA_real_); sum(a[ok] * b[ok]) / (sqrt(sum(a[ok]^2)) * sqrt(sum(b[ok]^2))) }
C_m <- centroid(cl_m); C_f <- centroid(cl_f); K <- cl_m$final_k
S <- matrix(NA_real_, K, K, dimnames = list(paste0("M_c", 1:K), paste0("F_c", 1:K)))
for (i in 1:K) for (j in 1:K) S[i, j] <- cosine(C_m[i, ], C_f[j, ])
fwrite(as.data.table(S, keep.rownames = "M_cluster"), file.path(RPT_DIR, "concordance_matrix.csv"))

# greedy matching: repeatedly take the most similar remaining male-female pair
remaining_m <- 1:K; remaining_f <- 1:K; matches <- list(); S_work <- S
while (length(remaining_m) > 0) {
  idx <- which(S_work == max(S_work, na.rm = TRUE), arr.ind = TRUE)[1, ]
  matches[[length(matches) + 1]] <- data.table(M_cluster = remaining_m[idx[1]], F_cluster = remaining_f[idx[2]],
                                               cosine = S_work[idx[1], idx[2]])
  S_work <- S_work[-idx[1], -idx[2], drop = FALSE]
  remaining_m <- remaining_m[-idx[1]]; remaining_f <- remaining_f[-idx[2]]
}
matched <- rbindlist(matches)[order(-cosine)]
matched[, rank := c("best", "2nd", "3rd", "4th")[seq_len(.N)]]
fwrite(matched, file.path(RPT_DIR, "matched_pairs.csv"))
