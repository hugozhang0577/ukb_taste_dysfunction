#!/usr/bin/env Rscript
# =============================================================================
# 01_cross_sex_concordance.R  (subtyping)
#
# Test whether the k=4 sex-stratified clusters in males and females correspond
# to the same subtypes. Because MOFA was fit independently per sex, factors are
# not directly comparable; concordance is computed in the shared clinical
# feature space (the original variables), not in factor space.
#
# Method: common feature set (union of top-20 discriminators across all 8
# clusters); within-sex z-score (each sex standardised to its own mean/sd, so
# sex-baseline features do not push every cross-sex pair orthogonal); per-cluster
# centroids; cosine similarity (4x4); greedy best matching of M to F clusters.
#
# Input:  subtyping/clusters/cluster_assignments_g1_{m,f}_k4.rds +
#         subtyping/reports/fig_tables/table_cluster_topfeatures_{m,f}_k4.csv +
#         data/ml_ready/group1_full.rds
# Output: subtyping/reports/concordance_matrix.csv + matched_pairs.csv + heatmap
# =============================================================================

rm(list = ls()); gc()
suppressPackageStartupMessages({ library(data.table); library(pheatmap) })

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR); set.seed(20260413)

P7_DIR  <- "output/subtyping"
RPT_DIR <- file.path(P7_DIR, "reports")
TBL_DIR <- file.path(RPT_DIR, "fig_tables")
FIG_DIR <- file.path(RPT_DIR, "fig")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
while (sink.number() > 0) sink()
sink(file.path(RPT_DIR, "concordance_log.txt"), split = TRUE)
on.exit(sink(), add = TRUE)

cat("=== M vs F cluster concordance (k=4) ===\n")

# ---- [1] load ---------------------------------------------------------------
cl_m <- readRDS(file.path(P7_DIR, "clusters", "cluster_assignments_g1_m_k4.rds"))
cl_f <- readRDS(file.path(P7_DIR, "clusters", "cluster_assignments_g1_f_k4.rds"))
cases <- as.data.table(readRDS("output/ml_ready/group1_full.rds"))[taste_2w_strict == 1]
eid_m <- as.integer(rownames(cl_m$Z_active)); eid_f <- as.integer(rownames(cl_f$Z_active))
m_lab <- setNames(cl_m$final_labels, eid_m); f_lab <- setNames(cl_f$final_labels, eid_f)
cat(sprintf("  male n=%d k=%d | female n=%d k=%d\n", length(eid_m), cl_m$final_k, length(eid_f), cl_f$final_k))

# ---- [2] common feature set -------------------------------------------------
top_m <- fread(file.path(TBL_DIR, "table_cluster_topfeatures_m_k4.csv"))
top_f <- fread(file.path(TBL_DIR, "table_cluster_topfeatures_f_k4.csv"))
common_feats <- intersect(unique(c(top_m$feature_id, top_f$feature_id)), names(cases))
cat(sprintf("[2] common features: %d\n", length(common_feats)))

# ---- [3] within-sex z-score -------------------------------------------------
M <- as.matrix(cases[, common_feats, with = FALSE]); storage.mode(M) <- "numeric"
sex_vec <- cases$sex; Mz <- M * NA_real_
for (s in c(0L, 1L)) {
  rows <- sex_vec == s; Ms <- M[rows, , drop = FALSE]
  mu <- apply(Ms, 2, mean, na.rm = TRUE); sd <- apply(Ms, 2, sd, na.rm = TRUE); sd[!is.finite(sd) | sd == 0] <- 1
  Mz[rows, ] <- sweep(sweep(Ms, 2, mu, "-"), 2, sd, "/")
}
Mz[!is.finite(Mz)] <- NA; rownames(Mz) <- cases$eid

# ---- [4] centroids ----------------------------------------------------------
centroid <- function(eids, lab, Mz, k) {
  out <- matrix(NA_real_, nrow = k, ncol = ncol(Mz), dimnames = list(paste0("c", 1:k), colnames(Mz)))
  for (cl in 1:k) {
    rows <- match(eids[lab == cl], as.integer(rownames(Mz))); rows <- rows[!is.na(rows)]
    out[cl, ] <- colMeans(Mz[rows, , drop = FALSE], na.rm = TRUE)
  }
  out
}
C_m <- centroid(eid_m, m_lab, Mz, cl_m$final_k); C_f <- centroid(eid_f, f_lab, Mz, cl_f$final_k)

# ---- [5] cosine similarity --------------------------------------------------
cosine <- function(a, b) { ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 5) return(NA_real_); sum(a[ok] * b[ok]) / (sqrt(sum(a[ok]^2)) * sqrt(sum(b[ok]^2))) }
K <- cl_m$final_k
S <- matrix(NA_real_, K, K, dimnames = list(paste0("M_c", 1:K), paste0("F_c", 1:K)))
for (i in 1:K) for (j in 1:K) S[i, j] <- cosine(C_m[i, ], C_f[j, ])
cat("[5] cosine similarity (M rows x F cols):\n"); print(round(S, 3))
fwrite(as.data.table(S, keep.rownames = "M_cluster"), file.path(RPT_DIR, "concordance_matrix.csv"))

# ---- [6] greedy best matching -----------------------------------------------
remaining_m <- 1:K; remaining_f <- 1:K; matches <- list(); S_work <- S
while (length(remaining_m) > 0) {
  idx <- which(S_work == max(S_work, na.rm = TRUE), arr.ind = TRUE)[1, ]
  matches[[length(matches) + 1]] <- list(M_cluster = remaining_m[idx[1]], F_cluster = remaining_f[idx[2]],
                                         cosine = S_work[idx[1], idx[2]])
  S_work <- S_work[-idx[1], -idx[2], drop = FALSE]
  remaining_m <- remaining_m[-idx[1]]; remaining_f <- remaining_f[-idx[2]]
  if (is.null(dim(S_work))) S_work <- matrix(S_work, 1, 1)
}
matched <- rbindlist(lapply(matches, as.data.table))[order(-cosine)]
matched[, rank := c("best", "2nd", "3rd", "4th")[seq_len(.N)]]
fwrite(matched, file.path(RPT_DIR, "matched_pairs.csv"))
cat("\n[6] best matched pairs:\n"); print(matched)

# ---- [7] heatmap ------------------------------------------------------------
png(file.path(FIG_DIR, "concordance_heatmap.png"), width = 1200, height = 1000, res = 180)
pheatmap(S, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = round(S, 2), number_color = "black",
         color = colorRampPalette(c("steelblue", "white", "firebrick"))(50), breaks = seq(-1, 1, length.out = 51),
         main = "Male vs Female cluster cosine similarity (k=4)", fontsize_number = 11)
dev.off()
cat("\n=== done ===\n")
