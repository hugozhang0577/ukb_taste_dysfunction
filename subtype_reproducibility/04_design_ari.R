# Agreement with the primary k = 4 solution, per sex: adjusted Rand index (ARI) for the 20 seed re-runs and the
# design cells (variance-balanced weighting, broadened feature set, both), with a 10,000-permutation null for
# the design cells; label retention after greedy cluster matching; clustering-step bootstrap (500 resamples).
library(data.table)
library(mclust)
source("unsupervised_subtyping/_subtype_map.R")
setwd(Sys.getenv("PROJECT_DIR", "."))
set.seed(20260413)

K <- 4L; B_BOOT <- 500L; N_PERM <- 10000L
CL <- "output/subtyping/clusters"; RPT_DIR <- "output/subtyping/reports"
DESIGNS <- c("balanced", "broadened", "broadened_balanced")

labels_of <- function(path) {
  o <- readRDS(path); lb <- o$final_labels
  if (is.null(names(lb))) names(lb) <- as.character(if (!is.null(o$eid)) o$eid else rownames(o$Z_active))
  lb
}
retention <- function(ref, cmp) {   # share of participants keeping their cluster after greedy matching
  ids <- intersect(names(ref), names(cmp)); a <- ref[ids]; b <- cmp[ids]
  tt <- table(b, a); map <- setNames(rep(NA_character_, nrow(tt)), rownames(tt))
  while (any(!is.na(tt)) && any(tt > 0, na.rm = TRUE)) {
    ij <- which(tt == max(tt, na.rm = TRUE), arr.ind = TRUE)[1, ]
    map[rownames(tt)[ij[1]]] <- colnames(tt)[ij[2]]
    tt[ij[1], ] <- NA; tt[, ij[2]] <- NA
  }
  mean(map[as.character(b)] == as.character(a), na.rm = TRUE)
}
bootstrap_ari <- function(Z, labels) vapply(seq_len(B_BOOT), function(b) {
  idx <- sample.int(nrow(Z), nrow(Z), replace = TRUE)
  km <- tryCatch(kmeans(Z[idx, , drop = FALSE], centers = K, nstart = 25, iter.max = 100), error = function(e) NULL)
  if (is.null(km)) NA_real_ else adjustedRandIndex(km$cluster, labels[idx])
}, numeric(1))
res_row <- function(sex, comparison, n = NA_integer_, ari = NA_real_, ari_lo = NA_real_, ari_hi = NA_real_,
                    retention = NA_real_, null_median = NA_real_, null_p95 = NA_real_, status = "ok")
  data.table(sex, comparison, n, ari, ari_lo, ari_hi, retention, null_median, null_p95, status)

rows <- list(); ret_rows <- list(); null_rows <- list()
for (sex in c("m", "f")) {
  main_path <- file.path(CL, sprintf("cluster_assignments_g1_%s_k%d.rds", sex, K))
  main_lab  <- labels_of(main_path)

  # clustering step only: k-means refitted on bootstrap resamples of the primary factor scores
  bs <- bootstrap_ari(readRDS(main_path)$Z_active, main_lab)
  rows[[length(rows) + 1]] <- res_row(sex, "bootstrap_ceiling_clustering_step", length(main_lab),
    median(bs, na.rm = TRUE), quantile(bs, 0.025, na.rm = TRUE), quantile(bs, 0.975, na.rm = TRUE))

  # seed re-runs
  keep <- numeric(0)
  for (s in other_seeds(sex)) {
    lb  <- labels_of(file.path(CL, sprintf("cluster_assignments_g1_%s_reseed%d_k%d.rds", sex, s, K)))
    ids <- intersect(names(main_lab), names(lb)); keep <- c(keep, retention(main_lab, lb))
    rows[[length(rows) + 1]] <- res_row(sex, paste0("seed_", s), length(ids),
                                        adjustedRandIndex(main_lab[ids], lb[ids]), retention = keep[length(keep)])
  }
  ret_rows[[sex]] <- data.table(sex = sex, comparison = "seed_refits", n_runs = length(keep),
                                retention_mean = mean(keep), retention_min = min(keep), retention_max = max(keep))

  # design cells, against a permutation null
  for (d in DESIGNS) {
    f <- file.path(CL, sprintf("cluster_assignments_g1_%s_%s_k%d.rds", sex, d, K))
    if (!file.exists(f)) { rows[[length(rows) + 1]] <- res_row(sex, d, status = "labels_not_available"); next }
    lb  <- labels_of(f); ids <- intersect(names(main_lab), names(lb))
    nul <- vapply(seq_len(N_PERM), function(i) adjustedRandIndex(sample(main_lab[ids]), lb[ids]), numeric(1))
    rows[[length(rows) + 1]] <- res_row(sex, d, length(ids), adjustedRandIndex(main_lab[ids], lb[ids]),
      retention = retention(main_lab, lb), null_median = median(nul), null_p95 = quantile(nul, 0.95))
    null_rows[[length(null_rows) + 1]] <- data.table(sex = sex, comparison = d, ari_null = nul)
  }
}
fwrite(rbindlist(rows, fill = TRUE), file.path(RPT_DIR, "reproducibility_ari_summary.csv"))
fwrite(rbindlist(ret_rows), file.path(RPT_DIR, "reproducibility_retention.csv"))
fwrite(rbindlist(null_rows), file.path(RPT_DIR, "reproducibility_null_distribution.csv"))
