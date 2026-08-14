#!/usr/bin/env Rscript
# =============================================================================
# 04_design_ari.R  (agreement between the main solution and each design cell)
# =============================================================================
#
# Puts every reproducibility comparison on one Adjusted Rand Index scale at
# fixed k, and reports, per sex:
#
#   seed re-runs      the whole pipeline re-run under a different random seed
#                     (02_reseed_stability.R). This is the reproducibility
#                     statement in the manuscript.
#   design cells      the same pipeline under a different MOFA design choice
#                     (03_design_sensitivity.R)
#
# Two reference quantities are computed so an ARI can be read against something:
#
#   permutation null  the ARI two unrelated partitions of the same sizes reach
#                     by chance. The ARI is chance-corrected, so this sits near
#                     zero; it is computed rather than assumed.
#   bootstrap ceiling ARI of the clustering step against itself under resampling
#                     of participants, with the factor space held fixed. It is
#                     an upper bound for the *clustering step only*: it re-runs
#                     k-means, not the factor decomposition, so it is blind to
#                     instability in the decomposition and must not be read as
#                     an upper bound on the seed re-runs.
#
# Also writes the per-participant retention rate each comparison implies — the
# share of participants whose subtype is unchanged after optimal relabelling —
# which is what Figure 5's reproducibility panel shows. Retention is NOT
# chance-corrected: with four subtypes its random baseline is about 25%, so it
# must not be read on the same scale as an ARI.
#
# Usage: Rscript 04_design_ari.R [k]          (k defaults to 4)
# Input:  output/subtyping/clusters/cluster_assignments_g1_{m,f}*_k{k}.rds
# Output: output/subtyping/reports/reproducibility_ari_summary.csv
#         output/subtyping/reports/reproducibility_retention.csv
#         output/subtyping/reports/reproducibility_null_distribution.csv
# =============================================================================

rm(list = ls()); gc()
suppressPackageStartupMessages({
  library(data.table)
  library(mclust)      # adjustedRandIndex
})

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR); set.seed(20260413)

ARGS <- commandArgs(trailingOnly = TRUE)
K <- if (length(ARGS) >= 1 && nzchar(ARGS[1])) as.integer(ARGS[1]) else 4L

B_BOOT <- 500L      # bootstrap replicates for the clustering-step ceiling
N_PERM <- 10000L    # permutation replicates for the chance null

CL      <- "output/subtyping/clusters"
RPT_DIR <- "output/subtyping/reports"
dir.create(RPT_DIR, showWarnings = FALSE, recursive = TRUE)
while (sink.number() > 0) sink()
sink(file.path(RPT_DIR, "reproducibility_ari_log.txt"), split = TRUE)
on.exit({ while (sink.number() > 0) sink() }, add = TRUE)

cat("=== reproducibility: ARI against the main solution (k =", K, ") ===\n")
cat("bootstrap replicates:", B_BOOT, "| permutation replicates:", N_PERM, "\n\n")

# ---- helpers ----------------------------------------------------------------
labels_of <- function(path) {
  o <- readRDS(path)
  lb <- o$final_labels
  if (is.null(names(lb))) {
    ids <- if (!is.null(o$eid)) o$eid else rownames(o$Z_active)
    if (is.null(ids)) stop("no participant ids in ", path)
    names(lb) <- as.character(ids)
  }
  lb
}

# Share of participants keeping their subtype after optimally relabelling the
# comparison partition onto the reference (Hungarian-style greedy on the
# contingency table, which is exact enough at k = 4).
retention <- function(ref, cmp) {
  ids <- intersect(names(ref), names(cmp))
  a <- ref[ids]; b <- cmp[ids]
  tab <- table(b, a)
  map <- setNames(rep(NA_character_, nrow(tab)), rownames(tab))
  tt <- tab
  while (any(!is.na(tt)) && any(tt > 0, na.rm = TRUE)) {
    ij <- which(tt == max(tt, na.rm = TRUE), arr.ind = TRUE)[1, ]
    map[rownames(tt)[ij[1]]] <- colnames(tt)[ij[2]]
    tt[ij[1], ] <- NA; tt[, ij[2]] <- NA
  }
  mean(map[as.character(b)] == as.character(a), na.rm = TRUE)
}

bootstrap_ceiling <- function(Z, labels, B, k) {
  n <- nrow(Z); out <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    km <- tryCatch(suppressWarnings(
      kmeans(Z[idx, , drop = FALSE], centers = k, nstart = 25, iter.max = 100)),
      error = function(e) NULL)
    out[b] <- if (is.null(km)) NA_real_ else adjustedRandIndex(km$cluster, labels[idx])
  }
  out
}

permutation_null <- function(a, b, N) {
  out <- numeric(N)
  for (i in seq_len(N)) out[i] <- adjustedRandIndex(sample(a), b)
  out
}

# ---- the comparisons, per sex -----------------------------------------------
# Design cells are discovered rather than listed, so a cell whose labels are not
# on disk is reported as missing instead of silently dropping out of the table.
DESIGNS <- c(balanced           = "balanced weighting",
             broadened          = "broadened feature set",
             broadened_balanced = "broadened + balanced")

rows <- list(); ret_rows <- list(); null_rows <- list()

for (sex in c("m", "f")) {
  main_path <- file.path(CL, sprintf("cluster_assignments_g1_%s_k%d.rds", sex, K))
  if (!file.exists(main_path)) { cat("[skip]", sex, "- main solution not found\n"); next }
  main_obj <- readRDS(main_path)
  main_lab <- labels_of(main_path)
  cat("\n--- ", toupper(sex), ": n = ", length(main_lab), " ---\n", sep = "")

  # bootstrap ceiling for the clustering step
  Z <- main_obj$Z_active
  bs <- bootstrap_ceiling(Z, main_lab, B_BOOT, K)
  cat(sprintf("  clustering-step bootstrap ceiling: median ARI %.3f  [%.3f, %.3f]\n",
              median(bs, na.rm = TRUE), quantile(bs, 0.025, na.rm = TRUE),
              quantile(bs, 0.975, na.rm = TRUE)))
  rows[[length(rows) + 1]] <- data.table(
    sex = sex, comparison = "bootstrap_ceiling_clustering_step", n = length(main_lab),
    ari = median(bs, na.rm = TRUE), ari_lo = quantile(bs, 0.025, na.rm = TRUE),
    ari_hi = quantile(bs, 0.975, na.rm = TRUE), retention = NA_real_,
    null_median = NA_real_, null_p95 = NA_real_, status = "ok")

  # seed re-runs
  seed_files <- list.files(CL, full.names = TRUE,
    pattern = sprintf("^cluster_assignments_g1_%s_reseed[0-9]+_k%d\\.rds$", sex, K))
  cat("  seed re-runs found:", length(seed_files), "\n")
  for (f in seed_files) {
    seed <- sub(".*reseed([0-9]+)_k.*", "\\1", basename(f))
    lb   <- labels_of(f)
    ids  <- intersect(names(main_lab), names(lb))
    rows[[length(rows) + 1]] <- data.table(
      sex = sex, comparison = paste0("seed_", seed), n = length(ids),
      ari = adjustedRandIndex(main_lab[ids], lb[ids]),
      ari_lo = NA_real_, ari_hi = NA_real_,
      retention = retention(main_lab, lb),
      null_median = NA_real_, null_p95 = NA_real_, status = "ok")
  }

  # design cells
  for (d in names(DESIGNS)) {
    f <- file.path(CL, sprintf("cluster_assignments_g1_%s_%s_k%d.rds", sex, d, K))
    if (!file.exists(f)) {
      cat("  [missing]", d, "- no labels at", basename(f), "\n")
      rows[[length(rows) + 1]] <- data.table(
        sex = sex, comparison = d, n = NA_integer_, ari = NA_real_,
        ari_lo = NA_real_, ari_hi = NA_real_, retention = NA_real_,
        null_median = NA_real_, null_p95 = NA_real_, status = "labels_not_available")
      next
    }
    lb  <- labels_of(f)
    ids <- intersect(names(main_lab), names(lb))
    ari <- adjustedRandIndex(main_lab[ids], lb[ids])
    nul <- permutation_null(main_lab[ids], lb[ids], N_PERM)
    ret <- retention(main_lab, lb)
    cat(sprintf("  %-20s ARI %.3f | retention %.1f%% | null median %.4f\n",
                DESIGNS[[d]], ari, 100 * ret, median(nul)))
    rows[[length(rows) + 1]] <- data.table(
      sex = sex, comparison = d, n = length(ids), ari = ari,
      ari_lo = NA_real_, ari_hi = NA_real_, retention = ret,
      null_median = median(nul), null_p95 = quantile(nul, 0.95), status = "ok")
    null_rows[[length(null_rows) + 1]] <- data.table(sex = sex, comparison = d, ari_null = nul)
  }

  # retention of the seed re-runs, pooled: the figure's central number
  if (length(seed_files) > 0) {
    keep <- vapply(seed_files, function(f) retention(main_lab, labels_of(f)), numeric(1))
    cat(sprintf("  seed re-runs: mean retention %.1f%% over %d seeds\n",
                100 * mean(keep), length(keep)))
    ret_rows[[length(ret_rows) + 1]] <- data.table(
      sex = sex, comparison = "seed_reruns", n_runs = length(keep),
      retention_mean = mean(keep), retention_min = min(keep), retention_max = max(keep))
  }
}

summary_dt <- rbindlist(rows, fill = TRUE)
fwrite(summary_dt, file.path(RPT_DIR, "reproducibility_ari_summary.csv"))
if (length(ret_rows))  fwrite(rbindlist(ret_rows),  file.path(RPT_DIR, "reproducibility_retention.csv"))
if (length(null_rows)) fwrite(rbindlist(null_rows), file.path(RPT_DIR, "reproducibility_null_distribution.csv"))
# Export and dx upload to RAP  (the reproducibility summary and its null)

cat("\n=== summary ===\n")
print(summary_dt[, .(sex, comparison, n, ari = round(ari, 3),
                     retention = round(retention, 3), status)])
cat("\nARI is chance-corrected; retention is not (random baseline about 1/k).\n")
