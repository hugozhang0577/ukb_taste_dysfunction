#!/usr/bin/env Rscript
# =============================================================================
# 06_umap_embedding.R  (two-dimensional view of the factor space)
# =============================================================================
#
# UMAP of the active MOFA factors, one embedding per sex, with each case carrying
# its k-means subtype. Figure 4B and 4C draw these coordinates.
#
# This is a display step, not an analysis step. The subtypes were found by
# k-means in the full active-factor space, NOT in these two dimensions, so the
# embedding is a way of looking at that space and not evidence about it. UMAP
# preserves local neighbourhoods and distorts global distances, so the gaps and
# the relative sizes of the visible islands should not be read quantitatively;
# separation in the plot is neither necessary nor sufficient for a real cluster.
# The quantitative statements live in 03_cluster_on_factors.R (silhouette,
# bootstrap Jaccard) and in ../subtype_reproducibility/.
#
# The embedding is seeded, so re-running reproduces the same picture. Changing
# the seed rotates and re-arranges it without changing the partition.
#
# Usage: Rscript 06_umap_embedding.R          (both sexes)
# Input:  output/subtyping/clusters/cluster_assignments_g1_{m,f}_k4.rds
#           uses $Z_active (active factors) and $final_labels (cluster id)
# Output: output/subtyping/reports/fig_tables/umap_coords_{m,f}_k4.csv
# =============================================================================

CODE_DIR <- Sys.getenv("CODE_DIR", unset = getwd())
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
suppressPackageStartupMessages({ library(data.table); library(uwot) })

source(file.path(CODE_DIR, "_subtype_map.R"))   # SUBTYPE_MAP

SEED <- 20260414L
set.seed(SEED)

CL      <- "output/subtyping/clusters"
OUT_DIR <- "output/subtyping/reports/fig_tables"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# n_neighbors / min_dist are set for a few thousand points: large enough that the
# embedding reflects broad structure rather than sampling noise.
N_NEIGHBORS <- 30
MIN_DIST    <- 0.30

embed_one_sex <- function(sex) {
  stopifnot(sex %in% c("m", "f"))
  f <- file.path(CL, sprintf("cluster_assignments_g1_%s_k4.rds", sex))
  if (!file.exists(f)) stop("cluster assignments not found: ", f)
  cl <- readRDS(f)

  Z      <- cl$Z_active
  labels <- cl$final_labels
  if (nrow(Z) != length(labels))
    stop(sprintf("sex=%s: %d factor rows but %d labels", sex, nrow(Z), length(labels)))
  cat(sprintf("[%s] n = %d, active factors = %d\n", sex, nrow(Z), ncol(Z)))

  um <- uwot::umap(Z, n_neighbors = N_NEIGHBORS, min_dist = MIN_DIST,
                   metric = "euclidean", n_components = 2,
                   seed = SEED, verbose = FALSE)

  letter <- SUBTYPE_MAP[[sex]][as.character(as.integer(labels))]
  if (anyNA(letter)) stop("unmapped cluster id -- check _subtype_map.R")

  dt <- data.table(eid     = rownames(Z),
                   cluster = as.integer(labels),
                   subtype = unname(letter),
                   UMAP1   = um[, 1],
                   UMAP2   = um[, 2])
  print(dt[, .N, by = subtype][order(subtype)])

  out <- file.path(OUT_DIR, sprintf("umap_coords_%s_k4.csv", sex))
  fwrite(dt, out)
  # Export and dx upload to RAP  (Figure 4B/4C read these coordinates)
  cat("  wrote", out, "\n")
  invisible(dt)
}

invisible(lapply(c("m", "f"), embed_one_sex))
cat("done\n")
