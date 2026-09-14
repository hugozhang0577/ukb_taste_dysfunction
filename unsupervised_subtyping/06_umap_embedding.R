#!/usr/bin/env Rscript
#> in   output/subtyping/clusters/cluster_assignments_g1_{}_k4.rds
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
  # Export and dx upload to RAP  (the embedding coordinates)
  cat("  wrote", out, "\n")
  invisible(dt)
}

invisible(lapply(c("m", "f"), embed_one_sex))
cat("done\n")
