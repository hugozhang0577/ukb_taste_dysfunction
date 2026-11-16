# Two-dimensional UMAP embedding of the active MOFA+ factors for each sex (30 neighbors, min_dist 0.3),
# with the subtype labels, for display.
library(data.table)
library(uwot)
source("unsupervised_subtyping/_subtype_map.R")
setwd(Sys.getenv("PROJECT_DIR", "."))
SEED <- 20260414L
set.seed(SEED)

for (sx in c("m", "f")) {
  cl  <- readRDS(sprintf("output/subtyping/clusters/cluster_assignments_g1_%s_k4.rds", sx))
  lab <- as.integer(cl$final_labels)
  um  <- umap(cl$Z_active, n_neighbors = 30, min_dist = 0.30, metric = "euclidean", n_components = 2,
              seed = SEED, verbose = FALSE)
  fwrite(data.table(eid = rownames(cl$Z_active), cluster = lab, subtype = unname(SUBTYPE_MAP[[sx]][as.character(lab)]),
                    UMAP1 = um[, 1], UMAP2 = um[, 2]),
         sprintf("output/subtyping/reports/fig_tables/umap_coords_%s_k4.csv", sx))
}
