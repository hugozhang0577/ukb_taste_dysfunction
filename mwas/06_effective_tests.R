# Correlation structure of the NMR panel in the discovery cohort: effective number of independent
# tests (Li and Ji 2005) and complete-linkage clustering of the measures at |r| > 0.80.
library(data.table)
source("mwas/mwas_config.R")
setwd(Sys.getenv("PROJECT_DIR", "."))

OUT_DIR     <- "output/mwas/effective_tests"
CORR_CUTOFF <- 0.80
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

nmr  <- fread(NMR_FILES[["group1"]])
cols <- setdiff(names(nmr), "eid")
x    <- as.matrix(nmr[, ..cols])
for (j in seq_len(ncol(x))) x[is.na(x[, j]), j] <- median(x[, j], na.rm = TRUE)  # median-impute
cor_mat <- cor(x, use = "pairwise.complete.obs")

# Li and Ji: M_eff = sum over eigenvalues of I(lambda >= 1) + (lambda - floor(lambda))
ev       <- eigen(cor_mat, only.values = TRUE)$values
ent_liji <- round(sum(ifelse(ev >= 1, 1 + (ev - floor(ev)), ev)))

cl_id <- cutree(hclust(as.dist(1 - abs(cor_mat)), method = "complete"), h = 1 - CORR_CUTOFF)
cl_n  <- table(cl_id)
clusters <- data.table(metabolite = cols, cluster_id = cl_id[cols],
                       cluster_size = as.integer(cl_n[as.character(cl_id[cols])]),
                       is_ce = cols %in% CE_37)

mwas <- fread("output/mwas/primary/primary_20cov.csv")
clusters <- merge(clusters, mwas[, .(protein, beta, se, or, pval, pval_fdr)],
                  by.x = "metabolite", by.y = "protein", all.x = TRUE)

# Membership-table flag: the CE-IVD member with the smallest P, else the member with the smallest P
rep_of <- function(m, ce, p) {
  ok <- !is.na(p)
  if (!any(ok)) return(m[1])
  pool <- if (any(ok & ce)) ok & ce else ok
  m[pool][which.min(p[pool])]
}
reps <- clusters[, .(metabolite = rep_of(metabolite, is_ce, pval)), by = cluster_id]
clusters[, is_representative := metabolite %in% reps$metabolite]
setorder(clusters, cluster_id, pval)
fwrite(clusters, file.path(OUT_DIR, "metabolite_clusters.csv"))

r <- cor_mat[upper.tri(cor_mat)]
writeLines(c(sprintf("NMR measures: %d", length(cols)),
             sprintf("Pairs with |r| > 0.8: %.1f%%", 100 * mean(abs(r) > 0.8)),
             sprintf("Effective number of tests (Li and Ji): %d", ent_liji),
             sprintf("Clusters at |r| > %.2f: %d", CORR_CUTOFF, max(cl_id))),
           file.path(OUT_DIR, "ent_report.txt"))
