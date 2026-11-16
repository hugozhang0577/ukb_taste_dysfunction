# One representative per correlation cluster that holds an FDR-significant measure: the CE-IVD
# certified member with the lead's direction and the smallest P, else the lead (smallest P).
library(data.table)
source("mwas/mwas_config.R")
setwd(Sys.getenv("PROJECT_DIR", "."))

FDR_ALPHA <- 0.05

# Volcano significance curve: -log10 P > 1.2 + 0.04 / (|log2 OR| - 0.02)
hyperbolic_threshold <- function(x, x0 = 0.02, cc = 0.04, ymin = 1.2)
  ifelse(abs(x) <= x0, Inf, ymin + cc / (abs(x) - x0))

raw      <- fread("output/mwas/primary/primary_20cov.csv")
clusters <- fread("output/mwas/effective_tests/metabolite_clusters.csv")
clusters <- unique(clusters[, .(metabolite, cluster_id, cluster_size)], by = "metabolite")
dt <- merge(clusters, raw, by.x = "metabolite", by.y = "protein")
dt[, pval_fdr_global := p.adjust(pval, method = "BH")]
dt[, is_sig_global := pval_fdr_global < FDR_ALPHA]

reps <- dt[!is.na(pval), {
  sig <- .SD[is_sig_global == TRUE]
  if (nrow(sig) == 0) NULL else {
    lead <- sig[which.min(pval)]
    ce   <- sig[metabolite %in% CE_37 & sign(log(or)) == sign(log(lead$or))]
    if (nrow(ce) > 0) ce[which.min(pval)] else lead
  }
}, by = cluster_id]

setnames(reps, "metabolite", "protein")
reps[, `:=`(neg_log10p = -log10(pval), log2or = log2(or), is_ce = protein %in% CE_37)]
reps[, sig_hyperbolic := neg_log10p > hyperbolic_threshold(log2or) &
       is.finite(hyperbolic_threshold(log2or))]
reps[, sig_fdr := is_sig_global]
front <- c("cluster_id", "protein", "is_ce", "or", "or_lower", "or_upper", "log2or", "pval",
           "neg_log10p", "sig_hyperbolic", "sig_fdr")
setcolorder(reps, c(front, setdiff(names(reps), front)))

# Copied to input/assoc_results/mwas_cluster_representatives.csv for feature selection
fwrite(reps, "output/mwas/effective_tests/cluster_representatives.csv")
