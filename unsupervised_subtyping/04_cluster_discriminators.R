# Each k = 4 cluster vs the other three on every candidate feature: Wilcoxon rank-sum test (continuous) or
# Fisher's exact test (binary), BH FDR within cluster and feature type; top-20 discriminators per cluster.
library(data.table)
setwd(Sys.getenv("PROJECT_DIR", "."))
set.seed(20260413)

SEX <- tolower(commandArgs(trailingOnly = TRUE)[1])   # Usage: Rscript unsupervised_subtyping/04_cluster_discriminators.R M|F
SUF <- sprintf("_%s_k4", SEX)
TBL_DIR <- "output/subtyping/reports/fig_tables"
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

cl_obj <- readRDS(sprintf("output/subtyping/clusters/cluster_assignments_g1%s.rds", SUF))
cases  <- as.data.table(readRDS("output/ml_ready/group1_full.rds"))
cases  <- cases[match(as.integer(rownames(cl_obj$Z_active)), cases$eid)][, cluster := factor(cl_obj$final_labels)]

manifest <- fread("output/feature_manifest/master_feature_manifest_final.csv")
classify_type <- function(src, vtype, dx) {
  if (!is.na(src) && src %in% c("PWAS", "MWAS")) return("continuous")
  if (!is.na(dx) && nzchar(as.character(dx))) return("binary")
  vt <- tolower(ifelse(is.na(vtype) | vtype == "", "unknown", vtype))
  if (vt %in% c("binary", "derived_binary")) "binary" else "continuous"
}
manifest[, ftype := mapply(classify_type, source_analysis, var_type, icd10)]
ftype_of  <- setNames(manifest$ftype, manifest$feature_id)
feats_all <- setdiff(intersect(manifest$feature_id, names(cases)),
                     manifest[var_role == "outcome_primary" | grepl("^PC[0-9]+$", feature_id), feature_id])

res <- list()
for (cl in levels(cases$cluster)) {
  in_cl <- cases$cluster == cl
  for (fid in feats_all) {
    x <- cases[[fid]]
    if (all(is.na(x))) next
    if (ftype_of[[fid]] == "continuous") {
      x1 <- x[in_cl]; x0 <- x[!in_cl]
      if (sum(!is.na(x1)) < 5 || sum(!is.na(x0)) < 5) next
      res[[length(res) + 1]] <- data.table(cluster = cl, feature_id = fid, type = "continuous",
        median_in = median(x1, na.rm = TRUE), median_out = median(x0, na.rm = TRUE),
        delta = median(x1, na.rm = TRUE) - median(x0, na.rm = TRUE),
        prev_in = NA_real_, prev_out = NA_real_, OR = NA_real_,
        p = tryCatch(wilcox.test(x1, x0, exact = FALSE)$p.value, error = function(e) NA_real_))
    } else {
      x <- as.integer(x); x1 <- x[in_cl]; x0 <- x[!in_cl]
      tab <- table(factor(as.integer(in_cl), levels = c(0, 1)), factor(x, levels = c(0, 1)))
      if (any(rowSums(tab) == 0) || any(colSums(tab) == 0)) next
      ft <- tryCatch(fisher.test(tab), error = function(e) NULL); if (is.null(ft)) next
      res[[length(res) + 1]] <- data.table(cluster = cl, feature_id = fid, type = "binary",
        median_in = NA_real_, median_out = NA_real_, delta = NA_real_,
        prev_in = mean(x1, na.rm = TRUE), prev_out = mean(x0, na.rm = TRUE),
        OR = unname(ft$estimate), p = ft$p.value)
    }
  }
}
prof <- rbindlist(res, fill = TRUE)
prof[, q := p.adjust(p, method = "BH"), by = .(cluster, type)]
prof <- prof[order(cluster, q)]
prof[type == "continuous", abs_effect := abs(delta)]
prof[type == "binary", abs_effect := abs(log2(pmax(OR, 1e-3)))]
fwrite(prof, sprintf("output/subtyping/reports/cluster_profiles%s.csv", SUF))

top <- prof[q < 0.05][order(cluster, -abs_effect), head(.SD, 20), by = cluster]
top <- merge(top, manifest[, .(feature_id, source_analysis, description)], by = "feature_id", all.x = TRUE)
setcolorder(top, c("cluster", "feature_id", "source_analysis", "type", "median_in", "median_out",
                   "delta", "prev_in", "prev_out", "OR", "p", "q", "abs_effect", "description"))
fwrite(top[order(cluster, -abs_effect)], file.path(TBL_DIR, sprintf("table_cluster_topfeatures%s.csv", SUF)))
