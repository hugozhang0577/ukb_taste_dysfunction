#!/usr/bin/env Rscript
# =============================================================================
# 04_cluster_discriminators.R  (subtyping)
#
# Take the cluster assignments from 03_cluster_on_factors and compute per-cluster descriptive stats +
# cluster-vs-rest discriminator tests across ALL features in the final manifest
# (not only the MOFA-input features):
#   continuous     : Wilcoxon rank-sum (cluster vs rest), median per cluster
#   binary/PheCode : Fisher exact (cluster vs rest), prevalence + OR
#   BH FDR within (cluster, feature_type)
#
# CLI: arg1 = ""/"M"/"F"; arg2 = forced k (optional).
#
# Input:  subtyping/clusters/cluster_assignments_g1*.rds + data/ml_ready/group1_full.rds
#         + manifest/master_feature_manifest_final.csv
# Output: subtyping/reports/cluster_profiles*.csv +
#         subtyping/reports/fig_tables/table_cluster_topfeatures*.csv + heatmap
# =============================================================================

rm(list = ls()); gc()
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(pheatmap) })

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR); set.seed(20260413)

ARGS    <- commandArgs(trailingOnly = TRUE)
SEX     <- if (length(ARGS) >= 1 && nzchar(ARGS[1])) toupper(ARGS[1]) else ""
FORCE_K <- if (length(ARGS) >= 2 && nzchar(ARGS[2])) as.integer(ARGS[2]) else NA_integer_
stopifnot(SEX %in% c("", "M", "F"))
SUF <- paste0(if (nzchar(SEX)) paste0("_", tolower(SEX)) else "", if (!is.na(FORCE_K)) paste0("_k", FORCE_K) else "")

ML_DIR  <- "output/ml_ready"
MF_DIR  <- "output/feature_manifest"
P7_DIR  <- "output/subtyping"
RPT_DIR <- file.path(P7_DIR, "reports")
TBL_DIR <- file.path(RPT_DIR, "fig_tables")
FIG_DIR <- file.path(RPT_DIR, "fig")
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE); dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
while (sink.number() > 0) sink()
sink(file.path(RPT_DIR, sprintf("characterization_log%s.txt", SUF)), split = TRUE)
on.exit(sink(), add = TRUE)

cat("=== cluster characterization ===\n")

# ---- [1] load + align -------------------------------------------------------
cl_obj <- readRDS(file.path(P7_DIR, "clusters", sprintf("cluster_assignments_g1%s.rds", SUF)))
labels <- cl_obj$final_labels; eids <- as.integer(rownames(cl_obj$Z_active))
stopifnot(length(labels) == length(eids))
cat(sprintf("  method=%s, k=%d\n", cl_obj$final_method, cl_obj$final_k)); print(table(labels))

cases <- as.data.table(readRDS(file.path(ML_DIR, "group1_full.rds")))[taste_2w_strict == 1]
if (nzchar(SEX)) cases <- cases[sex == if (SEX == "M") 1L else 0L]
stopifnot("eid" %in% names(cases))
m <- match(eids, cases$eid)
if (any(is.na(m))) stop(sprintf("%d cluster eids not in case data", sum(is.na(m))))
cases <- cases[m]; cases[, cluster := factor(labels)]
cat(sprintf("  aligned %d cases to %d clusters\n", nrow(cases), nlevels(cases$cluster)))

# ---- [2] classify features --------------------------------------------------
manifest <- fread(file.path(MF_DIR, "master_feature_manifest_final.csv"))
exclude_fids <- c(manifest[var_role == "outcome_primary", feature_id], manifest[grepl("^PC[0-9]+$", feature_id), feature_id])
classify_type <- function(src, vtype, phe) {
  if (!is.na(src) && src %in% c("PWAS", "MWAS")) return("continuous")
  if (!is.na(phe) && nzchar(as.character(phe))) return("binary")
  vt <- tolower(ifelse(is.na(vtype) | vtype == "", "unknown", vtype))
  if (vt %in% c("binary", "derived_binary")) return("binary")
  "continuous"
}
manifest[, ftype := mapply(classify_type, source_analysis, var_type, phecode)]
feats_all <- setdiff(intersect(manifest$feature_id, names(cases)), exclude_fids)
cat(sprintf("[2] testing %d features\n", length(feats_all)))

# ---- [3] cluster-vs-rest tests ----------------------------------------------
cat("[3] cluster-vs-rest tests\n")
cluster_levels <- levels(cases$cluster); results <- list()
for (cl in cluster_levels) {
  in_cl <- cases$cluster == cl; res_cl <- vector("list", length(feats_all))
  for (i in seq_along(feats_all)) {
    fid <- feats_all[i]; type <- manifest[feature_id == fid, ftype][1]; x <- cases[[fid]]
    if (all(is.na(x))) next
    if (type == "continuous") {
      x1 <- x[in_cl]; x0 <- x[!in_cl]
      if (sum(!is.na(x1)) < 5 || sum(!is.na(x0)) < 5) next
      pv <- tryCatch(wilcox.test(x1, x0, exact = FALSE)$p.value, error = function(e) NA_real_)
      res_cl[[i]] <- data.table(cluster = cl, feature_id = fid, type = "continuous",
        median_in = median(x1, na.rm = TRUE), median_out = median(x0, na.rm = TRUE),
        delta = median(x1, na.rm = TRUE) - median(x0, na.rm = TRUE),
        prev_in = NA_real_, prev_out = NA_real_, OR = NA_real_, p = pv)
    } else {
      x <- as.integer(x); x1 <- x[in_cl]; x0 <- x[!in_cl]
      tab <- table(factor(c(rep(1, sum(in_cl)), rep(0, sum(!in_cl))), levels = c(0, 1)), factor(x, levels = c(0, 1)))
      if (any(rowSums(tab) == 0) || any(colSums(tab) == 0)) next
      ft <- tryCatch(fisher.test(tab), error = function(e) NULL); if (is.null(ft)) next
      res_cl[[i]] <- data.table(cluster = cl, feature_id = fid, type = "binary",
        median_in = NA_real_, median_out = NA_real_, delta = NA_real_,
        prev_in = mean(x1, na.rm = TRUE), prev_out = mean(x0, na.rm = TRUE),
        OR = unname(ft$estimate), p = ft$p.value)
    }
  }
  results[[cl]] <- rbindlist(res_cl, fill = TRUE)
}
prof <- rbindlist(results, fill = TRUE)
prof[, q := p.adjust(p, method = "BH"), by = .(cluster, type)]
prof <- prof[order(cluster, q)]
prof[type == "continuous", abs_effect := abs(delta)]
prof[type == "binary", abs_effect := abs(log2(pmax(OR, 1e-3)))]
fwrite(prof, file.path(RPT_DIR, sprintf("cluster_profiles%s.csv", SUF)))

# ---- [4] top discriminators -------------------------------------------------
cat("[4] top discriminators per cluster (q<0.05)\n")
top_per_cluster <- prof[q < 0.05][order(cluster, -abs_effect), head(.SD, 20), by = cluster]
top_per_cluster <- merge(top_per_cluster, manifest[, .(feature_id, source_analysis, description)],
                         by = "feature_id", all.x = TRUE)
setcolorder(top_per_cluster, c("cluster", "feature_id", "source_analysis", "type", "median_in", "median_out",
  "delta", "prev_in", "prev_out", "OR", "p", "q", "abs_effect", "description"))
top_per_cluster <- top_per_cluster[order(cluster, -abs_effect)]
fwrite(top_per_cluster, file.path(TBL_DIR, sprintf("table_cluster_topfeatures%s.csv", SUF)))

# ---- [5] heatmap of top features --------------------------------------------
top_union <- unique(prof[q < 0.05][order(-abs_effect), head(feature_id, 30), by = cluster]$V1)
top_union <- intersect(top_union, names(cases))
if (length(top_union) >= 5) {
  M <- as.matrix(cases[, top_union, with = FALSE]); storage.mode(M) <- "numeric"
  M <- scale(M); M[is.na(M)] <- 0
  agg <- t(sapply(cluster_levels, function(cl) colMeans(M[cases$cluster == cl, , drop = FALSE], na.rm = TRUE)))
  rownames(agg) <- paste0("cluster_", cluster_levels); agg <- t(agg)
  agg[agg > 1.5] <- 1.5; agg[agg < -1.5] <- -1.5
  png(file.path(FIG_DIR, sprintf("heatmap_top_features%s.png", SUF)),
      width = max(800, 60 * ncol(agg) + 600), height = max(600, 14 * nrow(agg) + 200), res = 150)
  pheatmap(agg, cluster_cols = FALSE, color = colorRampPalette(c("steelblue", "white", "firebrick"))(50),
           breaks = seq(-1.5, 1.5, length.out = 51), main = "Top discriminators (union top-30 per cluster, z-mean)",
           fontsize_row = 7)
  dev.off()
} else cat("  skipped heatmap (too few significant features)\n")

# ---- [6] summary ------------------------------------------------------------
cat("\n[6] significant features per (cluster, type):\n")
print(prof[q < 0.05, .N, by = .(cluster, type)])
cat("\n=== done ===\n")
