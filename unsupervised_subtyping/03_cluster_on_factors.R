#!/usr/bin/env Rscript
# =============================================================================
# 03_cluster_on_factors.R  (subtype clustering)
#
# Take the MOFA+ factor scores (02_mofa_fit) and run k-means, PAM and HDBSCAN over a
# sweep of k. Select the final k by joint optimisation of silhouette, bootstrap
# Jaccard stability (Hennig 2007), and a >=5% cluster-size floor. Active factors
# = total variance >= 6% across views (from factor_activity.csv).
#
# CLI: arg1 = ""/"M"/"F" (sex); arg2 = forced k (optional).
#
# Input:  subtyping/mofa/mofa_factors_g1*.rds + subtyping/reports/factor_activity*.csv
# Output: subtyping/clusters/cluster_assignments_g1*.rds + subtyping/reports/cluster_sweep + final_k_decision + plots
# =============================================================================

rm(list = ls()); gc()
suppressPackageStartupMessages({ library(data.table); library(cluster); library(dbscan); library(fpc); library(ggplot2) })

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR); set.seed(20260413)

ARGS    <- commandArgs(trailingOnly = TRUE)
SEX     <- if (length(ARGS) >= 1 && nzchar(ARGS[1])) toupper(ARGS[1]) else ""
FORCE_K <- if (length(ARGS) >= 2 && nzchar(ARGS[2])) as.integer(ARGS[2]) else NA_integer_
stopifnot(SEX %in% c("", "M", "F"))
SEX_SUF <- if (nzchar(SEX)) paste0("_", tolower(SEX)) else ""
SUF     <- paste0(SEX_SUF, if (!is.na(FORCE_K)) paste0("_k", FORCE_K) else "")

P7_DIR  <- "output/subtyping"
RPT_DIR <- file.path(P7_DIR, "reports")
FIG_DIR <- file.path(RPT_DIR, "fig")
dir.create(file.path(P7_DIR, "clusters"), showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
while (sink.number() > 0) sink()
sink(file.path(RPT_DIR, sprintf("cluster_log%s.txt", SUF)), split = TRUE)
on.exit(sink(), add = TRUE)

cat("=== cluster on MOFA factors ===\n")

# ---- [1] load factors, select active ----------------------------------------
Z <- readRDS(file.path(P7_DIR, "mofa", sprintf("mofa_factors_g1%s.rds", SEX_SUF)))
activity <- fread(file.path(RPT_DIR, sprintf("factor_activity%s.csv", SEX_SUF)))
cat(sprintf("sex: %s | force_k: %s\n", if (nzchar(SEX)) SEX else "none", if (!is.na(FORCE_K)) FORCE_K else "auto"))
VAR_THRESH <- 6
factors_keep <- activity[total_var >= VAR_THRESH, factor]
cat(sprintf("  factors %d -> kept (total_var >= %g%%): %d (%s)\n",
            ncol(Z), VAR_THRESH, length(factors_keep), paste(factors_keep, collapse = ", ")))
Z_act <- Z[, factors_keep, drop = FALSE]; stopifnot(all(is.finite(Z_act)))

# ---- [2] sweep k for k-means + PAM ------------------------------------------
cat("[2] sweeping k=2..8 (k-means, PAM)\n")
K_RANGE <- 2:8; N <- nrow(Z_act); D <- dist(Z_act)
sweep_rows <- list(); labels_store <- list()
for (k in K_RANGE) {
  km <- kmeans(Z_act, centers = k, nstart = 25, iter.max = 100)
  sweep_rows[[length(sweep_rows) + 1]] <- data.table(method = "kmeans", k = k,
    silhouette = mean(silhouette(km$cluster, D)[, "sil_width"]), CH = calinhara(Z_act, km$cluster),
    min_cluster_pct = min(table(km$cluster)) / N)
  labels_store[[sprintf("kmeans_k%d", k)]] <- km$cluster
  pm <- pam(Z_act, k = k, metric = "euclidean", pamonce = 5)
  sweep_rows[[length(sweep_rows) + 1]] <- data.table(method = "pam", k = k,
    silhouette = pm$silinfo$avg.width, CH = calinhara(Z_act, pm$clustering),
    min_cluster_pct = min(table(pm$clustering)) / N)
  labels_store[[sprintf("pam_k%d", k)]] <- pm$clustering
}

# ---- [3] HDBSCAN ------------------------------------------------------------
cat("[3] HDBSCAN sweep (minPts {30,60,100,200})\n")
for (mp in c(30, 60, 100, 200)) {
  hd <- hdbscan(Z_act, minPts = mp); cl <- hd$cluster
  n_clusters <- length(unique(cl[cl > 0]))
  if (n_clusters >= 2) {
    valid <- cl > 0
    sil_hd <- if (sum(valid) > 10) mean(silhouette(cl[valid], dist(Z_act[valid, ]))[, "sil_width"]) else NA_real_
    ch_hd <- calinhara(Z_act[valid, ], cl[valid]); sz_hd <- min(table(cl[valid])) / N
  } else { sil_hd <- NA_real_; ch_hd <- NA_real_; sz_hd <- NA_real_ }
  sweep_rows[[length(sweep_rows) + 1]] <- data.table(method = sprintf("hdbscan_mp%d", mp),
    k = n_clusters, silhouette = sil_hd, CH = ch_hd, min_cluster_pct = sz_hd)
  labels_store[[sprintf("hdbscan_mp%d", mp)]] <- cl
}
sweep <- rbindlist(sweep_rows, fill = TRUE)

# ---- [4] bootstrap Jaccard for top candidates -------------------------------
cat("[4] bootstrap stability (B=100, 80% subsample) for k-means/PAM candidates\n")
cands <- sweep[method %in% c("kmeans", "pam") & min_cluster_pct >= 0.05 & silhouette > 0]
boot_jaccard <- function(method, k, B = 100, frac = 0.8) {
  jacc <- numeric(B)
  base <- if (method == "kmeans") kmeans(Z_act, centers = k, nstart = 25, iter.max = 100)$cluster
          else pam(Z_act, k = k, pamonce = 5)$clustering
  for (b in seq_len(B)) {
    idx <- sample(N, size = floor(frac * N), replace = FALSE)
    sub_lab <- if (method == "kmeans") kmeans(Z_act[idx, ], centers = k, nstart = 10, iter.max = 100)$cluster
               else pam(Z_act[idx, ], k = k, pamonce = 5)$clustering
    base_sub <- base[idx]
    j_per <- sapply(seq_len(k), function(cl) { a <- base_sub == cl
      max(sapply(seq_len(k), function(cl2) { b2 <- sub_lab == cl2; sum(a & b2) / max(1, sum(a | b2)) })) })
    jacc[b] <- mean(j_per)
  }
  jacc
}
cands[, bootstrap_jaccard := NA_real_]
for (i in seq_len(nrow(cands))) {
  jv <- boot_jaccard(cands$method[i], cands$k[i], B = 100)
  cands$bootstrap_jaccard[i] <- mean(jv)
  cat(sprintf("  %s k=%d  Jaccard=%.3f (sd=%.3f)\n", cands$method[i], cands$k[i], mean(jv), sd(jv)))
}
sweep <- merge(sweep, cands[, .(method, k, bootstrap_jaccard)], by = c("method", "k"), all.x = TRUE)
fwrite(sweep, file.path(RPT_DIR, sprintf("cluster_sweep%s.csv", SUF)))

# ---- [5] pick final ---------------------------------------------------------
cat("[5] selecting final solution\n")
elig <- sweep[!is.na(bootstrap_jaccard) & bootstrap_jaccard >= 0.6 & min_cluster_pct >= 0.05]
if (nrow(elig) == 0) { cat("  no solution meets bootstrap>=0.6 + size>=5%; falling back to best silhouette\n")
  elig <- sweep[method %in% c("kmeans", "pam") & min_cluster_pct >= 0.05] }
elig[, score := scale(silhouette)[, 1] + scale(bootstrap_jaccard)[, 1]]
final <- elig[which.max(score)]
if (!is.na(FORCE_K)) {
  forced <- sweep[method == "kmeans" & k == FORCE_K]
  if (nrow(forced) == 0) stop(sprintf("no kmeans entry for k=%d", FORCE_K))
  cat(sprintf("  overriding with forced k=%d (kmeans)\n", FORCE_K)); final <- forced
}
cat("  final pick:\n"); print(final)
final_lab <- labels_store[[sprintf("%s_k%d", final$method, final$k)]]

# ---- [6] save ---------------------------------------------------------------
saveRDS(list(factors_used = factors_keep, Z_active = Z_act, sweep = sweep, all_labels = labels_store,
             final_method = final$method, final_k = final$k, final_labels = final_lab,
             final_silhouette = final$silhouette, final_jaccard = final$bootstrap_jaccard),
        file.path(P7_DIR, "clusters", sprintf("cluster_assignments_g1%s.rds", SUF)))
writeLines(c("Final clustering decision",
  sprintf("Active factors (n=%d): %s", length(factors_keep), paste(factors_keep, collapse = ", ")),
  sprintf("Final method: %s | k: %d | silhouette: %.3f | bootstrap Jaccard: %.3f | min cluster: %.1f%%",
          final$method, final$k, final$silhouette, final$bootstrap_jaccard, 100 * final$min_cluster_pct),
  "", "Cluster sizes:", paste(capture.output(print(table(final_lab))), collapse = "\n")),
  file.path(RPT_DIR, sprintf("final_k_decision%s.txt", SUF)))

# ---- [7] figures ------------------------------------------------------------
p_sil <- ggplot(sweep[method %in% c("kmeans", "pam")], aes(x = k, y = silhouette, color = method, shape = method)) +
  geom_line() + geom_point(size = 3) + scale_x_continuous(breaks = K_RANGE) +
  labs(title = "Silhouette sweep (G1 cases on MOFA factors)", x = "k", y = "Average silhouette width") + theme_minimal()
ggsave(file.path(FIG_DIR, sprintf("silhouette_sweep%s.png", SUF)), p_sil, width = 7, height = 4.5, dpi = 180)
if (nrow(cands) > 0) {
  p_jc <- ggplot(cands, aes(x = k, y = bootstrap_jaccard, color = method, shape = method)) +
    geom_hline(yintercept = 0.6, linetype = 2, color = "grey50") +
    geom_hline(yintercept = 0.75, linetype = 2, color = "grey20") +
    geom_line() + geom_point(size = 3) + scale_x_continuous(breaks = K_RANGE) + ylim(0, 1) +
    labs(title = "Bootstrap Jaccard stability (B=100, 80% subsample)", x = "k", y = "Mean Jaccard") + theme_minimal()
  ggsave(file.path(FIG_DIR, sprintf("bootstrap_jaccard%s.png", SUF)), p_jc, width = 7, height = 4.5, dpi = 180)
}
cat("\n=== done ===\n")
