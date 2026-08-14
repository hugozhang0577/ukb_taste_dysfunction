#' =============================================================================
#' Effective number of tests and de-redundancy of the NMR panel
#' =============================================================================
#'
#' Purpose:
#'   1. PCA over the 327 NMR measures to estimate the effective number of tests
#'   2. hierarchical clustering of the correlation matrix to group collinear measures
#'   3. one representative per cluster, chosen using the primary MWAS model
#'   4. significance re-evaluated against the effective number of tests
#'
#' Input:
#'   - the preprocessed discovery-cohort NMR matrix (4xIQR trim, log1p, z-score)
#'   - the per-metabolite results of the primary MWAS model
#'
#' Output:
#'   - ent_report.txt            - effective-number-of-tests report
#'   - correlation_heatmap.pdf   - correlation matrix
#'   - pca_scree_plot.pdf        - PCA scree plot
#'   - cluster_dendrogram.pdf    - clustering dendrogram
#'   - metabolite_clusters.csv   - cluster membership per measure
#'   - deredundant_hits.csv      - de-redundant significant measures (main output)
#'
#' Requires: data.table, ggplot2, pheatmap, RColorBrewer
#' =============================================================================
# 0. Configuration ----
# =============================================================================

# --- paths ---
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)

NMR_FILE    <- file.path("input", "analysis_ready", "metabolomics_group1.csv")

# Per-metabolite results of the primary MWAS model
MWAS_RESULT <- file.path("output", "mwas", "primary", "primary_20cov.csv")

# Output directory
OUT_DIR     <- file.path("output", "mwas", "effective_tests")


# --- parameters ---
CORR_CUTOFF  <- 0.80    # clustering threshold; |r| > 0.80 puts two measures in the same cluster
VAR_95       <- 0.95     # components explaining 95% of the variance; one of the ENT estimates
VAR_99       <- 0.99     # components explaining 99% of the variance
FDR_ALPHA    <- 0.05     # FDR level
BONF_ALPHA   <- 0.05     # Bonferroni level

# --- Nightingale 37 CE-IVD Certified Biomarkers ---
# Source: Julkunen et al. (2021) eLife 10:e63033, Figure 2 source data
#   "37 biomarkers in the panel have been certified for diagnostics use"
# https://elifesciences.org/articles/63033
# See also Julkunen et al. (2023) Nat Commun 14:604, Fig 2a
# Names match the column names of the preprocessed NMR matrix.
CE_BIOMARKERS <- c(
  # --- lipoprotein lipids (5) ---
  "Total_C",            # total cholesterol
  "VLDL_C",             # VLDL cholesterol
  "LDL_C",              # LDL cholesterol
  "HDL_C",              # HDL cholesterol
  "Total_TG",           # total triglycerides
  
  # --- apolipoproteins (3) ---
  "ApoB",               # apolipoprotein B
  "ApoA1",              # apolipoprotein A-I
  "ApoB_by_ApoA1",      # ApoB / ApoA1 ratio
  
  # --- fatty acids, absolute concentration (7) ---
  "Total_FA",           # total fatty acids
  "Omega_3",            # omega-3 fatty acids
  "Omega_6",            # omega-6 fatty acids
  "PUFA",               # polyunsaturated fatty acids
  "MUFA",               # monounsaturated fatty acids
  "SFA",                # saturated fatty acids
  "DHA",                # docosahexaenoic acid
  
  # --- fatty acid percentages and ratios (8) ---
  "Omega_3_pct",        # omega-3 %
  "Omega_6_pct",        # omega-6 %
  "PUFA_pct",           # PUFA %
  "MUFA_pct",           # MUFA %
  "SFA_pct",            # SFA %
  "DHA_pct",            # DHA %
  "PUFA_by_MUFA",       # PUFA / MUFA
  "Omega_6_by_Omega_3", # omega-6 / omega-3
  
  # --- amino acids (9, including total BCAA) ---
  "Ala",                # alanine
  "Gly",                # glycine
  "His",                # histidine
  "Ile",                # isoleucine (BCAA)
  "Leu",                # leucine (BCAA)
  "Val",                # valine (BCAA)
  "Phe",                # phenylalanine
  "Tyr",                # tyrosine
  "Total_BCAA",         # total BCAA (Ile + Leu + Val)
  # Total_BCAA may be absent or differently named in some releases
  
  # --- glycolysis (2) ---
  "Glucose",            # glucose
  "Lactate",            # lactate
  
  # --- fluid balance and inflammation (3) ---
  "Creatinine",         # creatinine
  "Albumin",            # albumin
  "GlycA"               # glycoprotein acetyls (inflammation)
)
# 37 CE-IVD certified biomarkers in total.
#
# The following are biologically relevant but are not among the certified 37:
#   Remnant_C, Clinical_LDL_C, LA, LA_pct, Citrate, Pyruvate,
#   bOHbutyrate, Acetate, Acetoacetate, Acetone, Cholines,
#   sphingomyelins, glutamine (Gln)
# They are research-use grade and get no CE priority in the de-redundancy step.

# --- session ---
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("\n")
cat("======================================================================\n")
cat("  NMR panel: PCA and effective number of tests\n")
cat("======================================================================\n")
cat("  NMR data: ", NMR_FILE, "\n")
cat("  MWAS primary: ", MWAS_RESULT, "\n")
cat("  Output:   ", OUT_DIR, "\n")
cat("  Corr cutoff: |r| >", CORR_CUTOFF, "\n")
cat("======================================================================\n\n")


# 1. Load and check the NMR data ----
# =============================================================================

cat("[1] loading NMR data ...\n")
nmr <- fread(NMR_FILE)
cat("    dimensions:", nrow(nmr), "samples x", ncol(nmr), "columns\n")

# Identify the measure columns (everything but eid)
eid_col <- "eid"
metab_cols <- setdiff(names(nmr), eid_col)
cat("    measures:", length(metab_cols), "\n")

# Print the first ten names so the panel can be confirmed as Nightingale NMR
cat("    first ten measures:\n")
cat("      ", paste(head(metab_cols, 10), collapse = ", "), "\n")

# --- CE-IVD matching ---
ce_matched <- intersect(metab_cols, CE_BIOMARKERS)
ce_missing <- setdiff(CE_BIOMARKERS, metab_cols)

cat("\n    CE-IVD certified biomarkers (Julkunen et al. 2021 eLife):\n")
cat("      defined:", length(CE_BIOMARKERS), "\n")
cat("      matched: ", length(ce_matched), "\n")
if (length(ce_missing) > 0) {
  cat("      unmatched:", length(ce_missing), "\n")
  cat("        ", paste(ce_missing, collapse = ", "), "\n")
  cat("      [CHECK] the column name may differ, or the measure may be absent\n")
}
cat("      non-CE columns: ", length(metab_cols) - length(ce_matched), 
    " (lipoprotein subclass detail and similar)\n")

# Flag vector: TRUE = CE-IVD certified
is_ce <- metab_cols %in% ce_matched
names(is_ce) <- metab_cols

# Numeric matrix
nmr_mat <- as.matrix(nmr[, ..metab_cols])
cat("    numeric matrix:", nrow(nmr_mat), "×", ncol(nmr_mat), "\n")

# Missingness
n_na <- sum(is.na(nmr_mat))
cat("    missing:", n_na)
if (n_na > 0) {
  pct_na <- round(n_na / length(nmr_mat) * 100, 3)
  cat(" (", pct_na, "% - imputed with the column median so that PCA can run)\n")
  for (j in 1:ncol(nmr_mat)) {
    na_idx <- is.na(nmr_mat[, j])
    if (any(na_idx)) {
      nmr_mat[na_idx, j] <- median(nmr_mat[, j], na.rm = TRUE)
    }
  }
} else {
  cat(" (clean)\n")
}

# Check the standardisation state
# The z-score was computed on the pooled baseline sample, so within the
# discovery cohort the mean/SD drift from 0/1. That is expected and does not
col_means <- round(colMeans(nmr_mat), 3)
col_sds   <- round(apply(nmr_mat, 2, sd), 3)
cat("    column mean range: [", min(col_means), ",", max(col_means), "]\n")
cat("    column SD range: [", min(col_sds), ",", max(col_sds), "]\n")
if (max(abs(col_means)) > 1.0 || min(col_sds) < 0.3 || max(col_sds) > 3.0) {
  cat("    [WARNING] the data may not be standardised; check the preprocessing step\n")
} else if (max(abs(col_means)) > 0.2 || max(abs(col_sds - 1)) > 0.2) {
  cat("    [OK] drift is within the range expected for a subset of the pooled sample\n")
  cat("         PCA uses center = TRUE, scale. = TRUE and is unaffected\n")
} else {
  cat("    [OK] data are standardised (mean ~0, SD ~1)\n")
}

cat("\n")


# 2. Correlation matrix ----
# =============================================================================

cat("[2] computing the correlation matrix ...\n")
cor_mat <- cor(nmr_mat, use = "pairwise.complete.obs")
cat("    correlation matrix:", nrow(cor_mat), "×", ncol(cor_mat), "\n")

# Distribution of the correlations
upper_tri <- cor_mat[upper.tri(cor_mat)]
cat("    upper-triangle elements:", length(upper_tri), "\n")
cat("    |r| distribution:\n")
cat("      Mean |r|:   ", round(mean(abs(upper_tri)), 3), "\n")
cat("      Median |r|: ", round(median(abs(upper_tri)), 3), "\n")
cat("      |r| > 0.5:  ", sum(abs(upper_tri) > 0.5), " (",
    round(mean(abs(upper_tri) > 0.5) * 100, 1), "%)\n")
cat("      |r| > 0.8:  ", sum(abs(upper_tri) > 0.8), " (",
    round(mean(abs(upper_tri) > 0.8) * 100, 1), "%)\n")
cat("      |r| > 0.9:  ", sum(abs(upper_tri) > 0.9), " (",
    round(mean(abs(upper_tri) > 0.9) * 100, 1), "%)\n")
cat("      |r| > 0.95: ", sum(abs(upper_tri) > 0.95), " (",
    round(mean(abs(upper_tri) > 0.95) * 100, 1), "%)\n")

# If more than ~30% of pairs exceed |r| > 0.8, confirm that is expected
cat("\n")


# 3. PCA ----
# =============================================================================

cat("[3] running PCA ...\n")
pca_result <- prcomp(nmr_mat, center = TRUE, scale. = TRUE)

# Variance explained
eigenvalues <- pca_result$sdev^2
var_prop    <- eigenvalues / sum(eigenvalues)
var_cumsum  <- cumsum(var_prop)

# Key cumulative-variance points
n_pc_95  <- which(var_cumsum >= VAR_95)[1]
n_pc_99  <- which(var_cumsum >= VAR_99)[1]
n_pc_50  <- which(var_cumsum >= 0.50)[1]
n_pc_80  <- which(var_cumsum >= 0.80)[1]

cat("    components:", length(eigenvalues), "\n")
cat("    variance explained by PC1:", round(var_prop[1] * 100, 1), "%\n")
cat("    cumulative PC1-5:  ", round(var_cumsum[5] * 100, 1), "%\n")
cat("    cumulative PC1-10: ", round(var_cumsum[10] * 100, 1), "%\n")
cat("\n")
cat("    >>> cumulative variance <<<\n")
cat("    50% of variance: PC1-", n_pc_50, "\n")
cat("    80% of variance: PC1-", n_pc_80, "\n")
cat("    95% of variance: PC1-", n_pc_95, " ← ENT_PCA_95\n")
cat("    99% of variance: PC1-", n_pc_99, " ← ENT_PCA_99\n")

# --- Scree plot ---
pdf(file.path(OUT_DIR, "scree_plot.pdf"), width = 10, height = 6)

par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))

# left: scree (first 60)
n_show <- min(60, length(eigenvalues))
barplot(var_prop[1:n_show] * 100, 
        names.arg = 1:n_show,
        xlab = "Principal Component", 
        ylab = "Variance Explained (%)",
        main = "Scree Plot (NMR Metabolites)",
        col = ifelse(1:n_show <= n_pc_95, "#4292c6", "#d9d9d9"),
        border = NA, cex.names = 0.6)
abline(v = n_pc_95, col = "red", lty = 2, lwd = 2)
text(n_pc_95 + 1, max(var_prop[1:n_show] * 100) * 0.9, 
     paste0("PC", n_pc_95, "\n(95% var)"), col = "red", cex = 0.8, adj = 0)

# right: cumulative variance
plot(1:n_show, var_cumsum[1:n_show] * 100, type = "l", lwd = 2,
     xlab = "Number of PCs", ylab = "Cumulative Variance (%)",
     main = "Cumulative Variance Explained", col = "#2171b5")
abline(h = 95, col = "red", lty = 2)
abline(h = 99, col = "orange", lty = 2)
abline(v = n_pc_95, col = "red", lty = 3)
abline(v = n_pc_99, col = "orange", lty = 3)
points(n_pc_95, 95, pch = 19, col = "red", cex = 1.5)
points(n_pc_99, 99, pch = 19, col = "orange", cex = 1.5)
legend("bottomright", 
       legend = c(paste0("95% → ", n_pc_95, " PCs"),
                  paste0("99% → ", n_pc_99, " PCs")),
       col = c("red", "orange"), lty = 2, cex = 0.9)

dev.off()
cat("    Scree plot saved\n")


# 4. Effective number of tests: three estimators ----
# =============================================================================

cat("\n[4] effective number of tests ...\n")
cat("    three published estimators:\n\n")

# --- 4a. Nyholt (2004), eigenvalue-based ---
# M_eff = 1 + (M - 1) * (1 - Var(lambda_obs) / M)
# lambda_obs are the eigenvalues of the correlation matrix
eigen_cor   <- eigen(cor_mat, only.values = TRUE)$values
M           <- length(eigen_cor)
var_lambda  <- var(eigen_cor)  # standardised eigenvalues have mean 1
ent_nyholt  <- 1 + (M - 1) * (1 - var_lambda / M)
ent_nyholt  <- round(ent_nyholt)

cat("    [4a] Nyholt (2004):\n")
cat("         M_eff = 1 + (M-1) * (1 - Var(eigenvalues)/M)\n")
cat("         Var(eigenvalues) =", round(var_lambda, 3), "\n")
cat("         ENT_Nyholt =", ent_nyholt, "\n\n")

# --- 4b. Li & Ji (2005), refinement of Nyholt ---
# per eigenvalue x: f(x) = I(x >= 1) + (x - floor(x)) if x > 1
#                         x                               if x <= 1
# M_eff_LJ = sum(f(x))
f_lj <- function(x) {
  ifelse(x >= 1, 1 + (x - floor(x)), x)
}
ent_liji <- round(sum(f_lj(eigen_cor)))

cat("    [4b] Li & Ji (2005):\n")
cat("         M_eff = sum f(lambda_i), f(x) = I(x >= 1) + (x - floor(x))\n")
cat("         ENT_LiJi =", ent_liji, "\n\n")

# --- 4c. Galwey (2009) / PCA-based ---
# ENT = number of PCs needed to explain a given share of the variance
cat("    [4c] PCA-based (Galwey 2009):\n")
cat("         ENT_PCA_95 =", n_pc_95, " (95% of variance)\n")
cat("         ENT_PCA_99 =", n_pc_99, " (99% of variance)\n\n")

# --- summary ---
ent_summary <- data.table(
  Method = c("Nyholt (2004)", "Li & Ji (2005)", 
             "PCA 95%", "PCA 99%",
             "Raw (no correction)"),
  ENT = c(ent_nyholt, ent_liji, n_pc_95, n_pc_99, M),
  Bonf_Threshold = format(BONF_ALPHA / c(ent_nyholt, ent_liji, n_pc_95, n_pc_99, M),
                          scientific = TRUE, digits = 3),
  Description = c(
    "Eigenvalue variance-based",
    "Improved eigenvalue truncation",
    "PCs for 95% variance",
    "PCs for 99% variance",
    "Naive (assumes independence)"
  )
)

cat("    ┌────────────────────────────────────────────────────────────┐\n")
cat("    │              ENT ESTIMATION SUMMARY                       │\n")
cat("    ├──────────────────┬──────┬──────────────┬──────────────────┤\n")
cat(sprintf("    │ %-16s │ %4s │ %12s │ %-16s │\n",
            "Method", "ENT", "Bonf Thresh", "Note"))
cat("    ├──────────────────┼──────┼──────────────┼──────────────────┤\n")
for (i in 1:nrow(ent_summary)) {
  cat(sprintf("    │ %-16s │ %4d │ %12s │ %-16s │\n",
              ent_summary$Method[i], ent_summary$ENT[i],
              ent_summary$Bonf_Threshold[i],
              substr(ent_summary$Description[i], 1, 16)))
}
cat("    └──────────────────┴──────┴──────────────┴──────────────────┘\n")

# Preferred estimator
ent_recommended <- ent_liji  # Li & Ji is the most widely used estimator in the genetic-epidemiology literature
cat("\n    >>> reported estimator: Li & Ji ENT =", ent_recommended, "<<<\n")
cat("    It is the most widely cited estimator for multi-marker association\n")
cat("    studies and balances conservatism against power\n\n")


# 5. Hierarchical clustering (de-redundancy) ----
# =============================================================================

cat("[5] clustering (|r| >", CORR_CUTOFF, ") ...\n")

# Distance = 1 - |r|
dist_mat <- as.dist(1 - abs(cor_mat))
hc <- hclust(dist_mat, method = "complete")

# Cut at the threshold
cluster_ids <- cutree(hc, h = 1 - CORR_CUTOFF)
n_clusters  <- max(cluster_ids)

cat("    linkage: complete\n")
cat("    cut height: 1 - |r| =", 1 - CORR_CUTOFF, " (i.e. |r| >", CORR_CUTOFF, ")\n")
cat("    clusters: ", n_clusters, " (from ", M, " measures)\n")
cat("    compression: ", round(M / n_clusters, 1), "×\n")

# Cluster-size distribution
cluster_sizes <- table(cluster_ids)
cat("    cluster sizes:\n")
cat("      singletons:       ", sum(cluster_sizes == 1), "\n")
cat("      2-5 measures:     ", sum(cluster_sizes >= 2 & cluster_sizes <= 5), "\n")
cat("      6-10 measures:    ", sum(cluster_sizes >= 6 & cluster_sizes <= 10), "\n")
cat("      >10 measures:     ", sum(cluster_sizes > 10), "\n")

# Cluster membership table
cluster_dt <- data.table(
  metabolite   = metab_cols,
  cluster_id   = cluster_ids[metab_cols],
  cluster_size = as.integer(cluster_sizes[as.character(cluster_ids[metab_cols])]),
  is_ce        = is_ce[metab_cols]   # CE-IVD flag
)

# CE membership per cluster
ce_per_cluster <- cluster_dt[, .(n_ce = sum(is_ce), n_total = .N), by = cluster_id]
cat("    clusters containing a CE measure:", sum(ce_per_cluster$n_ce > 0), "/", n_clusters, "\n")
cat("    clusters with no CE measure:     ", sum(ce_per_cluster$n_ce == 0), "/", n_clusters, "\n")

# --- dendrogram ---
pdf(file.path(OUT_DIR, "cluster_dendrogram.pdf"), width = 20, height = 8)
par(mar = c(2, 4, 3, 1))
plot(hc, labels = FALSE, hang = -1,
     main = paste0("NMR Metabolite Clustering (", M, " → ", n_clusters, " clusters, |r| > ", CORR_CUTOFF, ")"),
     xlab = "", ylab = "1 - |r|", cex.main = 1.2)
abline(h = 1 - CORR_CUTOFF, col = "red", lty = 2, lwd = 2)
text(M * 0.85, 1 - CORR_CUTOFF + 0.02,
     paste0("cutoff: |r| = ", CORR_CUTOFF, " → ", n_clusters, " clusters"),
     col = "red", cex = 0.9)
dev.off()
cat("    Dendrogram saved\n")


# 6. Correlation heatmap, ordered by cluster ----
# =============================================================================

cat("\n[6] correlation heatmap ...\n")

# Order by cluster
order_idx <- order(cluster_ids)
cor_ordered <- cor_mat[order_idx, order_idx]

# Cluster annotation
annotation_row <- data.frame(
  Cluster = factor(cluster_ids[order_idx]),
  row.names = metab_cols[order_idx]
)

# Colour the larger clusters; small ones stay grey
big_clusters <- names(cluster_sizes[cluster_sizes >= 5])
n_big <- length(big_clusters)
if (n_big > 0) {
  palette_big <- colorRampPalette(brewer.pal(min(n_big, 12), "Set3"))(n_big)
  names(palette_big) <- big_clusters
}

pdf(file.path(OUT_DIR, "correlation_heatmap.pdf"), width = 14, height = 12)
pheatmap(cor_ordered,
         color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100),
         breaks = seq(-1, 1, length.out = 101),
         cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames = FALSE, show_colnames = FALSE,
         main = paste0("NMR Metabolite Correlation Matrix (n=", M,
                       ", ordered by ", n_clusters, " clusters)"),
         annotation_row = annotation_row,
         annotation_names_row = FALSE,
         legend_breaks = c(-1, -0.5, 0, 0.5, 1))
dev.off()
cat("    Heatmap saved\n")


# 7. Pick one representative per cluster using the MWAS results ----
# =============================================================================

cat("\n[7] reading the primary MWAS results ...\n")

if (!file.exists(MWAS_RESULT)) {
  cat("    [WARNING] MWAS result file not found:", MWAS_RESULT, "\n")
  cat("    falling back on PC1 loading to pick representatives\n")
  use_mwas <- FALSE
} else {
  mwas_dt <- fread(MWAS_RESULT)
  cat("    MWAS results:", nrow(mwas_dt), " measures\n")
  use_mwas <- TRUE
}

if (use_mwas) {
  
  # Join cluster membership to the MWAS results
  # The engine writes the measure name in the "protein" column
  if (!"protein" %in% names(mwas_dt)) {
    cat("    [CHECK] MWAS result columns:", paste(head(names(mwas_dt)), collapse = ", "), "\n")
    stop("Confirm which column carries the measure name")
  }
  
  merged <- merge(cluster_dt, mwas_dt, 
                  by.x = "metabolite", by.y = "protein", 
                  all.x = TRUE)
  
  cat("    matched:", sum(!is.na(merged$pval)), "/", nrow(merged), "\n")
  
  # =========================================================================
  # Representative selection: CE-IVD priority
  # =========================================================================
  # Rule:
  #   1. if the cluster contains a CE measure, take the CE measure with the smallest P
  #   2. otherwise take the cluster member with the smallest P
  # Rationale:
  #   - CE measures have defined clinical assays and reference ranges
  #   - downstream model explanations then map onto quantities clinicians use
  #   - they are comparable across studies (Nightingale IVD certification)
  #   - within a cluster |r| > 0.8, so little signal is lost by preferring the CE member
  # =========================================================================
  
  select_representative <- function(dt) {
    # dt: all measures of one cluster, carrying is_ce and pval
    dt_valid <- dt[!is.na(pval)]
    if (nrow(dt_valid) == 0) return(dt[1])  # All P values missing; take the first
    
    ce_subset <- dt_valid[is_ce == TRUE]
    if (nrow(ce_subset) > 0) {
      # Rule 1: smallest P among the CE members
      return(ce_subset[which.min(pval)])
    } else {
      # Rule 2: no CE member; smallest P in the cluster
      return(dt_valid[which.min(pval)])
    }
  }
  
  reps <- merged[, select_representative(.SD), by = cluster_id]
  reps <- reps[order(pval)]
  
  # Report the selection
  n_ce_reps     <- sum(reps$is_ce, na.rm = TRUE)
  n_non_ce_reps <- sum(!reps$is_ce, na.rm = TRUE)
  
  cat("\n    representatives:", nrow(reps), "\n")
  cat("    ┌────────────────────────────────────────────────────────┐\n")
  cat("    │  CE-PRIORITY REPRESENTATIVE SELECTION                  │\n")
  cat("    ├────────────────────────────────────────────────────────┤\n")
  cat(sprintf("    |  CE representative:     %3d (cluster had a CE member)   |\n", n_ce_reps))
  cat(sprintf("    |  non-CE representative: %3d (no CE member; min P)       |\n", n_non_ce_reps))
  cat(sprintf("    |  CE coverage:           %.1f%%                            |\n",
              n_ce_reps / nrow(reps) * 100))
  cat("    └────────────────────────────────────────────────────────┘\n")
  
  # How often does CE priority override the plain minimum-P choice?
  reps_blind <- merged[!is.na(pval), .SD[which.min(pval)], by = cluster_id]
  n_override <- sum(reps$metabolite != reps_blind$metabolite[match(reps$cluster_id, reps_blind$cluster_id)],
                    na.rm = TRUE)
  
  if (n_override > 0) {
    cat(sprintf("    [INFO] CE priority overrode the minimum-P choice in %d cluster(s)\n", n_override))
    
    # Show them, so the P-value cost of preferring CE can be judged
    override_clusters <- reps$cluster_id[reps$metabolite != 
                                           reps_blind$metabolite[match(reps$cluster_id, reps_blind$cluster_id)]]
    override_clusters <- override_clusters[!is.na(override_clusters)]
    
    if (length(override_clusters) > 0) {
      cat("    overridden clusters (first ten):\n")
      cat(sprintf("    %-6s  %-25s %-12s  %-25s %-12s  %-8s\n",
                  "Cluster", "CE-selected", "CE P-val", "Blind min-P", "Blind P-val", "P ratio"))
      
      n_show <- min(10, length(override_clusters))
      for (k in seq_len(n_show)) {
        cid   <- override_clusters[k]
        ce_row <- reps[cluster_id == cid]
        bl_row <- reps_blind[cluster_id == cid]
        p_ratio <- ce_row$pval / bl_row$pval
        cat(sprintf("    %-6d  %-25s %-12s  %-25s %-12s  %.2f×\n",
                    cid,
                    substr(ce_row$metabolite, 1, 25),
                    format(ce_row$pval, digits = 2, scientific = TRUE),
                    substr(bl_row$metabolite, 1, 25),
                    format(bl_row$pval, digits = 2, scientific = TRUE),
                    p_ratio))
      }
      cat("    P ratio below 2x means little is lost by preferring the CE measure;\n")
      cat("    above 5x the cluster is worth reviewing\n")
    }
  } else {
    cat("    [OK] CE priority never overrode minimum-P (the CE measure was always strongest)\n")
  }
  
} else {
  
  # Fallback: the measure with the largest PC1 loading, CE preferred
  loadings_pc1 <- abs(pca_result$rotation[, 1])
  cluster_dt[, pc1_loading := loadings_pc1[metabolite]]
  
  select_rep_loading <- function(dt) {
    ce_subset <- dt[is_ce == TRUE]
    if (nrow(ce_subset) > 0) return(ce_subset[which.max(pc1_loading)])
    return(dt[which.max(pc1_loading)])
  }
  
  reps <- cluster_dt[, select_rep_loading(.SD), by = cluster_id]
  reps[, pval := NA_real_]
  reps <- reps[order(-pc1_loading)]
  cat("    representatives (PC1 loading, CE preferred):", nrow(reps), "\n")
}


# 8. Significance under the effective number of tests ----
# =============================================================================

cat("\n[8] re-evaluating significance ...\n")

if (use_mwas) {
  
  # --- A: Bonferroni over the full panel, using ENT as the test count ---
  bonf_ent  <- BONF_ALPHA / ent_recommended
  bonf_raw  <- BONF_ALPHA / M
  
  mwas_dt[, sig_bonf_raw := pval < bonf_raw]
  mwas_dt[, sig_bonf_ent := pval < bonf_ent]
  
  # --- B: FDR over the cluster representatives ---
  reps[, pval_fdr_ent := p.adjust(pval, method = "BH")]
  reps[, pval_bonf_ent := p.adjust(pval, method = "bonferroni")]
  reps[, sig_fdr_ent := pval_fdr_ent < FDR_ALPHA]
  reps[, sig_bonf_ent := pval_bonf_ent < BONF_ALPHA]
  
  n_sig_fdr_ent  <- sum(reps$sig_fdr_ent, na.rm = TRUE)
  n_sig_bonf_ent <- sum(reps$sig_bonf_ent, na.rm = TRUE)
  
  cat("\n    ┌─────────────────────────────────────────────────────────┐\n")
  cat("    │           SIGNIFICANCE RE-ASSESSMENT                    │\n")
  cat("    ├─────────────────────────────────────────────────────────┤\n")
  cat(sprintf("    │  Raw 327 markers, naive FDR:       %3d FDR-sig       │\n",
              sum(mwas_dt$pval_fdr < 0.05, na.rm = TRUE)))
  cat(sprintf("    │  Raw 327 markers, naive Bonf:      %3d Bonf-sig      │\n",
              sum(mwas_dt$pval_bonf < 0.05, na.rm = TRUE)))
  cat("    │                                                         │\n")
  cat(sprintf("    │  Raw 327 markers, ENT-Bonf:        %3d ENT-Bonf-sig  │\n",
              sum(mwas_dt$sig_bonf_ent, na.rm = TRUE)))
  cat("    │                                                         │\n")
  cat(sprintf("    │  De-redundant %3d reps, FDR:       %3d FDR-sig       │\n",
              nrow(reps), n_sig_fdr_ent))
  cat(sprintf("    │  De-redundant %3d reps, Bonf:      %3d Bonf-sig      │\n",
              nrow(reps), n_sig_bonf_ent))
  cat("    └─────────────────────────────────────────────────────────┘\n")
  
  cat("\n    reported: FDR over the cluster representatives\n")
  cat("    one representative per cluster removes the collinearity that violates\n")
  cat("    the independence assumption behind FDR\n\n")
  
} else {
  cat("    [SKIP] no MWAS results; significance not re-evaluated\n")
  cat("    run the MWAS scan first, then re-run this script\n\n")
}


# 9. Genomic-control lambda before and after de-redundancy ----
# =============================================================================

cat("[9] lambda after de-redundancy ...\n")

if (use_mwas) {
  
  # Lambda over the full panel
  valid_p_all <- mwas_dt$pval[!is.na(mwas_dt$pval)]
  chi2_all    <- qchisq(1 - valid_p_all, df = 1)
  lambda_raw  <- median(chi2_all) / qchisq(0.5, df = 1)
  
  # Lambda over the representatives
  valid_p_reps <- reps$pval[!is.na(reps$pval)]
  chi2_reps    <- qchisq(1 - valid_p_reps, df = 1)
  lambda_ent   <- median(chi2_reps) / qchisq(0.5, df = 1)
  
  cat(sprintf("    lambda (full panel, %d measures):   %.3f\n", M, lambda_raw))
  cat(sprintf("    lambda (%d representatives):        %.3f\n", nrow(reps), lambda_ent))
  cat(sprintf("    reduction:                          %.1f%%\n",
              (1 - lambda_ent / lambda_raw) * 100))
  
  if (lambda_ent > 2) {
    cat("    [NOTE] lambda is still above 2 after de-redundancy\n")
    cat("           possible reasons: (1) genuinely widespread biological signal\n")
    cat("                             (2) residual collinearity\n")
    cat("                             (3) residual confounding\n")
    cat("           report both values and discuss\n")
  } else if (lambda_ent > 1.5) {
    cat("    [OK] lambda is between 1 and 2 after de-redundancy\n")
    cat("         the same order as the primary PWAS model (lambda 1.337)\n")
  } else {
    cat("    [OK] lambda is close to 1\n")
  }
  
} else {
  cat("    [SKIP] no MWAS results\n")
}

cat("\n")


# 10. Write outputs ----
# =============================================================================

cat("[10] writing outputs ...\n")

# --- 10a. cluster membership ---
if (use_mwas) {
  cluster_output <- merge(cluster_dt, 
                          mwas_dt[, .(protein, beta, se, or, pval, pval_fdr)],
                          by.x = "metabolite", by.y = "protein", all.x = TRUE)
  cluster_output[, is_representative := metabolite %in% reps$metabolite]
  setorder(cluster_output, cluster_id, pval)
} else {
  cluster_output <- cluster_dt
}
fwrite(cluster_output, file.path(OUT_DIR, "metabolite_clusters.csv"))
cat("    [OK] metabolite_clusters.csv\n")

# --- 10b. de-redundant significant list (main output) ---
if (use_mwas) {
  # Key columns
  output_cols <- intersect(
    c("cluster_id", "cluster_size", "metabolite", "is_ce", "beta", "se", "or", 
      "pval", "pval_fdr_ent", "pval_bonf_ent", "sig_fdr_ent", "sig_bonf_ent",
      "n_total", "n_case", "or_lower", "or_upper"),
    names(reps)
  )
  reps_output <- reps[, ..output_cols]
  setorder(reps_output, pval)
  fwrite(reps_output, file.path(OUT_DIR, "deredundant_hits.csv"))
  cat("    [OK] deredundant_hits.csv (", nrow(reps_output), " representatives)\n")
  
  # Significant subset
  if (n_sig_fdr_ent > 0) {
    sig_output <- reps_output[sig_fdr_ent == TRUE]
    fwrite(sig_output, file.path(OUT_DIR, "significant_metabolites.csv"))
    cat("    [OK] significant_metabolites.csv (", nrow(sig_output), " FDR-significant)\n")
  }
}

# --- 10c. ENT report ---
report_file <- file.path(OUT_DIR, "ent_report.txt")
sink(report_file)

cat("======================================================================\n")
cat("  NMR panel: effective tests and de-redundancy report\n")
cat("  Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("======================================================================\n\n")

cat("[DATA]\n")
cat("  NMR file:", NMR_FILE, "\n")
cat("  Samples:", nrow(nmr), "\n")
cat("  Metabolites:", M, "(CE-IVD:", length(ce_matched), ")\n\n")

cat("[CORRELATION STRUCTURE]\n")
cat("  Mean |r|:   ", round(mean(abs(upper_tri)), 3), "\n")
cat("  Median |r|: ", round(median(abs(upper_tri)), 3), "\n")
cat("  % pairs |r| > 0.5:", round(mean(abs(upper_tri) > 0.5) * 100, 1), "%\n")
cat("  % pairs |r| > 0.8:", round(mean(abs(upper_tri) > 0.8) * 100, 1), "%\n\n")

cat("[PCA]\n")
cat("  PC1 variance:", round(var_prop[1] * 100, 1), "%\n")
cat("  50% variance: PC1-", n_pc_50, "\n")
cat("  80% variance: PC1-", n_pc_80, "\n")
cat("  95% variance: PC1-", n_pc_95, "\n")
cat("  99% variance: PC1-", n_pc_99, "\n\n")

cat("[EFFECTIVE NUMBER OF TESTS]\n")
print(ent_summary)
cat("\n  Recommended: Li & Ji ENT =", ent_recommended, "\n\n")

cat("[CLUSTERING]\n")
cat("  Method: Complete linkage, cutoff |r| >", CORR_CUTOFF, "\n")
cat("  Clusters:", n_clusters, "\n")
cat("  Compression ratio:", round(M / n_clusters, 1), "×\n\n")

if (use_mwas) {
  cat("[CE-PRIORITY REPRESENTATIVE SELECTION]\n")
  cat("  CE representatives:     ", n_ce_reps, "\n")
  cat("  Non-CE representatives: ", n_non_ce_reps, "\n")
  cat("  CE coverage:            ", round(n_ce_reps / nrow(reps) * 100, 1), "%\n\n")
  
  cat("[SIGNIFICANCE AFTER CORRECTION]\n")
  cat("  Raw FDR-sig:                ", sum(mwas_dt$pval_fdr < 0.05, na.rm = TRUE), "\n")
  cat("  De-redundant FDR-sig:       ", n_sig_fdr_ent, "\n")
  cat("  De-redundant Bonf-sig:      ", n_sig_bonf_ent, "\n\n")
  
  cat("[LAMBDA]\n")
  cat("  Raw (327 markers):          ", round(lambda_raw, 3), "\n")
  cat("  De-redundant (", nrow(reps), " reps):    ", round(lambda_ent, 3), "\n\n")
  
  cat("[TOP 20 DE-REDUNDANT HITS]\n")
  top20 <- head(reps_output, 20)
  top20_cols <- intersect(c("metabolite", "is_ce", "cluster_id", "cluster_size", 
                            "or", "pval", "pval_fdr_ent"), names(top20))
  top20_print <- top20[, ..top20_cols]
  if ("or" %in% names(top20_print)) top20_print[, or := round(or, 3)]
  if ("pval_fdr_ent" %in% names(top20_print)) top20_print[, pval_fdr_ent := round(pval_fdr_ent, 4)]
  print(top20_print)
}

sink()
cat("    [OK] ent_report.txt\n")

# --- 10d. PCA variance table ---
pca_table <- data.table(
  PC = paste0("PC", 1:length(eigenvalues)),
  Eigenvalue = round(eigenvalues, 4),
  Var_Prop = round(var_prop * 100, 2),
  Var_Cumsum = round(var_cumsum * 100, 2)
)
fwrite(pca_table, file.path(OUT_DIR, "pca_variance.csv"))
# Export and dx upload to RAP  (metabolite_clusters.csv feeds 07; the remaining
# tables are the effective-test and de-redundancy diagnostics)
cat("    [OK] pca_variance.csv\n")

# --- 10e. correlation matrix (RDS) ---
saveRDS(cor_mat, file.path(OUT_DIR, "correlation_matrix.rds"))
cat("    [OK] correlation_matrix.rds\n")


# 11. Summary ----
# =============================================================================

cat("\n")
cat("======================================================================\n")
cat("  DE-REDUNDANCY COMPLETE\n")
cat("======================================================================\n")
cat("  measures in the panel:     ", M, "\n")
cat("    of which CE-IVD:         ", length(ce_matched), "\n")
cat("  effective number of tests: ", ent_recommended, " (Li & Ji)\n")
cat("  clusters:                  ", n_clusters, "\n")
if (use_mwas) {
  cat("  representatives:           ", nrow(reps), "\n")
  cat("    CE:                      ", n_ce_reps, " (",
      round(n_ce_reps / nrow(reps) * 100, 1), "%)\n")
  cat("    non-CE:                  ", n_non_ce_reps, "\n")
  cat("  FDR-significant:           ", n_sig_fdr_ent, "\n")
  cat("  Bonferroni-significant:    ", n_sig_bonf_ent, "\n")
  cat("  Lambda (raw → derep): ", round(lambda_raw, 3), " → ", round(lambda_ent, 3), "\n")
  
  # Share of the FDR-significant set that is CE certified
  if (n_sig_fdr_ent > 0) {
    sig_reps <- reps[sig_fdr_ent == TRUE]
    n_sig_ce <- sum(sig_reps$is_ce, na.rm = TRUE)
    cat("  CE among FDR-significant:  ", n_sig_ce, "/", n_sig_fdr_ent, " (",
        round(n_sig_ce / n_sig_fdr_ent * 100, 1), "%)\n")
  }
}
cat("\n")
cat("  [NEXT STEPS]\n")
cat("  1. review the significant measures in deredundant_hits.csv\n")
cat("     rows with is_ce = TRUE are the preferred candidates for the feature set\n")
cat("  2. confirm the drop in lambda is as expected\n")
cat("  3. if any representative is FDR-significant, continue to biological annotation\n")
cat("  4. if none is, consider relaxing CORR_CUTOFF to 0.7\n")
cat("  5. cross-omics check: which CE measure represents the HDL cluster that\n     the lead PWAS protein PON3 tracks\n")
cat("======================================================================\n")
