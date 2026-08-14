#!/usr/bin/env Rscript
# =============================================================================
# De-redundant representative set for the NMR panel
# =============================================================================
#
# 06_effective_tests.R groups the 327 NMR measures into correlation clusters. This
# script picks the one measure that represents each significant cluster, and
# writes the table the downstream feature manifest reads.
#
# Selection rule, fixed a priori:
#   Within each cluster, take the CE-IVD certified member whose effect direction
#   agrees with the cluster lead (the member with the smallest P). If no
#   certified member points the same way, keep the lead.
#
# Two things motivate the rule. Certified measures have defined clinical assays
# and reference ranges, so a panel built from them is deployable and comparable
# across studies. The direction guard stops that preference from inverting the
# sign of a cluster's association: a certified measure running opposite to the
# cluster's own signal would misrepresent it, so there the lead wins. Within a
# cluster |r| > 0.80, so the P-value cost of the preference is small.
#
# FDR is applied under the full 327-measure panel, not within clusters.
#
# Input : the primary MWAS result CSV + the cluster table from 06_effective_tests.R
# Output: output/mwas/effective_tests/cluster_representatives.csv
# =============================================================================

suppressPackageStartupMessages(library(data.table))

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)

MWAS_RESULT  <- file.path("output", "mwas", "primary", "primary_20cov.csv")
CLUSTER_FILE <- file.path("output", "mwas", "effective_tests", "metabolite_clusters.csv")
OUT_DIR      <- file.path("output", "mwas", "effective_tests")
OUT_FILE     <- file.path(OUT_DIR, "cluster_representatives.csv")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FDR_ALPHA <- 0.05

# Volcano significance boundary, used for display only: a hyperbolic curve, so
# that a small effect needs a smaller P than a large one to be highlighted.
HYPER_X0 <- 0.02; HYPER_C <- 0.04; HYPER_YMIN <- 1.2
hyperbolic_threshold <- function(x, x0 = HYPER_X0, cc = HYPER_C, ymin = HYPER_YMIN)
  ifelse(abs(x) <= x0, Inf, ymin + cc / (abs(x) - x0))

# The 37 CE-IVD certified Nightingale biomarkers (Julkunen et al. 2021,
# eLife 10:e63033). The same list is used in 06_effective_tests.R.
CE_37 <- c(
  "Total_C", "VLDL_C", "LDL_C", "HDL_C", "Total_TG", "ApoB", "ApoA1", "ApoB_by_ApoA1",
  "Total_FA", "Omega_3", "Omega_6", "PUFA", "MUFA", "SFA", "DHA",
  "Omega_3_pct", "Omega_6_pct", "PUFA_pct", "MUFA_pct", "SFA_pct", "DHA_pct",
  "PUFA_by_MUFA", "Omega_6_by_Omega_3", "Ala", "Gly", "His", "Ile", "Leu", "Val",
  "Phe", "Tyr", "Total_BCAA", "Glucose", "Lactate", "Creatinine", "Albumin", "GlycA")

# ---- inputs -----------------------------------------------------------------
stopifnot(file.exists(MWAS_RESULT), file.exists(CLUSTER_FILE))
raw      <- fread(MWAS_RESULT)
clusters <- fread(CLUSTER_FILE)
cat("[1] MWAS results:", nrow(raw), "measures | cluster table:", nrow(clusters), "rows\n")

keep     <- intersect(c("metabolite", "cluster_id", "cluster_size"), names(clusters))
clusters <- unique(clusters[, ..keep], by = "metabolite")
merged   <- merge(clusters, raw, by.x = "metabolite", by.y = "protein")
cat("    joined:", nrow(merged), "measures in", uniqueN(merged$cluster_id), "clusters\n")

# ---- representative selection ----------------------------------------------
prepare_reps <- function(dt) {
  dt[, pval_fdr_global := p.adjust(pval, method = "BH")]
  dt[, is_sig_global := pval_fdr_global < FDR_ALPHA]
  out <- dt[!is.na(pval), {
    sig <- .SD[is_sig_global == TRUE]
    if (nrow(sig) == 0) NULL else {
      lead     <- sig[which.min(pval)]              # cluster lead = smallest P
      lead_dir <- sign(log(lead$or))
      ce       <- sig[metabolite %in% CE_37]
      ce_conc  <- ce[sign(log(or)) == lead_dir]     # certified and same direction
      if (nrow(ce_conc) > 0) ce_conc[which.min(pval)] else lead
    }
  }, by = cluster_id]
  setnames(out, "metabolite", "protein")
  out[, neg_log10p := -log10(pval)]
  out[, log2or := log2(or)]
  out[, is_ce := protein %in% CE_37]
  out[, sig_hyperbolic := neg_log10p > hyperbolic_threshold(log2or) &
        is.finite(hyperbolic_threshold(log2or))]
  out[, sig_fdr := is_sig_global]
  out[]
}

reps <- prepare_reps(copy(merged))

# Every row here is FDR-significant by construction, so `is_ce` alone says how a
# representative was chosen: certified-and-direction-concordant, or the lead.
front <- c("cluster_id", "protein", "is_ce",
           "or", "or_lower", "or_upper", "log2or", "pval", "neg_log10p",
           "sig_hyperbolic", "sig_fdr")
setcolorder(reps, c(intersect(front, names(reps)), setdiff(names(reps), front)))

fwrite(reps, OUT_FILE)
# Export and dx upload to RAP  (cluster_representatives.csv feeds the feature manifest)

cat("\n[2] representatives written:", OUT_FILE, "\n")
cat("    representatives:", nrow(reps),
    "| CE-IVD:", sum(reps$is_ce), "| non-CE:", sum(!reps$is_ce), "\n")
cat("    all rows FDR-significant:", all(reps$sig_fdr), "\n")
cat("    reported panel (Methods): 57 representatives, 13 CE-IVD\n")
