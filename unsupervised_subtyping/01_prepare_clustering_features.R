#!/usr/bin/env Rscript
# =============================================================================
# 01_prepare_clustering_features.R  (subtyping entry point)
#
# Prepare case-only MOFA+ input:
#   1. load ml_ready/{g1,g2,g3}_full.rds (read-only)
#   2. filter to taste_2w_strict == 1 cases
#   3. partition the final manifest into 5 MOFA views:
#        olink         (Gaussian)   PWAS proteins
#        nmr           (Gaussian)   MWAS metabolites
#        clinical_cont (Gaussian)   continuous/ordinal covariates + ExWAS
#        clinical_bin  (Bernoulli)  binary ExWAS/covariate features
#        phecode       (Bernoulli)  PheCode binary indicators
#   4. z-score Gaussian views (params fit on G1 cases, reused for G2/G3)
#   5. keep binaries as 0/1 with NA preserved (MOFA+ native NA handling)
#
# Reads ml_ready/*.rds READ-ONLY; no imputation. Optional CLI arg "M"/"F" runs
# the sex-stratified analysis and suffixes outputs.
#
# Output: subtyping/inputs/clustering_input_g{1,2,3}{,_m,_f}.rds + view_partition.csv +
#         zscore_params_g1*.rds + feature_whitelist_g1*.rds; subtyping/reports/view_summary + prepare_clustering_log
# =============================================================================

rm(list = ls()); gc()
suppressPackageStartupMessages(library(data.table))

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
set.seed(20260413)

ARGS <- commandArgs(trailingOnly = TRUE)
SEX  <- if (length(ARGS) >= 1 && nzchar(ARGS[1])) toupper(ARGS[1]) else ""
stopifnot(SEX %in% c("", "M", "F"))
SUF      <- if (nzchar(SEX)) paste0("_", tolower(SEX)) else ""
SEX_CODE <- if (SEX == "M") 1L else if (SEX == "F") 0L else NA_integer_

ML_DIR  <- "output/ml_ready"
MF_DIR  <- "output/feature_manifest"
P7_DIR  <- "output/subtyping"
RPT_DIR <- "output/subtyping/reports"
dir.create(file.path(P7_DIR, "inputs"), showWarnings = FALSE, recursive = TRUE)
dir.create(RPT_DIR, showWarnings = FALSE, recursive = TRUE)
while (sink.number() > 0) sink()
sink(file.path(RPT_DIR, sprintf("prepare_clustering_log%s.txt", SUF)), split = TRUE)
on.exit(sink(), add = TRUE)

cat("=== prepare clustering features (MOFA+ input) ===\n")
cat(sprintf("sex stratification: %s\n", if (nzchar(SEX)) SEX else "none (all cases)"))

# ---- [1] build view partition -----------------------------------------------
cat("[1] building view partition\n")
manifest <- fread(file.path(MF_DIR, "master_feature_manifest_final.csv"))
stopifnot(all(c("feature_id", "source_analysis", "var_type", "var_role") %in% names(manifest)))

outcome_fids <- manifest[var_role == "outcome_primary", feature_id]
pc_fids      <- manifest[grepl("^PC[0-9]+$", feature_id), feature_id]
# smell_time / smell_extent are excluded from the clustering input and replaced
# by smell_2w_strict_OR, the composite that mirrors the taste_2w_strict
# OR-logic. That feature is derived once in 02_build_ml_ready.R, so it is
# already a manifest row and a column of the matrices read below; it lands in
# the clinical_cont view like any other continuous covariate.
# smell_time / smell_extent remain available for cluster characterization
# (04_cluster_discriminators) and as severity-orthogonality candidates in the severity-orthogonality analysis.
smell_severity_fids <- intersect(c("smell_time", "smell_extent"), manifest$feature_id)
excl_fids <- unique(c(outcome_fids, pc_fids, smell_severity_fids))

assign_view <- function(src, vtype, phe) {
  if (!is.na(src) && src == "PWAS") return("olink")
  if (!is.na(src) && src == "MWAS") return("nmr")
  if (!is.na(phe) && nzchar(as.character(phe))) return("phecode")
  vt <- tolower(ifelse(is.na(vtype) | vtype == "", "unknown", vtype))
  if (vt %in% c("binary", "derived_binary")) return("clinical_bin")
  "clinical_cont"
}
manifest[, view := mapply(assign_view, source_analysis, var_type, phecode)]
manifest[feature_id %in% excl_fids, view := "EXCLUDED"]
view_map <- manifest[view != "EXCLUDED", .(feature_id, view, source_analysis, var_type)]
if (!"smell_2w_strict_OR" %in% view_map$feature_id)
  stop("smell_2w_strict_OR is not in the manifest; run 02_build_ml_ready.R first")
fwrite(view_map, file.path(P7_DIR, "inputs", "view_partition.csv"))
cat("  view feature counts:\n"); print(view_map[, .N, by = view][order(-N)])

# ---- [2] load + filter cases + split into views -----------------------------
prepare_group <- function(grp, zparams = NULL, feature_whitelist = NULL) {
  cat(sprintf("[2] preparing %s\n", grp))
  rds_path <- file.path(ML_DIR, sprintf("%s_full.rds", grp))
  if (!file.exists(rds_path)) { cat("  SKIP (missing):", rds_path, "\n"); return(NULL) }
  dat <- as.data.table(readRDS(rds_path))
  if (!"taste_2w_strict" %in% names(dat)) stop("taste_2w_strict not found in ", rds_path)
  cases <- dat[taste_2w_strict == 1]
  if (!is.na(SEX_CODE)) {
    if (!"sex" %in% names(cases)) stop("sex column missing - cannot stratify")
    n_pre <- nrow(cases); cases <- cases[sex == SEX_CODE]
    cat(sprintf("  sex filter: kept %d/%d cases\n", nrow(cases), n_pre))
  }
  cat(sprintf("  total rows=%d, case rows=%d\n", nrow(dat), nrow(cases)))

  # smell_2w_strict_OR is derived once, in 02_build_ml_ready.R; deriving it a
  # second time here would be a copy of the same logic that could drift.
  if (!"smell_2w_strict_OR" %in% names(cases))
    stop("smell_2w_strict_OR missing from ", rds_path,
         " - run 02_build_ml_ready.R first")

  views <- split(view_map$feature_id, view_map$view); view_mats <- list()
  for (v in names(views)) {
    if (is.null(feature_whitelist)) {
      feats <- intersect(views[[v]], names(cases))
    } else {
      feats <- feature_whitelist[[v]]; if (is.null(feats) || length(feats) == 0) next
    }
    missing_feats <- setdiff(feats, names(cases))
    if (length(missing_feats) > 0) for (mf in missing_feats) cases[, (mf) := NA_real_]
    M <- as.matrix(cases[, feats, with = FALSE]); storage.mode(M) <- "numeric"
    if (is.null(feature_whitelist)) {  # filter only on G1; G2/G3 reuse the G1 set
      is_bern <- v %in% c("clinical_bin", "phecode")
      keep <- apply(M, 2, function(x) {
        if (sum(!is.na(x)) < 100) return(FALSE)
        s <- sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) return(FALSE)
        if (is_bern) { p <- mean(x, na.rm = TRUE); if (p < 0.01 || p > 0.99) return(FALSE) }
        TRUE
      })
      if (sum(!keep) > 0) cat(sprintf("    [%s] dropping %d feature(s) (obs<100 / zero-var / prev<1%%|>99%%)\n", v, sum(!keep)))
      M <- M[, keep, drop = FALSE]
    }
    if (ncol(M) == 0) next
    view_mats[[v]] <- M
  }

  # z-score Gaussian views (clinical_bin treated as Gaussian in MOFA -> z-scored too)
  gaussian_views <- intersect(c("olink", "nmr", "clinical_cont", "clinical_bin"), names(view_mats))
  if (is.null(zparams)) {
    zparams <- list()
    for (v in gaussian_views) {
      mu <- apply(view_mats[[v]], 2, mean, na.rm = TRUE); sd <- apply(view_mats[[v]], 2, sd, na.rm = TRUE)
      sd[!is.finite(sd) | sd == 0] <- 1; zparams[[v]] <- list(mu = mu, sd = sd)
    }
  }
  for (v in gaussian_views) {
    mu <- zparams[[v]]$mu[colnames(view_mats[[v]])]; sd <- zparams[[v]]$sd[colnames(view_mats[[v]])]
    Z <- sweep(sweep(view_mats[[v]], 2, mu, "-"), 2, sd, "/")
    Z[Z > 5] <- 5; Z[Z < -5] <- -5   # winsorise for numerical stability
    view_mats[[v]] <- Z
  }
  view_mats <- lapply(view_mats, t)   # features x samples (MOFA2 orientation)
  attr(view_mats, "eid") <- if ("eid" %in% names(cases)) cases$eid else seq_len(nrow(cases))
  attr(view_mats, "n_cases") <- nrow(cases); attr(view_mats, "zparams") <- zparams
  attr(view_mats, "group_name") <- grp
  view_mats
}

# ---- [3] run G1 (fits z-score params), apply to G2/G3 -----------------------
g1 <- prepare_group("group1")
zparams_g1 <- attr(g1, "zparams"); whitelist_g1 <- lapply(g1, rownames)
saveRDS(zparams_g1,   file.path(P7_DIR, "inputs", sprintf("zscore_params_g1%s.rds", SUF)))
saveRDS(whitelist_g1, file.path(P7_DIR, "inputs", sprintf("feature_whitelist_g1%s.rds", SUF)))
saveRDS(g1,           file.path(P7_DIR, "inputs", sprintf("clustering_input_g1%s.rds", SUF)))
for (grp in c("group2", "group3")) {
  gx <- prepare_group(grp, zparams = zparams_g1, feature_whitelist = whitelist_g1)
  if (!is.null(gx)) saveRDS(gx, file.path(P7_DIR, "inputs",
    sprintf("clustering_input_%s%s.rds", sub("group", "g", grp), SUF)))
}

# ---- [4] view-level summary -------------------------------------------------
summary_rows <- rbindlist(lapply(names(g1), function(v) {
  M <- g1[[v]]
  data.table(view = v, n_features = nrow(M), n_samples = ncol(M),
             pct_missing = round(100 * mean(is.na(M)), 2),
             likelihood = ifelse(v %in% c("olink", "nmr", "clinical_cont"), "gaussian", "bernoulli"))
}))
fwrite(summary_rows, file.path(RPT_DIR, sprintf("view_summary%s.csv", SUF)))
print(summary_rows)
cat("\n=== done ===\n")
