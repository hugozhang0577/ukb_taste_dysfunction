# MOFA+ input for one sex: cases split into five views (Olink, NMR, clinical continuous, clinical binary,
# prior diagnoses); features filtered and Gaussian views z-scored in group 1, and the same features and
# z-score parameters applied to groups 2 and 3.
library(data.table)
setwd(Sys.getenv("PROJECT_DIR", "."))
set.seed(20260413)

SEX      <- tolower(commandArgs(trailingOnly = TRUE)[1])   # Usage: Rscript unsupervised_subtyping/01_prepare_clustering_features.R M|F
SUF      <- paste0("_", SEX)
SEX_CODE <- if (SEX == "m") 1L else 0L
IN_DIR   <- "output/subtyping/inputs"
dir.create(IN_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- view partition (outcome, genetic PCs and smell duration/impact excluded)
manifest  <- fread("output/feature_manifest/master_feature_manifest_final.csv")
excl_fids <- unique(c(manifest[var_role == "outcome_primary", feature_id],
                      manifest[grepl("^PC[0-9]+$", feature_id), feature_id],
                      intersect(c("smell_time", "smell_extent"), manifest$feature_id)))
assign_view <- function(src, vtype, dx) {
  if (!is.na(src) && src == "PWAS") return("olink")
  if (!is.na(src) && src == "MWAS") return("nmr")
  if (!is.na(dx) && nzchar(as.character(dx))) return("prior_dx")
  vt <- tolower(ifelse(is.na(vtype) | vtype == "", "unknown", vtype))
  if (vt %in% c("binary", "derived_binary")) "clinical_bin" else "clinical_cont"
}
manifest[, view := mapply(assign_view, source_analysis, var_type, icd10)]
manifest[feature_id %in% excl_fids, view := "EXCLUDED"]
view_map <- manifest[view != "EXCLUDED", .(feature_id, view, source_analysis, var_type)]
fwrite(view_map, file.path(IN_DIR, "view_partition.csv"))

# ---- one cohort -> list of feature x sample matrices
prepare_group <- function(grp, zparams = NULL, whitelist = NULL) {
  cases <- as.data.table(readRDS(sprintf("output/ml_ready/%s_full.rds", grp)))[taste_2w_strict == 1 & sex == SEX_CODE]
  views <- split(view_map$feature_id, view_map$view); view_mats <- list()
  for (v in names(views)) {
    feats <- if (is.null(whitelist)) intersect(views[[v]], names(cases)) else whitelist[[v]]
    if (length(feats) == 0) next
    for (mf in setdiff(feats, names(cases))) cases[, (mf) := NA_real_]
    M <- as.matrix(cases[, feats, with = FALSE]); storage.mode(M) <- "numeric"
    if (is.null(whitelist)) {   # group 1: drop features observed in < 100, zero variance, or prevalence < 1% / > 99%
      is_bern <- v %in% c("clinical_bin", "prior_dx")
      keep <- apply(M, 2, function(x) {
        s <- sd(x, na.rm = TRUE); p <- mean(x, na.rm = TRUE)
        sum(!is.na(x)) >= 100 && is.finite(s) && s > 0 && (!is_bern || (p >= 0.01 && p <= 0.99))
      })
      M <- M[, keep, drop = FALSE]
    }
    if (ncol(M) > 0) view_mats[[v]] <- M
  }
  # z-score the Gaussian views (clinical_bin is Gaussian in MOFA+), winsorize at +/-5
  gaussian_views <- intersect(c("olink", "nmr", "clinical_cont", "clinical_bin"), names(view_mats))
  if (is.null(zparams)) zparams <- lapply(setNames(nm = gaussian_views), function(v) {
    s <- apply(view_mats[[v]], 2, sd, na.rm = TRUE); s[!is.finite(s) | s == 0] <- 1
    list(mu = apply(view_mats[[v]], 2, mean, na.rm = TRUE), sd = s)
  })
  for (v in gaussian_views) {
    cols <- colnames(view_mats[[v]])
    Z <- sweep(sweep(view_mats[[v]], 2, zparams[[v]]$mu[cols], "-"), 2, zparams[[v]]$sd[cols], "/")
    Z[Z > 5] <- 5; Z[Z < -5] <- -5
    view_mats[[v]] <- Z
  }
  view_mats <- lapply(view_mats, t)
  attr(view_mats, "eid") <- cases$eid
  attr(view_mats, "n_cases") <- nrow(cases); attr(view_mats, "zparams") <- zparams
  attr(view_mats, "group_name") <- grp
  view_mats
}

g1 <- prepare_group("group1")
saveRDS(g1, file.path(IN_DIR, sprintf("clustering_input_g1%s.rds", SUF)))
for (g in 2:3)
  saveRDS(prepare_group(paste0("group", g), attr(g1, "zparams"), lapply(g1, rownames)),
          file.path(IN_DIR, sprintf("clustering_input_g%d%s.rds", g, SUF)))
