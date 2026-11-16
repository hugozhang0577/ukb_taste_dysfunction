# Preprocess the baseline NMR measures once in the pooled sample (values beyond 4 IQR from the
# median set to missing, log1p, z-score), so effect sizes are comparable across cohorts; then
# split by cohort.
library(data.table)
setwd(Sys.getenv("PROJECT_DIR", "."))

QC_DIR   <- "output/metabolomics_qc"
OUT_DIR  <- "output/mwas/derive"
IQR_MULT <- 4
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Returns the transformed column, or NULL (column dropped) when fewer than 100 values remain
# or the SD is zero. Columns with negative values are shifted to a minimum of 1e-6 before log1p
# (the extended ratios add the two shift terms one after the other).
preprocess <- function(x, ext) {
  ok <- !is.na(x)
  if (sum(ok) < 100) return(NULL)
  med <- median(x[ok]); iqr <- IQR(x[ok])
  if (iqr > 0) x[ok & (x < med - IQR_MULT * iqr | x > med + IQR_MULT * iqr)] <- NA
  ok <- !is.na(x)
  if (sum(ok) < 100) return(NULL)
  mn <- min(x[ok])
  x <- if (mn >= 0) log1p(x) else if (ext) log1p(x + abs(mn) + 1e-6) else log1p(x + (abs(mn) + 1e-6))
  ok <- !is.na(x); s <- sd(x[ok])
  if (sum(ok) < 100 || is.na(s) || s == 0) return(NULL)
  (x - mean(x[ok])) / s
}

preprocess_cols <- function(dt, ext = FALSE) {
  for (v in setdiff(names(dt), "eid")) {
    z <- preprocess(dt[[v]], ext)
    if (is.null(z)) dt[, (v) := NULL] else set(dt, j = v, value = z)
  }
  dt
}

baseline <- function(f) {
  d <- fread(f)[visit_index == 0]
  d[, visit_index := NULL]
}

nmr <- baseline(file.path(QC_DIR, "nmr_biomarkers_qc.csv"))
nmr <- preprocess_cols(nmr[!duplicated(eid)])

# Extended ratios not already in the biomarker table
ext <- baseline(file.path(QC_DIR, "nmr_extended_ratios.csv"))
ext <- ext[, setdiff(names(ext), setdiff(names(nmr), "eid")), with = FALSE]
nmr <- merge(nmr, preprocess_cols(ext, ext = TRUE), by = "eid", all.x = TRUE)

# Copied to input/analysis_ready/ for the scans and the subtype characterization
for (g in c("group1", "group2", "group3")) {
  ids <- fread(sprintf("input/analysis_ready/phenotype_%s.csv", g), select = "eid")$eid
  fwrite(nmr[eid %in% ids], file.path(OUT_DIR, sprintf("metabolomics_%s.csv", g)))
}
