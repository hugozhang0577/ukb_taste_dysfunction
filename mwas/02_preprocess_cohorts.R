# =============================================================================
# UKB NMR Metabolomics — Baseline Extraction & Preprocessing
# =============================================================================
#
# The script is written to be run either end-to-end or section by section.
#
# Design:
#   Mirrors the Olink PWAS pipeline. Omics and phenotype data are stored
#   separately and joined on eid at run time, so the phenotype files are shared
#   with the PWAS pipeline rather than duplicated.
#
# Input:
#   1. nmr_biomarkers_qc.csv        - output of the ukbnmr QC step (required)
#   2. nmr_biomarker_qc_flags.csv   - per-measurement QC flags (optional)
#   3. nmr_extended_ratios.csv      - extended ratios (optional)
#   4. phenotype_group{1,2,3}.csv   - only the eid column is read, to split the
#                                     matrix into the three cohorts
#
# Output under output/mwas/derive/ (eid plus the NMR biomarker columns only):
#   metabolomics_group{1,2,3}.csv
#   nmr_biomarker_colnames.csv    - biomarker column list
#   nmr_scale_params.csv          - standardisation parameters
#   nmr_log1p_shift_info.csv      - shift applied before log1p
#   nmr_outlier_summary.csv       - outlier counts
#   nmr_biomarker_summary.csv     - descriptive statistics
#   nmr_preprocessing_report.txt  - run log
#
# Reference:
#   Julkunen et al. Nat Commun 14, 604 (2023)
#   Ritchie et al. Sci Data 10, 64 (2023)
# =============================================================================


# 0. Environment ----

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(datatable.print.nrows = 20)

library(data.table)

start_time <- Sys.time()
cat("start time:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")


# 0a. Paths ----

# Output directory of the ukbnmr QC step
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = "."); NMR_QC_DIR <- file.path(PROJECT_DIR, "output", "metabolomics_qc")

NMR_BIOMARKERS_FILE      <- file.path(NMR_QC_DIR, "nmr_biomarkers_qc.csv")
NMR_QC_FLAGS_FILE        <- file.path(NMR_QC_DIR, "nmr_biomarker_qc_flags.csv")
NMR_EXTENDED_RATIOS_FILE <- file.path(NMR_QC_DIR, "nmr_extended_ratios.csv")

# Cohort membership; only the eid column is read. These are the same per-cohort
# phenotype files every other scan uses, so the cohort definition has one source.
COHORT_DIR <- file.path(PROJECT_DIR, "input", "analysis_ready")

COHORT_FILES <- list(
  group1 = file.path(COHORT_DIR, "phenotype_group1.csv"),
  group2 = file.path(COHORT_DIR, "phenotype_group2.csv"),
  group3 = file.path(COHORT_DIR, "phenotype_group3.csv")
)

# Output paths. Derived matrices are written under output/; copy the three
# metabolomics_group*.csv into input/analysis_ready/ before running the batches
# (mwas_config.R reads them from there).
OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "mwas", "derive")


# 0b. Parameters ----

IQR_MULTIPLIER       <- 4       # Outlier threshold in IQRs (Julkunen: 4)
USE_QC_FLAGS         <- FALSE    # Whether to mask low-quality measurements using the QC flags
USE_EXTENDED_RATIOS  <- TRUE    # Whether to append the 76 extended ratios


# 0c. Initialise ----

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

log_lines <- character()

write_log <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  log_lines <<- c(log_lines, msg)
}

write_log(paste(rep("=", 70), collapse = ""))
write_log("UKB NMR Metabolomics Preprocessing Report")
write_log("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
write_log(paste(rep("=", 70), collapse = ""))


# 1. Load the NMR data and keep the baseline visit ----

write_log("")
write_log(">>> Step 1: load NMR data, keep baseline")
write_log(paste(rep("-", 50), collapse = ""))

stopifnot(file.exists(NMR_BIOMARKERS_FILE))

cat("[INFO] reading NMR data (this can take several minutes)...\n")
nmr_raw <- fread(NMR_BIOMARKERS_FILE)
write_log("  raw data: ", nrow(nmr_raw), " x ", ncol(nmr_raw))

stopifnot("eid" %in% names(nmr_raw))
stopifnot("visit_index" %in% names(nmr_raw))

# Visit-instance distribution, then
visit_tab <- table(nmr_raw$visit_index, useNA = "ifany")
for (v in names(visit_tab)) {
  write_log("  visit_index=", v, ": ", visit_tab[v])
}

# keep baseline
nmr <- nmr_raw[visit_index == 0]
nmr[, visit_index := NULL]
write_log("  baseline samples: ", nrow(nmr))

# De-duplicate on eid
n_dup <- sum(duplicated(nmr$eid))
if (n_dup > 0) {
  write_log("  [WARNING] duplicate eid: ", n_dup, " -> keeping the first record")
  nmr <- nmr[!duplicated(eid)]
}

# Identify the biomarker columns
biomarker_cols <- setdiff(names(nmr), "eid")
write_log("  biomarkers: ", length(biomarker_cols))

rm(nmr_raw); gc()

# >>> Checkpoint ----
cat("\n[CHECK] inspect interactively:\n")
cat("  dim(nmr)                  # dimensions\n")
cat("  head(biomarker_cols, 20)  # biomarker column names\n")
cat("  nmr[1:5, 1:8]             # first rows and columns\n\n")


# 2. QC-flag masking ----

write_log("")
write_log(">>> Step 2: QC-flag masking")
write_log(paste(rep("-", 50), collapse = ""))

if (USE_QC_FLAGS && file.exists(NMR_QC_FLAGS_FILE)) {
  
  qc_flags <- fread(NMR_QC_FLAGS_FILE)
  if ("visit_index" %in% names(qc_flags)) {
    qc_flags <- qc_flags[visit_index == 0]
    qc_flags[, visit_index := NULL]
  }
  
  qc_biomarker_cols <- intersect(biomarker_cols, setdiff(names(qc_flags), "eid"))
  write_log("  biomarkers matched: ", length(qc_biomarker_cols))
  
  n_masked <- 0L
  if (length(qc_biomarker_cols) > 0) {
    setkey(nmr, eid)
    setkey(qc_flags, eid)
    
    for (bc in qc_biomarker_cols) {
      flagged_eids <- qc_flags$eid[nchar(qc_flags[[bc]]) > 0 & !is.na(qc_flags[[bc]])]
      nf <- length(flagged_eids)
      if (nf > 0) {
        nmr[eid %in% flagged_eids, (bc) := NA_real_]
        n_masked <- n_masked + nf
      }
    }
  }
  write_log("  measurements masked: ", n_masked)
  rm(qc_flags); gc()
  
} else {
  write_log("  skipped (USE_QC_FLAGS = FALSE, or the flag file is absent)")
}


# 3. Preprocessing ----
# 
# Three steps: outlier trimming -> log1p -> SD standardisation.
#
# All three run on the pooled baseline sample and only then is the data
# split by cohort, so the cohorts share one mean/SD and their per-SD odds

write_log("")
write_log(">>> Step 3: preprocessing on the pooled baseline")
write_log(paste(rep("-", 50), collapse = ""))
write_log("  IQR multiplier: ", IQR_MULTIPLIER)

pre_miss <- mean(sapply(nmr[, ..biomarker_cols],
                        function(x) mean(is.na(x)))) * 100
write_log("  mean missingness before preprocessing: ", round(pre_miss, 2), "%")


# 3a. Outlier trimming ----

cat("[INFO] 3a. outlier trimming (> ", IQR_MULTIPLIER, " x IQR)...\n")

outlier_summary <- data.table(
  biomarker = character(), n_valid = integer(),
  n_outlier = integer(), pct_outlier = numeric()
)

for (bc in biomarker_cols) {
  vals <- nmr[[bc]]
  ok <- !is.na(vals)
  if (sum(ok) < 100) next
  
  med <- median(vals[ok])
  iqr <- IQR(vals[ok])
  if (iqr == 0) next
  
  lo <- med - IQR_MULTIPLIER * iqr
  hi <- med + IQR_MULTIPLIER * iqr
  bad <- ok & (vals < lo | vals > hi)
  n_out <- sum(bad)
  
  if (n_out > 0) set(nmr, which(bad), bc, NA_real_)
  
  outlier_summary <- rbind(outlier_summary, data.table(
    biomarker = bc, n_valid = sum(ok),
    n_outlier = n_out,
    pct_outlier = round(n_out / sum(ok) * 100, 3)
  ))
}

write_log("  outliers set to NA: ", sum(outlier_summary$n_outlier))
write_log("  mean outlier rate: ", round(mean(outlier_summary$pct_outlier), 3), "%")

setorder(outlier_summary, -pct_outlier)
write_log("  highest outlier rates:")
for (i in seq_len(min(5, nrow(outlier_summary)))) {
  write_log("    ", outlier_summary$biomarker[i], ": ",
            outlier_summary$n_outlier[i], " (",
            outlier_summary$pct_outlier[i], "%)")
}

fwrite(outlier_summary, file.path(OUTPUT_DIR, "nmr_outlier_summary.csv"))


# 3b. log1p transformation ----

cat("[INFO] 3b. log1p transformation...\n")

# Values can be negative after the ukbnmr correction (they are residuals),
# so such columns are shifted to non-negative before log1p; the shift is kept

shift_info <- data.table(biomarker = character(),
                         min_val = numeric(), shift = numeric())

for (bc in biomarker_cols) {
  vals <- nmr[[bc]]
  valid_vals <- vals[!is.na(vals)]
  if (length(valid_vals) < 100) next
  
  min_val <- min(valid_vals)
  if (min_val < 0) {
    shift <- abs(min_val) + 1e-6
    set(nmr, j = bc, value = log1p(vals + shift))
    shift_info <- rbind(shift_info, data.table(
      biomarker = bc,
      min_val = round(min_val, 4),
      shift = round(shift, 4)
    ))
  } else {
    set(nmr, j = bc, value = log1p(vals))
  }
}

write_log("  biomarkers requiring a shift: ", nrow(shift_info), " / ", length(biomarker_cols))
fwrite(shift_info, file.path(OUTPUT_DIR, "nmr_log1p_shift_info.csv"))


# 3c. SD standardisation ----

cat("[INFO] 3c. SD standardisation (z-score)...\n")

scale_params <- data.table(
  biomarker = character(), mean = numeric(), sd = numeric(),
  n_valid = integer(), n_missing = integer(), pct_missing = numeric()
)

skipped_cols <- character()

for (bc in biomarker_cols) {
  vals <- nmr[[bc]]
  ok <- !is.na(vals)
  nv <- sum(ok)
  
  if (nv < 100) { skipped_cols <- c(skipped_cols, bc); next }
  
  m <- mean(vals[ok])
  s <- sd(vals[ok])
  if (s == 0 || is.na(s)) { skipped_cols <- c(skipped_cols, bc); next }
  
  set(nmr, j = bc, value = (vals - m) / s)
  
  nm <- sum(!ok)
  scale_params <- rbind(scale_params, data.table(
    biomarker = bc, mean = round(m, 6), sd = round(s, 6),
    n_valid = nv, n_missing = nm,
    pct_missing = round(nm / nrow(nmr) * 100, 2)
  ))
}

write_log("  standardised: ", nrow(scale_params), " biomarkers")
if (length(skipped_cols) > 0) {
  write_log("  skipped (n < 100 or sd = 0): ", length(skipped_cols))
  write_log("    ", paste(head(skipped_cols, 10), collapse = ", "))
  nmr[, (skipped_cols) := NULL]
  biomarker_cols <- setdiff(biomarker_cols, skipped_cols)
}

write_log("  mean missingness after preprocessing: ", round(mean(scale_params$pct_missing), 2), "%")
fwrite(scale_params, file.path(OUTPUT_DIR, "nmr_scale_params.csv"))

# >>> Checkpoint ----
cat("\n[CHECK] after standardisation:\n")
cat("  head(scale_params)                 # standardisation parameters\n")
cat("  summary(nmr[[biomarker_cols[1]]])  # should be close to mean 0, sd 1\n\n")


# 4. Extended ratios (optional) ----

write_log("")
write_log(">>> Step 4: extended ratios")
write_log(paste(rep("-", 50), collapse = ""))

extended_cols <- character()

if (USE_EXTENDED_RATIOS && file.exists(NMR_EXTENDED_RATIOS_FILE)) {
  
  ext <- fread(NMR_EXTENDED_RATIOS_FILE)
  if ("visit_index" %in% names(ext)) {
    ext <- ext[visit_index == 0]
    ext[, visit_index := NULL]
  }
  
  extended_cols <- setdiff(names(ext), "eid")
  write_log("  columns in the extended file: ", length(extended_cols))
  
  # Drop columns that duplicate the core panel (ukbnmr re-emits the originals)
  dup_cols <- intersect(biomarker_cols, extended_cols)
  if (length(dup_cols) > 0) {
    write_log("  duplicate columns dropped: ", length(dup_cols))
    ext[, (dup_cols) := NULL]
    extended_cols <- setdiff(extended_cols, dup_cols)
  }
  write_log("  ratios actually added: ", length(extended_cols))
  
  # Same preprocessing: outlier trimming -> log1p -> standardisation
  ext_skip <- character()
  ext_scale <- data.table(
    biomarker = character(), mean = numeric(), sd = numeric(), n_valid = integer()
  )
  
  for (bc in extended_cols) {
    vals <- ext[[bc]]
    ok <- !is.na(vals)
    if (sum(ok) < 100) { ext_skip <- c(ext_skip, bc); next }
    
    # outliers
    med <- median(vals[ok]); iqr <- IQR(vals[ok])
    if (iqr > 0) {
      bad <- ok & (vals < med - IQR_MULTIPLIER * iqr |
                     vals > med + IQR_MULTIPLIER * iqr)
      if (sum(bad) > 0) set(ext, which(bad), bc, NA_real_)
    }
    
    # log1p
    vals <- ext[[bc]]; vv <- vals[!is.na(vals)]
    if (length(vv) < 100) { ext_skip <- c(ext_skip, bc); next }
    if (min(vv) < 0) {
      set(ext, j = bc, value = log1p(vals + abs(min(vv)) + 1e-6))
    } else {
      set(ext, j = bc, value = log1p(vals))
    }
    
    # standardise
    vals <- ext[[bc]]; ok <- !is.na(vals); nv <- sum(ok)
    if (nv < 100) { ext_skip <- c(ext_skip, bc); next }
    m <- mean(vals[ok]); s <- sd(vals[ok])
    if (s == 0 || is.na(s)) { ext_skip <- c(ext_skip, bc); next }
    set(ext, j = bc, value = (vals - m) / s)
    ext_scale <- rbind(ext_scale, data.table(
      biomarker = bc, mean = round(m, 6), sd = round(s, 6), n_valid = nv
    ))
  }
  
  if (length(ext_skip) > 0) {
    ext[, (ext_skip) := NULL]
    extended_cols <- setdiff(extended_cols, ext_skip)
  }
  
  write_log("  standardised: ", length(extended_cols))
  fwrite(ext_scale, file.path(OUTPUT_DIR, "nmr_extended_scale_params.csv"))
  
  nmr <- merge(nmr, ext, by = "eid", all.x = TRUE)
  write_log("  columns after merging: ", ncol(nmr))
  
  rm(ext); gc()
  
} else {
  write_log("  skipped")
}

all_nmr_cols <- c(biomarker_cols, extended_cols)
write_log("")
write_log("  NMR measures in total: ", length(all_nmr_cols),
          " (core=", length(biomarker_cols),
          ", extended=", length(extended_cols), ")")

# 5. Split by cohort ----

write_log("")
write_log(">>> Step 5: split the NMR data by cohort")
write_log(paste(rep("-", 50), collapse = ""))

cohort_stats <- data.table(
  group = character(), n_cohort = integer(),
  n_overlap = integer(), pct_coverage = numeric()
)

for (gname in names(COHORT_FILES)) {
  
  fpath <- COHORT_FILES[[gname]]
  write_log("")
  write_log("  -- ", toupper(gname), " --")
  
  if (!file.exists(fpath)) {
    write_log("  [ERROR] file not found: ", fpath)
    next
  }
  
  # Read the eid column only
  header <- names(fread(fpath, nrows = 0))
  eid_col <- intersect(c("eid", "f.eid", "ID", "id", "EID"), header)[1]
  
  if (is.na(eid_col)) {
    write_log("  [ERROR] no eid column found")
    next
  }
  
  cohort_eids <- fread(fpath, select = eid_col)[[1]]
  write_log("  cohort N: ", length(cohort_eids))
  
  # Subset to eid plus the NMR columns
  nmr_sub <- nmr[eid %in% cohort_eids]
  n_overlap <- nrow(nmr_sub)
  pct <- round(n_overlap / length(cohort_eids) * 100, 1)
  write_log("  matched in NMR: ", n_overlap, " (", pct, "% of cohort)")
  
  if (n_overlap == 0) {
    write_log("  [WARNING] no matching participants; skipping")
    next
  }
  
  # Mean NMR completeness
  avg_cov <- mean(sapply(
    nmr_sub[, ..all_nmr_cols],
    function(x) mean(!is.na(x))
  )) * 100
  write_log("  mean NMR completeness: ", round(avg_cov, 1), "%")
  
  # Write
  # Name matches what mwas_config.R reads back; a mismatch here leaves every
  # batch script without its input.
  out_file <- file.path(OUTPUT_DIR, paste0("metabolomics_", gname, ".csv"))
  fwrite(nmr_sub, out_file)
  fsize <- round(file.size(out_file) / 1024 / 1024, 1)
  write_log("  written: ", basename(out_file), " (",
            nrow(nmr_sub), " x ", ncol(nmr_sub), ", ", fsize, " MB)")
  
  cohort_stats <- rbind(cohort_stats, data.table(
    group = gname, n_cohort = length(cohort_eids),
    n_overlap = n_overlap, pct_coverage = pct
  ))
}


# 6. Column list and CE-IVD flag ----

write_log("")
write_log(">>> Step 6: biomarker column list")
write_log(paste(rep("-", 50), collapse = ""))

# The 37 CE-IVD certified biomarkers; the same list is used in
# 06_effective_tests.R. Source: Julkunen et al. (2021) eLife
# 10:e63033 ("37 biomarkers in the panel have been certified for diagnostics
# use"); see also Julkunen 2023 Nat Commun. Linoleic acid (LA/LA_pct),
# Clinical_LDL_C, Gln and the ketone/glycolysis intermediates are research-use
# grade, are not among the certified 37, and are therefore not flagged here.
CE_BIOMARKERS <- c(
  # lipoprotein lipids (5)
  "Total_C", "VLDL_C", "LDL_C", "HDL_C", "Total_TG",
  # apolipoproteins (3)
  "ApoB", "ApoA1", "ApoB_by_ApoA1",
  # fatty acids, absolute concentration (7)
  "Total_FA", "Omega_3", "Omega_6", "PUFA", "MUFA", "SFA", "DHA",
  # fatty acids, percentages and ratios (8)
  "Omega_3_pct", "Omega_6_pct", "PUFA_pct", "MUFA_pct", "SFA_pct", "DHA_pct",
  "PUFA_by_MUFA", "Omega_6_by_Omega_3",
  # amino acids (9, including total BCAA)
  "Ala", "Gly", "His", "Ile", "Leu", "Val", "Phe", "Tyr", "Total_BCAA",
  # glycolysis (2)
  "Glucose", "Lactate",
  # fluid balance and inflammation (3)
  "Creatinine", "Albumin", "GlycA"
)

ce_matched <- intersect(all_nmr_cols, CE_BIOMARKERS)

colnames_dt <- data.table(
  biomarker = all_nmr_cols,
  type = ifelse(all_nmr_cols %in% biomarker_cols, "core", "extended_ratio"),
  is_ce_certified = all_nmr_cols %in% ce_matched
)

fwrite(colnames_dt, file.path(OUTPUT_DIR, "nmr_biomarker_colnames.csv"))

write_log("  core: ", sum(colnames_dt$type == "core"))
write_log("  extended_ratio: ", sum(colnames_dt$type == "extended_ratio"))
write_log("  CE certified (matched): ", sum(colnames_dt$is_ce_certified))
if (length(ce_matched) > 0) {
  write_log("  CE matched: ", paste(head(ce_matched, 15), collapse = ", "),
            if (length(ce_matched) > 15) " ..." else "")
}

# >>> Checkpoint ----
cat("\n[CHECK] CE-IVD matches (37 expected, Julkunen 2021 eLife):\n")
cat("  colnames_dt[is_ce_certified == TRUE]  # matched CE biomarkers\n")
cat("  unmatched CE names:", paste(setdiff(CE_BIOMARKERS, ce_matched), collapse = ", "), "\n")
cat("  if any are missing, check the column names and re-run this section\n\n")


# 7. Descriptive statistics ----

write_log("")
write_log(">>> Step 7: descriptive statistics (standardised, pooled baseline)")
write_log(paste(rep("-", 50), collapse = ""))

summary_dt <- rbindlist(lapply(all_nmr_cols, function(bc) {
  vals <- nmr[[bc]]
  ok <- !is.na(vals)
  nv <- sum(ok)
  if (nv == 0) return(NULL)
  qs <- quantile(vals[ok], probs = c(0, 0.25, 0.5, 0.75, 1))
  data.table(
    biomarker = bc, n_valid = nv,
    pct_missing = round(sum(!ok) / nrow(nmr) * 100, 2),
    mean = round(mean(vals[ok]), 4), sd = round(sd(vals[ok]), 4),
    min = round(qs[1], 4), q25 = round(qs[2], 4),
    median = round(qs[3], 4), q75 = round(qs[4], 4),
    max = round(qs[5], 4)
  )
}))

fwrite(summary_dt, file.path(OUTPUT_DIR, "nmr_biomarker_summary.csv"))
write_log("  written: nmr_biomarker_summary.csv (", nrow(summary_dt), " biomarkers)")


# 8. Summary ----

write_log("")
write_log(paste(rep("=", 70), collapse = ""))
write_log("                          Summary")
write_log(paste(rep("=", 70), collapse = ""))

write_log("")
write_log("Design:")
write_log("  omics and phenotype stored separately, as in the Olink PWAS pipeline")
write_log("  preprocessing on the pooled baseline, so cohort effect sizes are comparable")
write_log("  the output files are passed directly as the engine's feature matrix")

write_log("")
write_log("preprocessing:")
write_log("  outliers: > ", IQR_MULTIPLIER, " x IQR from median -> NA")
write_log("  transformation: log1p (negative columns shifted first)")
write_log("  standardisation: z-score (mean 0, sd 1)")

write_log("")
write_log("NMR measures: ", length(all_nmr_cols),
          " (core=", length(biomarker_cols),
          ", extended=", length(extended_cols), ")")

write_log("")
write_log("cohort coverage:")
write_log(sprintf("  %-10s %10s %10s %10s",
                  "Group", "Cohort_N", "NMR_N", "Coverage"))
write_log("  ", paste(rep("-", 45), collapse = ""))
for (i in seq_len(nrow(cohort_stats))) {
  r <- cohort_stats[i]
  write_log(sprintf("  %-10s %10d %10d %9.1f%%",
                    r$group, r$n_cohort, r$n_overlap, r$pct_coverage))
}

write_log("")
write_log("output files:")
# Export and dx upload to RAP  (the three per-cohort matrices and the
# preprocessing parameter tables). Copy metabolomics_group{1,2,3}.csv into
# input/analysis_ready/ before running the batch scripts.
for (f in sort(list.files(OUTPUT_DIR, full.names = TRUE))) {
  fsize <- round(file.size(f) / 1024 / 1024, 1)
  write_log("  ", basename(f), " (", fsize, " MB)")
}

end_time <- Sys.time()
dur <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 1)
write_log("")
write_log("elapsed: ", dur, " min")

# Write the log file (define LOG_FILE before use)
LOG_FILE <- file.path(OUTPUT_DIR, "nmr_preprocessing_report.txt")
writeLines(log_lines, LOG_FILE)
cat("\n")
cat("======================================================================\n")
cat("  Done\n")
cat("======================================================================\n\n")
cat("  output directory: ", OUTPUT_DIR, "\n")
cat("  log:              ", LOG_FILE, "\n\n")
cat("  downstream usage (same as the Olink PWAS pipeline):\n\n")
cat("  Rscript ../mwas_glm_analysis.R \\\n")
cat("    --olink      metabolomics_group1.csv \\\n")
cat("    --phenotype  group1_pure_white_british_full_med_oral_clinical.csv \\\n")
cat("    --outcome    taste_2w_strict \\\n")
cat("    --method     glm \\\n")
cat("    --covariates \"age_baseline,sex,PC1,PC2,PC3,PC4\" \\\n")
cat("    --output     mwas_N0_minimal.csv\n\n")

rm(nmr); gc()
