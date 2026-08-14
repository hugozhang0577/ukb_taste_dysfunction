#!/usr/bin/env Rscript
# =============================================================================
# UK Biobank NMR Metabolomics QC Pipeline
# Using ukbnmr package (Version 3 algorithm)
# =============================================================================
# 
# Input : the merged UK Biobank Category 220 / 221 / 222 extract (CSV or TSV)
# Output: NMR metabolomics data with technical variation removed
#
# Reference:
# Ritchie S.C. et al., Quality control and removal of technical variation of 
# NMR metabolic biomarker data in ~120,000 UK Biobank participants, 
# Sci Data 10, 64 (2023). doi: 10.1038/s41597-023-01949-y
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Environment
# -----------------------------------------------------------------------------

# Clear the workspace
rm(list = ls())
gc()

# Session options
options(stringsAsFactors = FALSE)
options(datatable.print.nrows = 20)

# Record the start time
start_time <- Sys.time()
cat("=== UK Biobank NMR QC Pipeline ===\n")
cat("start time:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n\n")

# -----------------------------------------------------------------------------
# 1. Packages
# -----------------------------------------------------------------------------

cat(">>> Step 1: loading packages...\n")

# Dependencies are checked, not installed: a run should not silently change the
# environment it is reproduced in. ukbnmr carries the version-3 biomarker QC
# algorithm (Ritchie et al. 2023) and is on CRAN.
missing_pkgs <- setdiff(c("ukbnmr", "data.table"),
                        rownames(installed.packages()))
if (length(missing_pkgs) > 0)
  stop("missing package(s): ", paste(missing_pkgs, collapse = ", "),
       "\n  install them first, e.g. install.packages(c(",
       paste(sprintf('"%s"', missing_pkgs), collapse = ", "), "))")

suppressPackageStartupMessages({
  library(ukbnmr)
  library(data.table)
})

cat("ukbnmr version:", as.character(packageVersion("ukbnmr")), "\n")
cat("data.table version:", as.character(packageVersion("data.table")), "\n\n")

# -----------------------------------------------------------------------------
# 2. Configuration
# -----------------------------------------------------------------------------

cat(">>> Step 2: configuration...\n")

# ============== paths (set these before running) ==============
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
INPUT_FILE <- file.path(PROJECT_DIR, "input", "metabolomics", "nmr_merged_all.csv")
OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "metabolomics_qc")

# ============== optional parameters ==============
# ukbnmr QC algorithm version (1, 2 or 3)
# Version 3 is the current algorithm and is the one intended for the full ~500k release
QC_VERSION <- 3L

# Whether to drop samples on outlier plates
REMOVE_OUTLIER_PLATES <- TRUE

# Field separator of the input file
FILE_SEP <- ","  # "," for CSV, "\t" for TSV

# ===========================================

# Create the output directory
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
  cat("created output directory:", OUTPUT_DIR, "\n")
}

cat("input file:", INPUT_FILE, "\n")
cat("output directory:", OUTPUT_DIR, "\n")
cat("QC algorithm version:", QC_VERSION, "\n")
cat("remove outlier plates:", REMOVE_OUTLIER_PLATES, "\n\n")

# -----------------------------------------------------------------------------
# 3. Read the data
# -----------------------------------------------------------------------------

cat(">>> Step 3: reading data...\n")

# Check that the input exists
if (!file.exists(INPUT_FILE)) {
  stop("input file does not exist - ", INPUT_FILE)
}

# Read
cat("reading (this can take several minutes)...\n")
exported <- fread(INPUT_FILE, sep = FILE_SEP, header = TRUE)

cat("dimensions:", nrow(exported), "rows x", ncol(exported), "columns\n")
cat("samples:", nrow(exported), "\n\n")

# -----------------------------------------------------------------------------
# 4. Validate the extract
# -----------------------------------------------------------------------------

cat(">>> Step 4: validating the extract...\n")

# Check for the eid column
if (!"eid" %in% names(exported)) {
  stop("the data has no 'eid' column")
}
cat("[OK] eid column present\n")

# Check the column-name format
sample_cols <- head(names(exported)[names(exported) != "eid"], 5)
cat("example column names:", paste(sample_cols, collapse = ", "), "\n")

# Category 222 holds the sample-processing fields the QC correction needs
cat("\nchecking sample-processing fields (Category 222):\n")
required_222_fields <- c(
  "p23649" = "Shipment Plate",
  "p23650" = "Spectrometer",
  "p23651" = "Measurement Quality Flagged",
  "p23652" = "High Lactate",
  "p23653" = "High Pyruvate",
  "p23654" = "Low Glucose",
  "p23655" = "Low Protein",
  "p23658" = "Sample Measured Date and Time",
  "p23659" = "Sample Prepared Date and Time",
  "p23660" = "Well position within plate",
  "p20282" = "Processing batch",
  "p20283" = "Resolved plate swaps"
)

missing_fields <- c()
for (field_id in names(required_222_fields)) {
  field_name <- required_222_fields[field_id]
  col_i0 <- paste0(field_id, "_i0")
  
  if (any(grepl(paste0("^", field_id), names(exported)))) {
    cat("  [OK]", field_id, "(", field_name, ")\n")
  } else {
    cat("  [FAIL]", field_id, "(", field_name, ") - MISSING\n")
    missing_fields <- c(missing_fields, field_id)
  }
}

if (length(missing_fields) > 0) {
  warning("some sample-processing fields are missing; the QC correction will be incomplete")
}

# Check the Category 220 biomarker fields
cat("\nchecking biomarker fields (Category 220):\n")
# 3-Hydroxybutyrate (Field 23474) is used as the probe field
biomarker_example <- "p23474"
if (any(grepl(paste0("^", biomarker_example), names(exported)))) {
  cat("  [OK] biomarker fields present (probe: ", biomarker_example, ")\n")
} else {
  stop("biomarker fields are missing")
}

# Count the instance-0 and instance-1 fields
n_i0 <- sum(grepl("_i0$", names(exported)))
n_i1 <- sum(grepl("_i1$", names(exported)))
cat("\n  Instance 0 (baseline) fields:", n_i0, "\n")
cat("  Instance 1 (repeat) fields:", n_i1, "\n\n")

# -----------------------------------------------------------------------------
# 5. Remove technical variation
# -----------------------------------------------------------------------------

cat(">>> Step 5: removing technical variation...\n")
cat("this step can take over an hour...\n")
cat("memory required: at least 32 GB RAM\n\n")

qc_start <- Sys.time()

# The ukbnmr QC call
processed <- remove_technical_variation(
  x = exported,
  remove.outlier.plates = REMOVE_OUTLIER_PLATES,
  version = QC_VERSION
)

qc_end <- Sys.time()
qc_duration <- difftime(qc_end, qc_start, units = "mins")
cat("\nQC complete; elapsed:", round(qc_duration, 1), "min\n\n")

# -----------------------------------------------------------------------------
# 6. Output structure
# -----------------------------------------------------------------------------

cat(">>> Step 6: output structure...\n\n")

cat("the returned object has these components:\n")
cat("  $biomarkers              - biomarker concentrations after QC correction\n")
cat("  $biomarker_qc_flags      - per-biomarker QC flags\n")
cat("  $sample_processing       - sample-processing information\n")
cat("  $log_offset              - offset applied before log transformation\n")
cat("  $outlier_plate_detection - outlier-plate detection\n\n\n")

# Overview of the biomarker table
cat("biomarker table dimensions:", 
    nrow(processed$biomarkers), "rows x", 
    ncol(processed$biomarkers), "columns\n")

# Distribution of visit_index
visit_table <- table(processed$biomarkers$visit_index)
cat("  Visit 0 (Baseline):", visit_table["0"], "samples\n")
if ("1" %in% names(visit_table)) {
  cat("  Visit 1 (Repeat):", visit_table["1"], "samples\n")
}

# First few biomarker column names
biomarker_cols <- names(processed$biomarkers)[!names(processed$biomarkers) %in% c("eid", "visit_index")]
cat("\nexample biomarker columns:\n")
cat(" ", paste(head(biomarker_cols, 10), collapse = ", "), "...\n\n")

# -----------------------------------------------------------------------------
# 7. Write results
# -----------------------------------------------------------------------------

cat(">>> Step 7: writing results...\n")

# Output file paths
output_files <- list(
  biomarkers = file.path(OUTPUT_DIR, "nmr_biomarkers_qc.csv"),
  biomarker_qc_flags = file.path(OUTPUT_DIR, "nmr_biomarker_qc_flags.csv"),
  sample_processing = file.path(OUTPUT_DIR, "nmr_sample_processing.csv"),
  log_offset = file.path(OUTPUT_DIR, "nmr_log_offset.csv"),
  outlier_plates = file.path(OUTPUT_DIR, "nmr_outlier_plate_detection.csv")
)

# Write each component
# Export and dx upload to RAP  (the QC'd biomarker table and its QC companions)
cat("writing biomarkers...\n")
fwrite(processed$biomarkers, output_files$biomarkers)
cat("  ->", output_files$biomarkers, "\n")

cat("writing biomarker_qc_flags...\n")
fwrite(processed$biomarker_qc_flags, output_files$biomarker_qc_flags)
cat("  ->", output_files$biomarker_qc_flags, "\n")

cat("writing sample_processing...\n")
fwrite(processed$sample_processing, output_files$sample_processing)
cat("  ->", output_files$sample_processing, "\n")

cat("writing log_offset...\n")
fwrite(processed$log_offset, output_files$log_offset)
cat("  ->", output_files$log_offset, "\n")

cat("writing outlier_plate_detection...\n")
fwrite(processed$outlier_plate_detection, output_files$outlier_plates)
cat("  ->", output_files$outlier_plates, "\n\n")

# -----------------------------------------------------------------------------
# 8. QC summary report
# -----------------------------------------------------------------------------

cat(">>> Step 8: QC summary report...\n\n")

# Build the summary
summary_report <- list()

# Basic counts
summary_report$total_samples <- nrow(processed$biomarkers)
summary_report$total_biomarkers <- length(biomarker_cols)
summary_report$visit0_samples <- as.numeric(visit_table["0"])
summary_report$visit1_samples <- ifelse("1" %in% names(visit_table), as.numeric(visit_table["1"]), 0)

# Missingness
missing_stats <- sapply(processed$biomarkers[, ..biomarker_cols], function(x) sum(is.na(x)))
summary_report$avg_missing_rate <- mean(missing_stats) / nrow(processed$biomarkers) * 100

# Outlier plates
if (!is.null(processed$outlier_plate_detection)) {
  n_outlier_plates <- sum(processed$outlier_plate_detection$outlier_plate, na.rm = TRUE)
  summary_report$outlier_plates <- n_outlier_plates
}

# Print the summary
cat("========== QC summary ==========\n")
cat("total samples:", summary_report$total_samples, "\n")
cat("  - Baseline (Visit 0):", summary_report$visit0_samples, "\n")
cat("  - Repeat (Visit 1):", summary_report$visit1_samples, "\n")
cat("biomarkers:", summary_report$total_biomarkers, "\n")
cat("mean missingness:", round(summary_report$avg_missing_rate, 2), "%\n")
if (!is.null(summary_report$outlier_plates)) {
  cat("outlier plates detected:", summary_report$outlier_plates, "\n")
}
cat("QC algorithm version:", QC_VERSION, "\n")
cat("==================================\n\n")

# Write the summary
summary_file <- file.path(OUTPUT_DIR, "qc_summary_report.txt")
sink(summary_file)
cat("UK Biobank NMR QC Summary Report\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=====================================\n\n")
cat("Input file:", INPUT_FILE, "\n")
cat("QC algorithm version:", QC_VERSION, "\n")
cat("Remove outlier plates:", REMOVE_OUTLIER_PLATES, "\n\n")
cat("Results:\n")
cat("  Total samples:", summary_report$total_samples, "\n")
cat("  Baseline samples:", summary_report$visit0_samples, "\n")
cat("  Repeat samples:", summary_report$visit1_samples, "\n")
cat("  Total biomarkers:", summary_report$total_biomarkers, "\n")
cat("  Average missing rate:", round(summary_report$avg_missing_rate, 2), "%\n")
if (!is.null(summary_report$outlier_plates)) {
  cat("  Outlier plates detected:", summary_report$outlier_plates, "\n")
}
cat("\nOutput files:\n")
for (name in names(output_files)) {
  cat(" ", name, ":", output_files[[name]], "\n")
}
sink()
cat("summary written:", summary_file, "\n\n")

# -----------------------------------------------------------------------------
# 9. Optional: extended biomarker ratios
# -----------------------------------------------------------------------------

cat(">>> Step 9: computing extended biomarker ratios (optional)...\n")

# Compute the 76 additional extended ratios
extended_ratios <- compute_extended_ratios(processed$biomarkers)
extended_ratio_qc <- compute_extended_ratio_qc_flags(processed$biomarker_qc_flags)

# Write the extended ratios
extended_file <- file.path(OUTPUT_DIR, "nmr_extended_ratios.csv")
fwrite(extended_ratios, extended_file)
cat("extended ratios written:", extended_file, "\n")

extended_qc_file <- file.path(OUTPUT_DIR, "nmr_extended_ratio_qc_flags.csv")
fwrite(extended_ratio_qc, extended_qc_file)
cat("extended-ratio QC flags written:", extended_qc_file, "\n\n")

# -----------------------------------------------------------------------------
# 10. Done
# -----------------------------------------------------------------------------

end_time <- Sys.time()
total_duration <- difftime(end_time, start_time, units = "mins")

cat("===========================================\n")
cat("QC pipeline complete\n")
cat("total elapsed:", round(total_duration, 1), "min\n")
cat("end time:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("===========================================\n\n")

cat("output files:\n")
for (f in list.files(OUTPUT_DIR, full.names = TRUE)) {
  file_size <- file.size(f) / 1024 / 1024  # MB
  cat("  ", basename(f), "-", round(file_size, 1), "MB\n")
}

cat("\nnext steps:\n")
cat("1. review qc_summary_report.txt\n")
cat("2. use nmr_biomarkers_qc.csv for the downstream analysis\n")
cat("3. to additionally adjust for biological covariates, see the ukbnmr vignette\n")
cat("4. cite: Ritchie SC et al. Sci Data 10, 64 (2023)\n\n")

# Release memory
rm(exported)
gc()

cat("Done!\n")

