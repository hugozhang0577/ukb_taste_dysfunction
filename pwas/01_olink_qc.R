################################################################################
#                    UKB Olink Proteomics QC Pipeline
#                              
#  Purpose: Comprehensive quality control for UK Biobank Olink proteomics data
#           with missingness mechanism diagnosis and stratified imputation
#  
#  Key Features:
#    1. Distinguishes True Missing (row doesn't exist) vs Below-LOD (result < LOD)
#    2. Diagnoses missingness mechanism (MCAR vs MNAR)
#    3. Stratified imputation based on missingness type:
#       - Below-LOD (MNAR, left-censored): MinProb method
#       - True Missing (likely MAR): Protein-wise median
#
#  QC Strategy:
#    1. Filter samples to phenotype cohort first
#    2. Remove proteins with >30% missing
#    3. Remove samples with >50% missing (primarily Batch 7 half-panel)
#    4. Diagnose missingness mechanism
#    5. Stratified imputation
#
#  References:
#    - Lind et al., Eur J Prev Cardiol 2024
#    - Jin et al., Sci Rep 2021 (imputation comparison)
#    - Lazar et al., J Proteome Res 2016 (MCAR/MNAR)
#    - Wei et al., PLoS One 2020 (GSimp/missForest comparison for Olink)
################################################################################

# ============================================================================
# 1. Setup and Configuration
# ============================================================================
cat("================================================================================\n")
cat("UKB Olink Proteomics QC Pipeline\n")
cat("================================================================================\n\n")
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gridExtra)
})

# ============================================================================
# 2. Configuration Parameters
# ============================================================================
# One cohort per run: set COHORT to group1 (discovery), group2 or group3 and run
# again. Imputation is fitted within the cohort being processed, so the cohorts
# are deliberately not pooled here.
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
COHORT      <- Sys.getenv("COHORT", unset = "group1")
if (!COHORT %in% c("group1", "group2", "group3"))
  stop("COHORT must be group1, group2 or group3; got: ", COHORT)

config <- list(
  cohort = COHORT,

  # Input files
  olink_file = file.path(PROJECT_DIR, "input", "proteomics", "olink_full_dataset.csv"),
  pheno_file = file.path(PROJECT_DIR, "input", "analysis_ready",
                         sprintf("phenotype_%s.csv", COHORT)),

  # QC thresholds
  protein_miss_threshold = 0.30,   # Remove proteins with >30% missing
  sample_miss_threshold = 0.50,    # Remove samples with >50% missing
  
  # Analysis settings
  instance = 0,                    # Use baseline instance only
  
  # Imputation settings
  mnar_method = "MinProb",         # For below-LOD values (left-censored MNAR)
  mar_method = "median",           # For true missing values (likely MAR)
  minprob_q = 0.01,                # Quantile for MinProb imputation
  minprob_sd_scale = 0.3,          # SD scaling factor for MinProb
  
  # Random seed for reproducibility
  seed = 42,
  
  # Output. output_prefix carries the cohort, so the QC'd matrix is written
  # under the name 02_run_pwas.R reads back (proteomics_<cohort>.csv) and the
  # three cohorts cannot overwrite one another.
  output_dir = file.path(PROJECT_DIR, "output", "pwas", "olink_qc"),
  output_prefix = sprintf("proteomics_%s", COHORT)
)

# Create output directory
dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

# Set seed for reproducibility
set.seed(config$seed)

cat("[Configuration]\n")
cat(sprintf("  Olink data: %s\n", config$olink_file))
cat(sprintf("  Phenotype file: %s\n", config$pheno_file))
cat(sprintf("  Protein threshold: <=%.0f%% missing\n", config$protein_miss_threshold * 100))
cat(sprintf("  Sample threshold: <=%.0f%% missing\n", config$sample_miss_threshold * 100))
cat(sprintf("  Below-LOD imputation: %s (q=%.2f, sd_scale=%.1f)\n", 
            config$mnar_method, config$minprob_q, config$minprob_sd_scale))
cat(sprintf("  True missing imputation: %s\n", config$mar_method))
cat(sprintf("  Random seed: %d\n", config$seed))
cat(sprintf("  Output directory: %s\n", config$output_dir))
cat("================================================================================\n\n")

# ============================================================================
# 3. Initialize QC Report
# ============================================================================
qc_report <- list(
  start_time = Sys.time(),
  config = config,
  steps = list()
)

# ============================================================================
# 4. Load Data
# ============================================================================
cat("========== Step 1: Loading Data ==========\n")

# Load phenotype file first
cat("Loading phenotype data...\n")
pheno <- fread(config$pheno_file)
pheno_eids <- unique(pheno$eid)
cat(sprintf("  Phenotype samples: %s\n", format(length(pheno_eids), big.mark = ",")))

# Load Olink data
cat("Loading Olink data (this may take a few minutes)...\n")
prot_data <- fread(config$olink_file, 
                   drop = "V1",
                   na.strings = c("", "NA", "NaN"))
cat(sprintf("  Raw data: %s rows\n", format(nrow(prot_data), big.mark = ",")))

# Filter to instance 0 only
prot_data <- prot_data[ins_index == config$instance]
cat(sprintf("  Instance %d: %s rows\n", config$instance, format(nrow(prot_data), big.mark = ",")))

qc_report$steps$raw_data <- list(
  total_rows = nrow(prot_data),
  total_samples = uniqueN(prot_data$eid),
  total_proteins = uniqueN(prot_data$protein_id)
)

# ============================================================================
# 5. Filter to Phenotype Cohort
# ============================================================================
cat("\n========== Step 2: Filter to Phenotype Cohort ==========\n")

n_samples_before <- uniqueN(prot_data$eid)
prot_data <- prot_data[eid %in% pheno_eids]
n_samples_after <- uniqueN(prot_data$eid)
n_samples_not_in_olink <- length(pheno_eids) - n_samples_after

cat(sprintf("  Samples before filter: %s\n", format(n_samples_before, big.mark = ",")))
cat(sprintf("  Samples after filter: %s\n", format(n_samples_after, big.mark = ",")))
cat(sprintf("  Samples removed (not in phenotype): %s\n", format(n_samples_before - n_samples_after, big.mark = ",")))
cat(sprintf("  Phenotype samples not in Olink: %d\n", n_samples_not_in_olink))

qc_report$steps$pheno_filter <- list(
  samples_before = n_samples_before,
  samples_after = n_samples_after,
  samples_removed = n_samples_before - n_samples_after,
  pheno_not_in_olink = n_samples_not_in_olink
)

# ============================================================================
# 6. Analyze Sample Characteristics (Batch Distribution)
# ============================================================================
cat("\n========== Step 3: Sample Characteristics ==========\n")

# Count proteins per sample
sample_stats <- prot_data[, .(
  n_proteins = .N,
  Batch = unique(Batch)[1]
), by = eid]

cat("\n[Proteins per Sample]\n")
print(summary(sample_stats$n_proteins))

cat("\n[Batch Distribution]\n")
batch_dist <- sample_stats[, .(
  n_samples = .N,
  mean_proteins = round(mean(n_proteins), 0),
  min_proteins = min(n_proteins),
  max_proteins = max(n_proteins)
), by = Batch][order(Batch)]
print(batch_dist)

# Identify Batch 7 (half-panel) issue
batch7_samples <- sample_stats[Batch == 7, .N]
batch7_mean_prot <- sample_stats[Batch == 7, mean(n_proteins)]
cat(sprintf("\n[Note] Batch 7 has %d samples with mean %.0f proteins (half-panel)\n", 
            batch7_samples, batch7_mean_prot))

qc_report$steps$batch_analysis <- list(
  batch_distribution = batch_dist,
  batch7_samples = batch7_samples,
  batch7_mean_proteins = batch7_mean_prot
)

# ============================================================================
# 7. Analyze Below-LOD Values (BEFORE converting to wide format)
# ============================================================================
cat("\n========== Step 4: Below-LOD Analysis ==========\n")

# Calculate below-LOD statistics per protein
lod_stats <- prot_data[, .(
  n_total = .N,
  n_below_lod = sum(result < LOD, na.rm = TRUE),
  below_lod_rate = sum(result < LOD, na.rm = TRUE) / .N,
  mean_result = mean(result, na.rm = TRUE),
  median_lod = median(LOD, na.rm = TRUE)
), by = protein_id]

cat("\n[Below-LOD Rate Distribution Across Proteins]\n")
print(summary(lod_stats$below_lod_rate))

cat("\n[Below-LOD Breakdown]\n")
cat(sprintf("  Proteins with 0%% below-LOD: %d\n", sum(lod_stats$below_lod_rate == 0)))
cat(sprintf("  Proteins with <5%% below-LOD: %d\n", sum(lod_stats$below_lod_rate < 0.05)))
cat(sprintf("  Proteins with <10%% below-LOD: %d\n", sum(lod_stats$below_lod_rate < 0.10)))
cat(sprintf("  Proteins with ≥10%% below-LOD: %d\n", sum(lod_stats$below_lod_rate >= 0.10)))
cat(sprintf("  Proteins with ≥20%% below-LOD: %d\n", sum(lod_stats$below_lod_rate >= 0.20)))

total_below_lod <- sum(lod_stats$n_below_lod)
total_measurements <- sum(lod_stats$n_total)
cat(sprintf("\n  Total below-LOD values: %s (%.2f%% of all measurements)\n", 
            format(total_below_lod, big.mark = ","),
            total_below_lod / total_measurements * 100))

qc_report$steps$lod_analysis <- list(
  lod_rate_summary = summary(lod_stats$below_lod_rate),
  total_below_lod = total_below_lod,
  total_measurements = total_measurements,
  pct_below_lod = total_below_lod / total_measurements * 100
)

# ============================================================================
# 8. Create Complete Matrix and Identify True Missing
# ============================================================================
cat("\n========== Step 5: Identify True Missing Values ==========\n")

all_eids <- unique(prot_data$eid)
all_proteins <- unique(prot_data$protein_id)

# Expected complete combinations
n_expected <- length(all_eids) * length(all_proteins)
n_existing <- nrow(prot_data)
n_true_missing <- n_expected - n_existing

cat(sprintf("  Expected measurements (complete): %s\n", format(n_expected, big.mark = ",")))
cat(sprintf("  Existing measurements: %s\n", format(n_existing, big.mark = ",")))
cat(sprintf("  True missing (row doesn't exist): %s (%.2f%%)\n", 
            format(n_true_missing, big.mark = ","),
            n_true_missing / n_expected * 100))

# Calculate true missing per protein
existing_counts <- prot_data[, .(n_measured = .N), by = protein_id]
true_miss_per_protein <- merge(
  data.table(protein_id = all_proteins),
  existing_counts,
  by = "protein_id",
  all.x = TRUE
)
true_miss_per_protein[is.na(n_measured), n_measured := 0]
true_miss_per_protein[, n_true_missing := length(all_eids) - n_measured]
true_miss_per_protein[, true_missing_rate := n_true_missing / length(all_eids)]

cat("\n[True Missing Rate Distribution Across Proteins]\n")
print(summary(true_miss_per_protein$true_missing_rate))

qc_report$steps$true_missing <- list(
  n_expected = n_expected,
  n_existing = n_existing,
  n_true_missing = n_true_missing,
  true_missing_rate = n_true_missing / n_expected
)

# ============================================================================
# 9. Combine Missing Statistics and Diagnose Mechanism
# ============================================================================
cat("\n========== Step 6: Diagnose Missingness Mechanism ==========\n")

# Merge LOD stats with true missing stats
diagnosis <- merge(lod_stats, true_miss_per_protein[, .(protein_id, n_true_missing, true_missing_rate)],
                   by = "protein_id", all = TRUE)

# Fill NAs
diagnosis[is.na(n_below_lod), n_below_lod := 0]
diagnosis[is.na(n_true_missing), n_true_missing := 0]
diagnosis[is.na(below_lod_rate), below_lod_rate := 0]
diagnosis[is.na(true_missing_rate), true_missing_rate := 0]

# Calculate combined unreliable rate
diagnosis[, `:=`(
  n_unreliable = n_below_lod + n_true_missing,
  unreliable_rate = below_lod_rate + true_missing_rate
)]

# Classify missingness mechanism
diagnosis[, likely_mechanism := fcase(
  n_true_missing == 0 & n_below_lod == 0, "No Issues",
  n_true_missing == 0 & n_below_lod > 0, "Below-LOD Only (MNAR)",
  n_true_missing > 0 & n_below_lod == 0, "True Missing Only (MAR)",
  n_below_lod / (n_true_missing + n_below_lod) > 0.8, "Predominantly MNAR",
  n_true_missing / (n_true_missing + n_below_lod) > 0.8, "Predominantly MAR",
  default = "Mixed (MAR + MNAR)"
)]

cat("\n[Missingness Mechanism Classification]\n")
mechanism_table <- diagnosis[, .N, by = likely_mechanism][order(-N)]
print(mechanism_table)

# MCAR Test: correlation between mean intensity and unreliable rate
# If negative correlation: evidence of MNAR (low-intensity proteins have more issues)
mcar_test_data <- diagnosis[!is.na(mean_result) & unreliable_rate > 0]

if (nrow(mcar_test_data) >= 10) {
  mcar_cor_test <- cor.test(mcar_test_data$mean_result, 
                            mcar_test_data$unreliable_rate, 
                            method = "spearman")
  
  mcar_interpretation <- if (mcar_cor_test$p.value < 0.05 && mcar_cor_test$estimate < -0.2) {
    "Strong evidence of MNAR (low-intensity proteins have more missing/below-LOD values)"
  } else if (mcar_cor_test$p.value < 0.05 && mcar_cor_test$estimate < 0) {
    "Moderate evidence of MNAR"
  } else if (mcar_cor_test$p.value >= 0.05) {
    "No significant correlation - consistent with MCAR"
  } else {
    "Inconclusive"
  }
  
  cat(sprintf("\n[MCAR Test: Intensity vs Unreliable Rate]\n"))
  cat(sprintf("  Spearman rho: %.3f\n", mcar_cor_test$estimate))
  cat(sprintf("  P-value: %.2e\n", mcar_cor_test$p.value))
  cat(sprintf("  Interpretation: %s\n", mcar_interpretation))
  
  qc_report$steps$mcar_test <- list(
    rho = mcar_cor_test$estimate,
    p_value = mcar_cor_test$p.value,
    interpretation = mcar_interpretation
  )
} else {
  cat("\n[MCAR Test] Insufficient data for test\n")
  qc_report$steps$mcar_test <- list(
    rho = NA,
    p_value = NA,
    interpretation = "Insufficient data"
  )
}

qc_report$steps$mechanism_diagnosis <- list(
  mechanism_table = mechanism_table,
  total_below_lod = sum(diagnosis$n_below_lod),
  total_true_missing = sum(diagnosis$n_true_missing)
)

# ============================================================================
# 10. Convert to Wide Format
# ============================================================================
cat("\n========== Step 7: Convert to Wide Format ==========\n")

# Create wide format for result values
prot_wide <- dcast(prot_data[, .(eid, protein_id, result)], 
                   eid ~ protein_id, 
                   value.var = "result")

# Also create wide format for LOD values (needed for imputation)
lod_wide <- dcast(prot_data[, .(eid, protein_id, LOD)], 
                  eid ~ protein_id, 
                  value.var = "LOD",
                  fun.aggregate = function(x) median(x, na.rm = TRUE))

n_samples_wide <- nrow(prot_wide)
n_proteins_wide <- ncol(prot_wide) - 1

cat(sprintf("  Wide format: %s samples × %d proteins\n", 
            format(n_samples_wide, big.mark = ","), n_proteins_wide))

# ============================================================================
# 11. Calculate Overall Missing Rates for QC
# ============================================================================
cat("\n========== Step 8: Calculate Missing Rates for QC ==========\n")

protein_cols <- setdiff(names(prot_wide), "eid")

# Total missing = True Missing + Below-LOD (for QC purposes, we count both)
# But for wide format, NA only represents True Missing
# Below-LOD values are present but < LOD

# Protein-level: combine true missing rate with below-LOD rate
protein_total_miss <- merge(
  diagnosis[, .(protein_id, below_lod_rate, true_missing_rate)],
  data.table(protein_id = protein_cols),
  by = "protein_id"
)
protein_total_miss[, total_issue_rate := below_lod_rate + true_missing_rate]

cat("\n[Total Issue Rate (Missing + Below-LOD) per Protein]\n")
print(summary(protein_total_miss$total_issue_rate))

# For QC filtering, we use true missing rate (since below-LOD can still be used)
protein_miss_rate <- sapply(prot_wide[, ..protein_cols], function(x) mean(is.na(x)))

cat("\n[True Missing Rate per Protein (for QC)]\n")
print(summary(protein_miss_rate))

# Sample-level missing rate
sample_miss_rate <- rowMeans(is.na(prot_wide[, ..protein_cols]))

cat("\n[Sample Missing Rate Distribution]\n")
print(summary(sample_miss_rate))

# ============================================================================
# 12. Apply QC Filters
# ============================================================================
cat("\n========== Step 9: Apply QC Filters ==========\n")

# Step A: Filter proteins (based on true missing rate)
keep_proteins <- names(protein_miss_rate)[protein_miss_rate <= config$protein_miss_threshold]
removed_proteins <- names(protein_miss_rate)[protein_miss_rate > config$protein_miss_threshold]

cat(sprintf("\n[Protein QC]\n"))
cat(sprintf("  Threshold: <=%.0f%% true missing\n", config$protein_miss_threshold * 100))
cat(sprintf("  Before: %d proteins\n", length(protein_cols)))
cat(sprintf("  Removed: %d proteins\n", length(removed_proteins)))
cat(sprintf("  Retained: %d proteins (%.1f%%)\n", 
            length(keep_proteins), length(keep_proteins)/length(protein_cols)*100))

# Step B: Recalculate sample missing rate based on retained proteins
prot_filtered <- prot_wide[, c("eid", keep_proteins), with = FALSE]
lod_filtered <- lod_wide[, c("eid", keep_proteins), with = FALSE]
sample_miss_rate_new <- rowMeans(is.na(prot_filtered[, ..keep_proteins]))

# Step C: Filter samples
keep_sample_idx <- sample_miss_rate_new <= config$sample_miss_threshold
keep_samples <- prot_filtered$eid[keep_sample_idx]
removed_samples <- prot_filtered$eid[!keep_sample_idx]

cat(sprintf("\n[Sample QC]\n"))
cat(sprintf("  Threshold: <=%.0f%% missing\n", config$sample_miss_threshold * 100))
cat(sprintf("  Before: %d samples\n", nrow(prot_filtered)))
cat(sprintf("  Removed: %d samples\n", length(removed_samples)))
cat(sprintf("  Retained: %d samples (%.1f%%)\n", 
            length(keep_samples), length(keep_samples)/nrow(prot_filtered)*100))

# Apply sample filter
prot_qc <- prot_filtered[eid %in% keep_samples]
lod_qc <- lod_filtered[eid %in% keep_samples]

cat(sprintf("\n[QC Result]\n"))
cat(sprintf("  Final data: %d samples × %d proteins\n", nrow(prot_qc), length(keep_proteins)))

# Analyze removed samples by batch
removed_sample_batch <- sample_stats[eid %in% removed_samples, .N, by = Batch][order(Batch)]
cat("\n[Removed Samples by Batch]\n")
print(removed_sample_batch)

# Update diagnosis to only include kept proteins
diagnosis_qc <- diagnosis[protein_id %in% keep_proteins]

qc_report$steps$qc_filter <- list(
  protein_threshold = config$protein_miss_threshold,
  sample_threshold = config$sample_miss_threshold,
  proteins_before = length(protein_cols),
  proteins_removed = length(removed_proteins),
  proteins_retained = length(keep_proteins),
  samples_before = nrow(prot_wide),
  samples_removed = length(removed_samples),
  samples_retained = length(keep_samples),
  removed_by_batch = removed_sample_batch
)

# ============================================================================
# 13. Identify Values to Impute (Post-QC)
# ============================================================================
cat("\n========== Step 10: Identify Values to Impute ==========\n")

# For each cell in the QC'd matrix, determine if it needs imputation and why:
# 1. True Missing (NA in wide format) - row didn't exist in original data
# 2. Below-LOD - value exists but is below LOD

# Count true missing (NA)
n_true_missing_qc <- sum(is.na(prot_qc[, ..keep_proteins]))

# Count below-LOD (value < LOD, not NA)
n_below_lod_qc <- 0
for (col in keep_proteins) {
  result_vec <- prot_qc[[col]]
  lod_vec <- lod_qc[[col]]
  # Below-LOD: not NA and result < LOD
  below_lod_idx <- !is.na(result_vec) & !is.na(lod_vec) & (result_vec < lod_vec)
  n_below_lod_qc <- n_below_lod_qc + sum(below_lod_idx)
}

total_cells_qc <- nrow(prot_qc) * length(keep_proteins)

cat(sprintf("  Total data points: %s\n", format(total_cells_qc, big.mark = ",")))
cat(sprintf("  True missing (NA): %s (%.2f%%)\n", 
            format(n_true_missing_qc, big.mark = ","),
            n_true_missing_qc / total_cells_qc * 100))
cat(sprintf("  Below-LOD: %s (%.2f%%)\n", 
            format(n_below_lod_qc, big.mark = ","),
            n_below_lod_qc / total_cells_qc * 100))
cat(sprintf("  Total to impute: %s (%.2f%%)\n", 
            format(n_true_missing_qc + n_below_lod_qc, big.mark = ","),
            (n_true_missing_qc + n_below_lod_qc) / total_cells_qc * 100))

qc_report$steps$imputation_targets <- list(
  total_cells = total_cells_qc,
  n_true_missing = n_true_missing_qc,
  n_below_lod = n_below_lod_qc,
  pct_true_missing = n_true_missing_qc / total_cells_qc * 100,
  pct_below_lod = n_below_lod_qc / total_cells_qc * 100
)

# ============================================================================
# 14. Stratified Imputation
# ============================================================================
cat("\n========== Step 11: Stratified Imputation ==========\n")

cat(sprintf("  Below-LOD (MNAR) method: %s\n", config$mnar_method))
cat(sprintf("  True missing (MAR) method: %s\n", config$mar_method))

# Create a copy for imputation
prot_imputed <- copy(prot_qc)

# Track imputation details
imputation_log <- data.table(
  protein_id = keep_proteins,
  n_true_missing_imputed = 0L,
  n_below_lod_imputed = 0L,
  method_true_missing = config$mar_method,
  method_below_lod = config$mnar_method
)

cat("\n  Processing proteins...\n")
pb_interval <- ceiling(length(keep_proteins) / 10)

for (i in seq_along(keep_proteins)) {
  col <- keep_proteins[i]
  result_vec <- prot_imputed[[col]]
  lod_vec <- lod_qc[[col]]
  
  # Identify below-LOD values (not NA, but < LOD)
  below_lod_idx <- !is.na(result_vec) & !is.na(lod_vec) & (result_vec < lod_vec)
  n_below_lod <- sum(below_lod_idx)
  
  # Identify true missing (NA)
  true_missing_idx <- is.na(result_vec)
  n_true_missing <- sum(true_missing_idx)
  
  # ===== Impute Below-LOD values (MNAR - left-censored) =====
  if (n_below_lod > 0) {
    observed <- result_vec[!is.na(result_vec) & !below_lod_idx]
    
    if (length(observed) >= 5) {
      if (config$mnar_method == "MinProb") {
        # MinProb: sample from left tail of distribution
        q_val <- quantile(observed, config$minprob_q, na.rm = TRUE)
        sigma <- sd(observed, na.rm = TRUE) * config$minprob_sd_scale
        imputed_values <- rnorm(n_below_lod, mean = q_val, sd = sigma)
        prot_imputed[[col]][below_lod_idx] <- imputed_values
      } else {
        # Fallback: LOD/2
        lod_median <- median(lod_vec, na.rm = TRUE)
        prot_imputed[[col]][below_lod_idx] <- lod_median / 2
      }
    } else {
      # Not enough observed values, use LOD/2
      lod_median <- median(lod_vec, na.rm = TRUE)
      prot_imputed[[col]][below_lod_idx] <- lod_median / 2
    }
    
    imputation_log[i, n_below_lod_imputed := n_below_lod]
  }
  
  # ===== Impute True Missing values (MAR) =====
  if (n_true_missing > 0) {
    if (config$mar_method == "median") {
      # Use values after below-LOD imputation for calculating median
      current_values <- prot_imputed[[col]][!true_missing_idx]
      impute_value <- median(current_values, na.rm = TRUE)
    } else {
      current_values <- prot_imputed[[col]][!true_missing_idx]
      impute_value <- mean(current_values, na.rm = TRUE)
    }
    
    prot_imputed[[col]][true_missing_idx] <- impute_value
    imputation_log[i, n_true_missing_imputed := n_true_missing]
  }
  
  # Progress indicator
  if (i %% pb_interval == 0) {
    cat(sprintf("    Processed %d/%d proteins (%.0f%%)\n", 
                i, length(keep_proteins), i/length(keep_proteins)*100))
  }
}

# Calculate totals
imputation_log[, n_total_imputed := n_true_missing_imputed + n_below_lod_imputed]
imputation_log[, pct_imputed := round(n_total_imputed / nrow(prot_imputed) * 100, 2)]

# Verify no missing values remain
n_missing_after <- sum(is.na(prot_imputed[, ..keep_proteins]))

cat(sprintf("\n[Imputation Summary]\n"))
cat(sprintf("  Below-LOD values imputed (%s): %s\n", 
            config$mnar_method,
            format(sum(imputation_log$n_below_lod_imputed), big.mark = ",")))
cat(sprintf("  True missing values imputed (%s): %s\n", 
            config$mar_method,
            format(sum(imputation_log$n_true_missing_imputed), big.mark = ",")))
cat(sprintf("  Total values imputed: %s\n", 
            format(sum(imputation_log$n_total_imputed), big.mark = ",")))
cat(sprintf("  Missing after imputation: %d\n", n_missing_after))

qc_report$steps$imputation <- list(
  mnar_method = config$mnar_method,
  mar_method = config$mar_method,
  n_below_lod_imputed = sum(imputation_log$n_below_lod_imputed),
  n_true_missing_imputed = sum(imputation_log$n_true_missing_imputed),
  n_total_imputed = sum(imputation_log$n_total_imputed),
  missing_after = n_missing_after
)

# ============================================================================
# 15. Update Phenotype File
# ============================================================================
cat("\n========== Step 12: Update Phenotype File ==========\n")

pheno_qc <- pheno[eid %in% keep_samples]
cat(sprintf("  Original phenotype samples: %d\n", nrow(pheno)))
cat(sprintf("  QC'd phenotype samples: %d\n", nrow(pheno_qc)))

# ============================================================================
# 16. Save Outputs
# ============================================================================
cat("\n========== Step 13: Save Outputs ==========\n")

# 1. QC'd wide-format data (imputed)
output_data_file <- file.path(config$output_dir, 
                              paste0(config$output_prefix, ".csv"))
fwrite(prot_imputed, output_data_file)
cat(sprintf("  Saved: %s\n", output_data_file))

# 2. QC'd phenotype file
output_pheno_file <- file.path(config$output_dir, 
                               paste0(config$output_prefix, "_pheno.csv"))
fwrite(pheno_qc, output_pheno_file)
cat(sprintf("  Saved: %s\n", output_pheno_file))

# 3. Removed proteins list
protein_miss_df <- data.table(
  protein_id = names(protein_miss_rate),
  true_missing_rate = as.numeric(protein_miss_rate)
)
protein_miss_df <- merge(protein_miss_df, 
                         diagnosis[, .(protein_id, below_lod_rate, likely_mechanism, mean_result)],
                         by = "protein_id", all.x = TRUE)
removed_proteins_df <- protein_miss_df[protein_id %in% removed_proteins][order(-true_missing_rate)]
output_removed_prot <- file.path(config$output_dir, paste0(config$output_prefix, "_removed_proteins.csv"))
fwrite(removed_proteins_df, output_removed_prot)
cat(sprintf("  Saved: %s\n", output_removed_prot))

# 4. Removed samples list
sample_miss_df <- data.table(
  eid = prot_wide$eid,
  miss_rate = sample_miss_rate
)
removed_samples_df <- merge(
  sample_miss_df[eid %in% removed_samples],
  sample_stats[eid %in% removed_samples, .(eid, n_proteins, Batch)],
  by = "eid"
)[order(-miss_rate)]
output_removed_samp <- file.path(config$output_dir, paste0(config$output_prefix, "_removed_samples.csv"))
fwrite(removed_samples_df, output_removed_samp)
cat(sprintf("  Saved: %s\n", output_removed_samp))

# 5. Imputation details
output_impute <- file.path(config$output_dir, paste0(config$output_prefix, "_imputation_details.csv"))
fwrite(imputation_log, output_impute)
cat(sprintf("  Saved: %s\n", output_impute))

# 6. Retained proteins list with diagnosis
retained_proteins_df <- protein_miss_df[protein_id %in% keep_proteins][order(true_missing_rate)]
output_retained_prot <- file.path(config$output_dir, paste0(config$output_prefix, "_retained_proteins.csv"))
fwrite(retained_proteins_df, output_retained_prot)
cat(sprintf("  Saved: %s\n", output_retained_prot))

# 7. Full missingness diagnosis
output_diagnosis <- file.path(config$output_dir, paste0(config$output_prefix, "_missingness_diagnosis.csv"))
fwrite(diagnosis, output_diagnosis)
# Export and dx upload to RAP  (the QC'd protein matrix, the matched phenotype
# file and the QC logs). Copy proteomics_<cohort>.csv into input/analysis_ready/
# before running 02_run_pwas.R, or point OLINK_DIR here.
cat(sprintf("  Saved: %s\n", output_diagnosis))

# ============================================================================
# 17. Generate QC Report
# ============================================================================
cat("\n========== Step 14: Generate QC Report ==========\n")

qc_report$end_time <- Sys.time()
qc_report$duration <- difftime(qc_report$end_time, qc_report$start_time, units = "mins")

# Create comprehensive report
report_text <- sprintf("
################################################################################
              UKB OLINK PROTEOMICS QC REPORT (RIGOROUS VERSION)
################################################################################

Date: %s
Duration: %.1f minutes
Random Seed: %d

================================================================================
1. CONFIGURATION
================================================================================
Input Files:
  - Olink data: %s
  - Phenotype file: %s

QC Thresholds:
  - Protein missing threshold: <=%.0f%%
  - Sample missing threshold: <=%.0f%%
  
Imputation Strategy (Stratified by Missingness Mechanism):
  - Below-LOD values (MNAR, left-censored): %s (q=%.2f, sd_scale=%.1f)
  - True missing values (MAR): %s

================================================================================
2. DATA FLOW SUMMARY
================================================================================

                    Samples         Proteins
                    -------         --------
Raw Olink data:     %s              %d
After pheno filter: %s              %d
After protein QC:   %s              %d
After sample QC:    %s              %d (FINAL)

================================================================================
3. MISSINGNESS MECHANISM ANALYSIS
================================================================================

[Definition]
  - True Missing: Sample-protein combination does not exist in data
    (row absent from long-format data)
  - Below-LOD: Measurement exists but result < limit of detection
    (left-censored, MNAR mechanism)

[Pre-QC Statistics]
  Total measurements expected: %s
  Actual measurements: %s
  True missing values: %s (%.2f%%)
  Below-LOD values: %s (%.2f%%)

[Missingness Mechanism per Protein]
%s

[MCAR Test: Mean Intensity vs Unreliable Rate]
  Spearman rho: %.3f
  P-value: %.2e
  Interpretation: %s

  Note: A significant negative correlation indicates MNAR (missing not at random),
  where low-intensity proteins are more likely to have missing/below-LOD values
  due to detection limits. This supports the use of left-censored imputation
  methods (MinProb) for below-LOD values.

================================================================================
4. REMOVED PROTEINS (%d total)
================================================================================
Reason: True missing rate > %.0f%%

Top 10 removed proteins:
%s

Full list saved to the removed-proteins log.

================================================================================
5. REMOVED SAMPLES (%d total)
================================================================================
Reason: Missing rate > %.0f%% (primarily Batch 7 half-panel samples)

Removed samples by Batch:
%s

Note: Batch 7 used Olink Explore 1536 panel (~1,460 proteins) while other 
batches used Olink Explore 3072 (~2,920 proteins). The 50%% threshold 
effectively separates these incompatible panel versions.

Full list saved to the removed-samples log.

================================================================================
6. IMPUTATION DETAILS
================================================================================
Stratified imputation was performed based on missingness mechanism:

[Below-LOD Values (MNAR - Left-censored)]
  Method: %s
  - For each protein, values below LOD were replaced with random draws
    from N(q%.0f, SD*%.1f), where q%.0f is the %.0fth percentile of 
    observed values and SD is the standard deviation of observed values.
  Values imputed: %s

[True Missing Values (MAR)]
  Method: Protein-wise %s
  - Missing values were imputed using the %s of observed values
    for each protein (calculated after below-LOD imputation).
  Values imputed: %s

[Total]
  Total values imputed: %s (%.2f%% of data matrix)
  Missing after imputation: %d

Imputation details per protein saved to the imputation log.

================================================================================
7. FINAL DATASET
================================================================================
Dimensions: %d samples × %d proteins

Output files:
  - %s.csv: QC'd and imputed protein matrix (wide format; the engine input)
  - %s_pheno.csv: Matched phenotype data
  - ..._removed_proteins.csv: Excluded proteins with missingness statistics
  - ..._removed_samples.csv: Excluded samples with batch information
  - ..._retained_proteins.csv: Retained proteins with missingness diagnosis
  - ..._imputation_details.csv: Imputation counts per protein by type
  - ..._missingness_diagnosis.csv: Full missingness mechanism diagnosis

================================================================================
8. METHODS SECTION (for manuscript)
================================================================================

Proteomics Data Processing and Quality Control

Plasma proteomics data were generated using the Olink Explore 3072 platform as 
part of the UK Biobank Pharma Proteomics Project. Normalized Protein eXpression 
(NPX) values, representing relative protein abundance on a log2 scale, were used 
for all analyses. Quality control was restricted to instance 0 (baseline visit) 
measurements from %s participants with available phenotype data.

Sample quality control identified that Batch 7 samples (N=%d) were assayed using 
an earlier panel version (Olink Explore 1536) with approximately half the protein 
coverage (~1,460 proteins vs ~2,920 proteins). To ensure data homogeneity, samples 
with >%.0f%% missing proteins were excluded (N=%d), effectively removing 
incompatible batch samples.

Protein quality control excluded assays with >%.0f%% missing values across samples 
(N=%d proteins removed), retaining %d proteins (%.1f%%) for downstream analysis.

Missingness mechanism was evaluated by examining the correlation between mean 
protein intensity and the proportion of unreliable values (missing + below limit 
of detection). A significant negative Spearman correlation (rho=%.3f, p=%.2e) 
indicated that low-abundance proteins were more likely to have missing or 
below-detection values, consistent with a missing-not-at-random (MNAR) mechanism 
for below-LOD values.

A stratified imputation approach was employed based on missingness mechanism:
(1) Below-LOD values (%.2f%% of retained data), representing left-censored MNAR 
data, were imputed using the MinProb method, which draws random values from a 
truncated normal distribution centered at the %.0fth percentile of observed values 
with a compressed standard deviation (scale factor=%.1f).
(2) True missing values (%.2f%% of retained data), likely missing-at-random (MAR), 
were imputed using protein-wise %s values.

After quality control and imputation, the final dataset comprised %d samples and 
%d proteins. Technical variables (PlateID, Batch) were retained for potential 
inclusion as covariates in downstream analyses.

================================================================================
9. REFERENCES
================================================================================

1. Sun BB, et al. Plasma proteomic associations with genetics and health in the 
   UK Biobank. Nature. 2023;622(7982):329-338.

2. Jin L, et al. A comparative study of evaluating missing value imputation 
   methods in label-free proteomics. Sci Rep. 2021;11(1):1760.

3. Lazar C, et al. Accounting for the Multiple Natures of Missing Values in 
   Label-Free Quantitative Proteomics Data Sets to Compare Imputation Strategies. 
   J Proteome Res. 2016;15(4):1116-1125.

4. Wei R, et al. Missing Value Imputation Approach for Mass Spectrometry-based 
   Metabolomics Data. Sci Rep. 2018;8(1):663.

================================================================================
",
                       # Header
                       format(Sys.time(), "%%Y-%%m-%%d %%H:%%M:%%S"),
                       as.numeric(qc_report$duration),
                       config$seed,
                       
                       # Configuration
                       config$olink_file,
                       config$pheno_file,
                       config$protein_miss_threshold * 100,
                       config$sample_miss_threshold * 100,
                       config$mnar_method, config$minprob_q, config$minprob_sd_scale,
                       config$mar_method,
                       
                       # Data flow
                       format(qc_report$steps$raw_data$total_samples, big.mark = ","),
                       qc_report$steps$raw_data$total_proteins,
                       format(qc_report$steps$pheno_filter$samples_after, big.mark = ","),
                       qc_report$steps$raw_data$total_proteins,
                       format(qc_report$steps$pheno_filter$samples_after, big.mark = ","),
                       qc_report$steps$qc_filter$proteins_retained,
                       format(qc_report$steps$qc_filter$samples_retained, big.mark = ","),
                       qc_report$steps$qc_filter$proteins_retained,
                       
                       # Missingness analysis
                       format(qc_report$steps$true_missing$n_expected, big.mark = ","),
                       format(qc_report$steps$true_missing$n_existing, big.mark = ","),
                       format(qc_report$steps$true_missing$n_true_missing, big.mark = ","),
                       qc_report$steps$true_missing$true_missing_rate * 100,
                       format(qc_report$steps$lod_analysis$total_below_lod, big.mark = ","),
                       qc_report$steps$lod_analysis$pct_below_lod,
                       paste(capture.output(print(mechanism_table)), collapse = "\n"),
                       
                       # MCAR test
                       ifelse(is.na(qc_report$steps$mcar_test$rho), NA, qc_report$steps$mcar_test$rho),
                       ifelse(is.na(qc_report$steps$mcar_test$p_value), NA, qc_report$steps$mcar_test$p_value),
                       qc_report$steps$mcar_test$interpretation,
                       
                       # Removed proteins
                       length(removed_proteins),
                       config$protein_miss_threshold * 100,
                       paste(capture.output(print(head(removed_proteins_df, 10))), collapse = "\n"),
                       
                       # Removed samples
                       length(removed_samples),
                       config$sample_miss_threshold * 100,
                       paste(capture.output(print(removed_sample_batch)), collapse = "\n"),
                       
                       # Imputation details
                       config$mnar_method,
                       config$minprob_q * 100, config$minprob_sd_scale,
                       config$minprob_q * 100, config$minprob_q * 100,
                       format(qc_report$steps$imputation$n_below_lod_imputed, big.mark = ","),
                       config$mar_method,
                       config$mar_method,
                       format(qc_report$steps$imputation$n_true_missing_imputed, big.mark = ","),
                       format(qc_report$steps$imputation$n_total_imputed, big.mark = ","),
                       qc_report$steps$imputation$n_total_imputed / total_cells_qc * 100,
                       qc_report$steps$imputation$missing_after,
                       
                       # Final dataset
                       nrow(prot_imputed),
                       length(keep_proteins),
                       config$output_prefix,
                       config$output_prefix,
                       
                       # Methods section
                       format(qc_report$steps$pheno_filter$samples_after, big.mark = ","),
                       batch7_samples,
                       config$sample_miss_threshold * 100,
                       length(removed_samples),
                       config$protein_miss_threshold * 100,
                       length(removed_proteins),
                       length(keep_proteins),
                       length(keep_proteins) / length(protein_cols) * 100,
                       ifelse(is.na(qc_report$steps$mcar_test$rho), 0, qc_report$steps$mcar_test$rho),
                       ifelse(is.na(qc_report$steps$mcar_test$p_value), 1, qc_report$steps$mcar_test$p_value),
                       qc_report$steps$imputation_targets$pct_below_lod,
                       config$minprob_q * 100,
                       config$minprob_sd_scale,
                       qc_report$steps$imputation_targets$pct_true_missing,
                       config$mar_method,
                       nrow(prot_imputed),
                       length(keep_proteins)
)

# Save report
output_report <- file.path(config$output_dir, "QC_Report.txt")
writeLines(report_text, output_report)
cat(sprintf("  Saved: %s\n", output_report))

# ============================================================================
# 18. Generate Diagnostic Plots
# ============================================================================
cat("\n========== Step 15: Generate Diagnostic Plots ==========\n")

pdf(file.path(config$output_dir, "QC_Diagnostic_Plots.pdf"), width = 12, height = 8)

# Plot 1: Sample protein count distribution
p1 <- ggplot(sample_stats, aes(x = n_proteins)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.8) +
  geom_vline(xintercept = 2923 * (1 - config$sample_miss_threshold), 
             color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 2923 * (1 - config$sample_miss_threshold) + 50, 
           y = Inf, label = "50% threshold", color = "red", vjust = 2, hjust = 0) +
  labs(title = "Distribution of Proteins per Sample (Before QC)",
       subtitle = "Red line indicates 50% threshold separating full-panel from half-panel samples",
       x = "Number of Proteins",
       y = "Number of Samples") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))
print(p1)

# Plot 2: Protein count by Batch
p2 <- ggplot(sample_stats, aes(x = factor(Batch), y = n_proteins, fill = factor(Batch))) +
  geom_boxplot(alpha = 0.8) +
  geom_hline(yintercept = 2923 * (1 - config$sample_miss_threshold), 
             color = "red", linetype = "dashed") +
  labs(title = "Proteins per Sample by Batch",
       subtitle = "Batch 7 shows distinct half-panel pattern (~1,460 proteins)",
       x = "Batch",
       y = "Number of Proteins") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))
print(p2)

# Plot 3: True Missing vs Below-LOD per protein
p3 <- ggplot(diagnosis, aes(x = true_missing_rate * 100, y = below_lod_rate * 100)) +
  geom_point(aes(color = likely_mechanism), alpha = 0.6, size = 2) +
  geom_vline(xintercept = config$protein_miss_threshold * 100, 
             color = "red", linetype = "dashed") +
  scale_color_brewer(palette = "Set1") +
  labs(title = "True Missing Rate vs Below-LOD Rate per Protein",
       subtitle = "Points right of red line were removed in QC",
       x = "True Missing Rate (%)",
       y = "Below-LOD Rate (%)",
       color = "Mechanism") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")
print(p3)

# Plot 4: MCAR test visualization
p4 <- ggplot(diagnosis[unreliable_rate > 0 & !is.na(mean_result)], 
             aes(x = mean_result, y = unreliable_rate * 100)) +
  geom_point(alpha = 0.4, color = "steelblue") +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(title = "MCAR Test: Mean Intensity vs Unreliable Rate",
       subtitle = sprintf("Spearman rho = %.3f, p = %.2e", 
                          ifelse(is.na(qc_report$steps$mcar_test$rho), NA, qc_report$steps$mcar_test$rho),
                          ifelse(is.na(qc_report$steps$mcar_test$p_value), NA, qc_report$steps$mcar_test$p_value)),
       x = "Mean Protein Intensity (NPX)",
       y = "Unreliable Rate (Missing + Below-LOD) (%)") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))
print(p4)

# Plot 5: Missingness mechanism distribution
mechanism_df <- as.data.frame(mechanism_table)
names(mechanism_df) <- c("Mechanism", "Count")
p5 <- ggplot(mechanism_df, aes(x = reorder(Mechanism, Count), y = Count, fill = Mechanism)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = Count), hjust = -0.1) +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Missingness Mechanism Classification",
       subtitle = "Per-protein classification based on missing type proportions",
       x = "",
       y = "Number of Proteins") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))
print(p5)

# Plot 6: Imputation summary
impute_long <- melt(imputation_log[n_total_imputed > 0, 
                                   .(protein_id, 
                                     `Below-LOD (MinProb)` = n_below_lod_imputed,
                                     `True Missing (Median)` = n_true_missing_imputed)],
                    id.vars = "protein_id",
                    variable.name = "Type",
                    value.name = "Count")
impute_long <- impute_long[Count > 0]

p6 <- ggplot(impute_long, aes(x = Count, fill = Type)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("Below-LOD (MinProb)" = "coral", 
                               "True Missing (Median)" = "steelblue")) +
  labs(title = "Distribution of Imputed Values per Protein by Type",
       subtitle = sprintf("Total: %s below-LOD + %s true missing = %s imputed",
                          format(sum(imputation_log$n_below_lod_imputed), big.mark = ","),
                          format(sum(imputation_log$n_true_missing_imputed), big.mark = ","),
                          format(sum(imputation_log$n_total_imputed), big.mark = ",")),
       x = "Number of Imputed Values",
       y = "Number of Proteins",
       fill = "Imputation Type") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")
print(p6)

# Plot 7: QC flow summary
qc_flow <- data.frame(
  Stage = factor(c("Raw", "Pheno Filter", "Protein QC", "Sample QC"),
                 levels = c("Raw", "Pheno Filter", "Protein QC", "Sample QC")),
  Samples = c(qc_report$steps$raw_data$total_samples,
              qc_report$steps$pheno_filter$samples_after,
              qc_report$steps$pheno_filter$samples_after,
              qc_report$steps$qc_filter$samples_retained),
  Proteins = c(qc_report$steps$raw_data$total_proteins,
               qc_report$steps$raw_data$total_proteins,
               qc_report$steps$qc_filter$proteins_retained,
               qc_report$steps$qc_filter$proteins_retained)
)

p7a <- ggplot(qc_flow, aes(x = Stage, y = Samples)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = format(Samples, big.mark = ",")), vjust = -0.5) +
  labs(title = "Sample Count Through QC Pipeline",
       y = "Number of Samples") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))

p7b <- ggplot(qc_flow, aes(x = Stage, y = Proteins)) +
  geom_bar(stat = "identity", fill = "coral", alpha = 0.8) +
  geom_text(aes(label = format(Proteins, big.mark = ",")), vjust = -0.5) +
  labs(title = "Protein Count Through QC Pipeline",
       y = "Number of Proteins") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))

grid.arrange(p7a, p7b, ncol = 2)

dev.off()
cat(sprintf("  Saved: QC_Diagnostic_Plots.pdf\n"))

# ============================================================================
# 19. Final Summary
# ============================================================================
cat("\n")
cat("################################################################################\n")
cat("#                    QC PIPELINE (RIGOROUS) COMPLETED                         #\n")
cat("################################################################################\n")
cat("\n")
cat(sprintf("Duration: %.1f minutes\n", as.numeric(qc_report$duration)))
cat(sprintf("Output directory: %s\n", config$output_dir))
cat("\n")
cat("=== FINAL DATASET ===\n")
cat(sprintf("Samples: %d\n", nrow(prot_imputed)))
cat(sprintf("Proteins: %d\n", length(keep_proteins)))
cat("\n")
cat("=== MISSINGNESS SUMMARY ===\n")
cat(sprintf("Below-LOD imputed (MinProb): %s\n", 
            format(qc_report$steps$imputation$n_below_lod_imputed, big.mark = ",")))
cat(sprintf("True missing imputed (Median): %s\n", 
            format(qc_report$steps$imputation$n_true_missing_imputed, big.mark = ",")))
cat("\n")
cat("=== FILES GENERATED ===\n")
cat(sprintf("  1. %s.csv - QC'd protein matrix (engine input)\n", config$output_prefix))
cat(sprintf("  2. %s_pheno.csv - Matched phenotype data\n", config$output_prefix))
cat(sprintf("  3. %s_removed_proteins.csv - Excluded proteins\n", config$output_prefix))
cat(sprintf("  4. %s_removed_samples.csv - Excluded samples\n", config$output_prefix))
cat(sprintf("  5. %s_retained_proteins.csv - Retained proteins with diagnosis\n", config$output_prefix))
cat(sprintf("  6. %s_imputation_details.csv - Per-protein imputation counts\n", config$output_prefix))
cat(sprintf("  7. %s_missingness_diagnosis.csv - Full missingness diagnosis\n", config$output_prefix))
cat("  8. QC_Report.txt - Comprehensive QC report\n")
cat("  9. QC_Diagnostic_Plots.pdf - Diagnostic visualizations\n")
cat("\n")
cat("================================================================================\n")

# Clean up
rm(prot_data, prot_wide, prot_filtered, lod_wide, lod_filtered)
gc()

cat("\nDone!\n")
