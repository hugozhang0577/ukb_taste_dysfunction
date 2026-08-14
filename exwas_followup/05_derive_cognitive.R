#!/usr/bin/env Rscript
# ============================================================================
# Follow-up ExWAS — cognitive function (Category 116)
# ============================================================================
#
# Operates on the in-memory `raw` (+ `pheno_list`) loaded by
# 01_derive_mental_digestive.R.
#
# Input:   `raw` data.table (from shared framework Steps 1-3)
# Output:  followup_exwas_cognitive_function_{group}.csv
#          followup_exwas_cognitive_function_variable_dict.csv
#
# Strategy (Option A — individual tests, no composite):
#   Each cognitive test enters ExWAS as an independent continuous exposure.
#   Log-transformation applied to time-based and error measures per
#   Fawns-Ritchie & Deary (2020). No Z-standardization or reverse-coding
#   (unnecessary when tests are not combined into a composite).
#
# Instance selection:
#   i1 (2021 follow-up) preferred over i0 (2014) where available,
#   as i1 has 10-15% lower missingness and is temporally closer to
#   the taste outcome (2022-2023). Both precede the outcome.
#
# Key references:
#   - Fawns-Ritchie & Deary 2020 PLOS ONE — reliability/validity + transforms
#   - Hagenaars et al. 2016 Mol Psychiatry — PCA data reduction + stability
#   - Lazar et al. 2023 Front Psychol — three-factor SEM structure
# ============================================================================

cat("\n================================================================\n")
cat("  Cognitive function (Category 116)\n")
cat("================================================================\n\n")

# ============================================================================
# 1. FIELD REGISTRY
# ============================================================================

cog_tests <- list(
  
  pairs_matching = list(
    field = 20132,
    description = "Pairs matching (number of errors)",
    cols_i0 = c("p20132_i0_a0", "p20132_i0_a1", "p20132_i0_a2"),
    cols_i1 = NULL,  # Not available in i1 (replaced by Matrix in 2021)
    domain = "Visual memory",
    transform = "log(x+1)",
    direction = "lower = better"
  ),
  
  tmt_a = list(
    field = 20156,
    description = "Trail Making Test part A (seconds)",
    cols_i0 = "p20156_i0",
    cols_i1 = "p20156_i1",
    domain = "Processing speed",
    transform = "log(x)",
    direction = "lower = better"
  ),
  
  tmt_b = list(
    field = 20157,
    description = "Trail Making Test part B (seconds)",
    cols_i0 = "p20157_i0",
    cols_i1 = "p20157_i1",
    domain = "Processing speed + executive function",
    transform = "log(x)",
    direction = "lower = better"
  ),
  
  symbol_digit = list(
    field = 20159,
    description = "Symbol Digit Substitution (number correct)",
    cols_i0 = "p20159_i0",
    cols_i1 = "p20159_i1",
    domain = "Complex processing speed",
    transform = "none",
    direction = "higher = better"
  ),
  
  fluid_iq = list(
    field = 20191,
    description = "Fluid Intelligence / verbal-numerical reasoning (correct out of 13)",
    cols_i0 = "p20191_i0",
    cols_i1 = "p20191_i1",
    domain = "Verbal-numerical reasoning",
    transform = "none",
    direction = "higher = better"
  ),
  
  matrix_reasoning = list(
    field = 20760,
    description = "Matrix pattern completion (number of puzzles correctly solved)",
    cols_i0 = NULL,  # Not available in i0
    cols_i1 = "p20760_i1",
    domain = "Non-verbal reasoning",
    transform = "none",
    direction = "higher = better"
  )
)

cat("  Registered tests:", length(cog_tests), "\n")
register_fields("Cognitive function",
                vapply(cog_tests, function(x) x$field, numeric(1)))

# Verify columns present
all_cog_cols <- unlist(lapply(cog_tests, function(t) c(t$cols_i0, t$cols_i1)))
all_cog_cols <- all_cog_cols[!is.null(all_cog_cols)]
cog_missing <- all_cog_cols[!all_cog_cols %in% names(raw)]
if (length(cog_missing) > 0) {
  cat("  WARNING: Missing columns:", paste(cog_missing, collapse = ", "), "\n")
} else {
  cat("  All declared columns found in data [OK]\n")
}

# ============================================================================
# 2. WORKING COPY + INSTANCE SELECTION
# ============================================================================

cat("\n--- Step 2: Instance selection (i1 preferred, i0 fallback) ---\n")

cog <- raw[, c("eid", all_cog_cols), with = FALSE]

# --- Pairs Matching: aggregate across arrays (a0, a1, a2) ---
# Score = total number of incorrect attempts across all rounds
# Only i0 available; take sum of a0+a1+a2 (each is error count per round)
pm_cols <- cog_tests$pairs_matching$cols_i0
pm_present <- intersect(pm_cols, names(cog))

if (length(pm_present) > 0) {
  cog[, cog_pairs_errors_raw := rowSums(.SD, na.rm = FALSE),
     .SDcols = pm_present]
  n_pm <- sum(!is.na(cog$cog_pairs_errors_raw))
  cat("  Pairs Matching (i0, sum of", length(pm_present), "arrays): n_valid=",
      n_pm, "\n")
} else {
  cog[, cog_pairs_errors_raw := NA_real_]
  cat("  Pairs Matching: no columns found\n")
}

# --- TMT-A: i1 preferred, fallback to i0 ---
has_tmt_a_i1 <- "p20156_i1" %in% names(cog)
has_tmt_a_i0 <- "p20156_i0" %in% names(cog)
if (has_tmt_a_i1 && has_tmt_a_i0) {
  cog[, cog_tmt_a_raw := fifelse(!is.na(p20156_i1),
                                as.numeric(p20156_i1),
                                as.numeric(p20156_i0))]
} else if (has_tmt_a_i1) {
  cog[, cog_tmt_a_raw := as.numeric(p20156_i1)]
} else if (has_tmt_a_i0) {
  cog[, cog_tmt_a_raw := as.numeric(p20156_i0)]
} else {
  cog[, cog_tmt_a_raw := NA_real_]
}
n_tmt_a <- sum(!is.na(cog$cog_tmt_a_raw))
n_from_i1 <- if (has_tmt_a_i1) sum(!is.na(cog[["p20156_i1"]])) else 0L
cat("  TMT-A: n_valid=", n_tmt_a, "(", n_from_i1, "from i1)\n")

# --- TMT-B: i1 preferred, fallback to i0 ---
has_tmt_b_i1 <- "p20157_i1" %in% names(cog)
has_tmt_b_i0 <- "p20157_i0" %in% names(cog)
if (has_tmt_b_i1 && has_tmt_b_i0) {
  cog[, cog_tmt_b_raw := fifelse(!is.na(p20157_i1),
                                as.numeric(p20157_i1),
                                as.numeric(p20157_i0))]
} else if (has_tmt_b_i1) {
  cog[, cog_tmt_b_raw := as.numeric(p20157_i1)]
} else if (has_tmt_b_i0) {
  cog[, cog_tmt_b_raw := as.numeric(p20157_i0)]
} else {
  cog[, cog_tmt_b_raw := NA_real_]
}
n_tmt_b <- sum(!is.na(cog$cog_tmt_b_raw))
cat("  TMT-B: n_valid=", n_tmt_b, "\n")

# --- Symbol Digit: i1 preferred, fallback to i0 ---
has_sd_i1 <- "p20159_i1" %in% names(cog)
has_sd_i0 <- "p20159_i0" %in% names(cog)
if (has_sd_i1 && has_sd_i0) {
  cog[, cog_symbol_digit := fifelse(
    !is.na(p20159_i1),
    suppressWarnings(as.numeric(p20159_i1)),
    suppressWarnings(as.numeric(p20159_i0))
  )]
} else if (has_sd_i1) {
  cog[, cog_symbol_digit := suppressWarnings(as.numeric(p20159_i1))]
} else if (has_sd_i0) {
  cog[, cog_symbol_digit := suppressWarnings(as.numeric(p20159_i0))]
} else {
  cog[, cog_symbol_digit := NA_real_]
}
n_sd <- sum(!is.na(cog$cog_symbol_digit))
cat("  Symbol Digit: n_valid=", n_sd, "\n")

# --- Fluid Intelligence: i1 preferred, fallback to i0 ---
# Character type in data (replace mode); extract numeric
has_fi_i1 <- "p20191_i1" %in% names(cog)
has_fi_i0 <- "p20191_i0" %in% names(cog)
if (has_fi_i1 && has_fi_i0) {
  cog[, cog_fluid_iq := fifelse(
    !is.na(p20191_i1),
    suppressWarnings(as.numeric(p20191_i1)),
    suppressWarnings(as.numeric(p20191_i0))
  )]
} else if (has_fi_i1) {
  cog[, cog_fluid_iq := suppressWarnings(as.numeric(p20191_i1))]
} else if (has_fi_i0) {
  cog[, cog_fluid_iq := suppressWarnings(as.numeric(p20191_i0))]
} else {
  cog[, cog_fluid_iq := NA_real_]
}
n_fi <- sum(!is.na(cog$cog_fluid_iq))
cat("  Fluid Intelligence: n_valid=", n_fi, "\n")

# --- Matrix pattern completion: i1 only ---
if ("p20760_i1" %in% names(cog)) {
  cog[, cog_matrix_reasoning := as.numeric(p20760_i1)]
} else {
  cog[, cog_matrix_reasoning := NA_real_]
}
n_mr <- sum(!is.na(cog$cog_matrix_reasoning))
cat("  Matrix pattern completion (i1 only): n_valid=", n_mr, "\n")

# ============================================================================
# 3. TRANSFORMATIONS
# ============================================================================

cat("\n--- Step 3: Log-transformations (Fawns-Ritchie & Deary 2020) ---\n")

# Pairs Matching: log(errors + 1) — handles zero errors
cog[, cog_pairs_errors := log(cog_pairs_errors_raw + 1)]
cat("  cog_pairs_errors = log(raw + 1)\n")

# TMT-A: log(seconds)
cog[, cog_tmt_a := log(cog_tmt_a_raw)]
# Handle invalid values (zero or negative times → NA)
cog[!is.finite(cog_tmt_a), cog_tmt_a := NA_real_]
cat("  cog_tmt_a = log(seconds)\n")

# TMT-B: log(seconds)
cog[, cog_tmt_b := log(cog_tmt_b_raw)]
cog[!is.finite(cog_tmt_b), cog_tmt_b := NA_real_]
cat("  cog_tmt_b = log(seconds)\n")

# TMT B-A: executive function index (Ankudowich et al. 2019)
# Computed on raw seconds, THEN log-transformed
cog[, cog_tmt_ba := {
  diff <- cog_tmt_b_raw - cog_tmt_a_raw
  # B-A should be positive (B is harder); negative values → data quality issue
  fifelse(diff > 0, log(diff), NA_real_)
}]
n_ba <- sum(!is.na(cog$cog_tmt_ba))
cat("  cog_tmt_ba = log(TMT-B_raw - TMT-A_raw): n_valid=", n_ba, "\n")

# Symbol Digit and Fluid IQ: no transformation (count data, already ~normal)
# Matrix pattern completion: no transformation (count of puzzles correct)

# ============================================================================
# 4. OUTLIER TRUNCATION
# ============================================================================
#
# Per Ankudowich et al. (2019): truncate extreme values
# TMT-A/B: cap at 5th/95th percentile (after log)
# Symbol Digit: cap at [1, 40] (per UKB documentation)
# ============================================================================

cat("\n--- Step 4: Outlier truncation ---\n")

# Symbol Digit: cap at 1-40 per UKB documentation
n_sd_trunc <- sum(cog$cog_symbol_digit < 1 | cog$cog_symbol_digit > 40,
                  na.rm = TRUE)
cog[cog_symbol_digit < 1, cog_symbol_digit := 1]
cog[cog_symbol_digit > 40, cog_symbol_digit := 40]
cat("  Symbol Digit truncated to [1, 40]:", n_sd_trunc, "values\n")

# TMT: no hard truncation (log-transform already compresses extremes)
# But flag implausible values (< 5 seconds or > 600 seconds raw)
n_tmt_implausible <- sum(cog$cog_tmt_a_raw < 5 | cog$cog_tmt_a_raw > 600,
                         na.rm = TRUE)
cat("  TMT-A implausible (<5s or >600s):", n_tmt_implausible, "values\n")

# ============================================================================
# 5. ASSEMBLE FINAL EXPOSURE MATRIX
# ============================================================================

cat("\n--- Step 5: Assembly ---\n")

cog_analysis_vars <- c(
  "cog_tmt_a",          # log(TMT-A seconds)
  "cog_tmt_b",          # log(TMT-B seconds)
  "cog_tmt_ba",         # log(TMT-B - TMT-A)
  "cog_symbol_digit",   # raw count correct
  "cog_pairs_errors",   # log(errors + 1)
  "cog_fluid_iq",       # raw count correct (0-13)
  "cog_matrix_reasoning"  # raw count of puzzles correct
)

cog_out <- cog[, c("eid", cog_analysis_vars), with = FALSE]

cat("  Output:", nrow(cog_out), "rows x", length(cog_analysis_vars), "variables\n")

cog_respondents <- cog_out[, .(has_data = any(!is.na(.SD))),
                         .SDcols = cog_analysis_vars, by = eid]
cat("  Respondents:", format(sum(cog_respondents$has_data), big.mark = ","), "\n")

# ============================================================================
# 6. PER-GROUP OUTPUT
# ============================================================================

cat("\n--- Step 6: Per-group output ---\n")

cog_var_dict <- rbindlist(list(
  dict_entry("cog_tmt_a", 20156, "Cognitive function", "Cat 116",
             "continuous",
             "log(seconds); i1 preferred, i0 fallback; lower=better",
             "Processing speed; shares attentional resources with gustatory processing (Fawns-Ritchie & Deary 2020)"),
  dict_entry("cog_tmt_b", 20157, "Cognitive function", "Cat 116",
             "continuous",
             "log(seconds); i1 preferred, i0 fallback; lower=better",
             "Processing speed + executive function; cognitive flexibility for sensory integration"),
  dict_entry("cog_tmt_ba", "20157-20156", "Cognitive function", "Cat 116",
             "derived_continuous",
             "log(TMT-B_raw - TMT-A_raw); NA if B-A ≤ 0",
             "Executive function index (Ankudowich et al. 2019); set-shifting capacity"),
  dict_entry("cog_symbol_digit", 20159, "Cognitive function", "Cat 116",
             "continuous",
             "Number correct [1-40]; i1 preferred, i0 fallback; higher=better",
             "Complex processing speed; highest loading on speed factor (Fawns-Ritchie 2020, loading=0.80)"),
  dict_entry("cog_pairs_errors", 20132, "Cognitive function", "Cat 116",
             "continuous",
             "log(sum of errors across rounds + 1); i0 only; lower=better",
             "Visual episodic memory; lowest test-retest ICC=0.17 (Hagenaars 2016)"),
  dict_entry("cog_fluid_iq", 20191, "Cognitive function", "Cat 116",
             "continuous",
             "Number correct (0-13); i1 preferred, i0 fallback; higher=better",
             "Verbal-numerical reasoning; Cronbach alpha=0.62 (Hagenaars 2016)"),
  dict_entry("cog_matrix_reasoning", 20760, "Cognitive function", "Cat 116",
             "continuous",
             "Number of puzzles correctly solved; i1 only; higher=better",
             "Non-verbal reasoning (matrix pattern completion)")
))

fwrite(cog_var_dict, file.path(OUTPUT_DIR, "followup_exwas_cognitive_function_variable_dict.csv"))
cat("  Saved variable dictionary:", nrow(cog_var_dict), "variables\n")

for (gname in names(pheno_list)) {
  cat("  --- Group:", toupper(gname), "---\n")
  gpheno <- pheno_list[[gname]]
  gdata <- merge(cog_out, gpheno[, .(eid, taste_2w_strict)], by = "eid")
  gdata <- gdata[!is.na(taste_2w_strict)]
  n_cases <- sum(gdata$taste_2w_strict == 1, na.rm = TRUE)
  cat("    N=", nrow(gdata), " cases=", n_cases, "\n")
  
  out_file <- file.path(OUTPUT_DIR,
                        paste0("followup_exwas_cognitive_function_", gname, ".csv"))
  fwrite(gdata[, !"taste_2w_strict"], out_file)
  cat("    Saved:", out_file, "\n")
  
  miss_g <- missingness_report(gdata, cog_analysis_vars)
  fwrite(miss_g, file.path(OUTPUT_DIR,
                           paste0("followup_exwas_cognitive_function_", gname, "_missingness.csv")))
}

cat("\n  Cognitive-function cleaning complete.\n")
