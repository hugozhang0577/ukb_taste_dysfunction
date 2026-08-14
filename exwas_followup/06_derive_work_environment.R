#!/usr/bin/env Rscript
# ============================================================================
# Follow-up ExWAS — work environment (Category 130)
# ============================================================================
#
# Operates on the in-memory `raw` (+ `pheno_list`) loaded by
# 01_derive_mental_digestive.R.
#
# Input:   `raw` data.table (from shared framework Steps 1-3)
# Output:  followup_exwas_work_environment_{group}.csv
#          followup_exwas_work_environment_variable_dict.csv
#
# Strategy:
#   Cat 130 (Employment History) collected per-job exposure ratings as
#   array fields (a0 through a39 representing up to 40 lifetime jobs).
#   We collapse across all jobs to derive lifetime "ever exposed" binary
#   indicators per Dimakakou et al. (2020), who established the precedent
#   for using UKB self-reported workplace exposure variables.
#
#   Each exposure field uses a Likert scale (Rarely/Never, Sometimes,
#   Often, Always). "Ever exposed" is defined as ANY job rated
#   "Sometimes" or higher across the lifetime work history.
#
# Variables (11 total):
#   Airborne particulates:  work_noisy, work_dusty
#   Chemical exposures:     work_chemical_fumes, work_cigarette_smoke,
#                           work_paints_thinners, work_pesticides, work_diesel
#   Composite:              work_chemical_score (sum of 4 chemical binaries)
#   Respiratory symptoms:   work_breathing_problems
#   Shift work:             work_shift_work_ever, work_night_shifts_ever
#
# Key references:
#   - UKB Work History Questionnaire (2015)
#   - Dimakakou et al. 2020 IJERPH — UKB self-reported workplace exposures
#   - Sadhra et al. 2020 Thorax — UKB occupational AOP analysis
#   - De Matteis et al. 2022 Thorax — UKB lifetime occupational exposures
# ============================================================================

cat("\n================================================================\n")
cat("  Work environment (Category 130)\n")
cat("================================================================\n\n")

# ============================================================================
# 1. FIELD REGISTRY
# ============================================================================

work_exposures <- list(
  noisy            = list(fid = 22606, label = "Workplace had a lot of noise"),
  dusty            = list(fid = 22609, label = "Workplace very dusty"),
  chemical_fumes   = list(fid = 22610, label = "Workplace full of chemical/other fumes"),
  cigarette_smoke  = list(fid = 22611, label = "Workplace had a lot of cigarette smoke"),
  paints_thinners  = list(fid = 22613, label = "Used paints, glues, thinners or strippers"),
  pesticides       = list(fid = 22614, label = "Used pesticides"),
  diesel           = list(fid = 22615, label = "Workplace had a lot of diesel exhaust"),
  breathing_probs  = list(fid = 22616, label = "Job ever caused breathing problems"),
  shift_work       = list(fid = 22620, label = "Job involved shift work"),
  night_shifts     = list(fid = 22650, label = "Job involved night shifts")
)

cat("  Registered exposure fields:", length(work_exposures), "\n")
register_fields("Work environment",
                vapply(work_exposures, function(x) x$fid, numeric(1)))

# ============================================================================
# 2. HELPER: Find all array columns for a field
# ============================================================================

find_array_cols <- function(fid, dt, n_max = 40) {
  # Try standard array naming: p{fid}_a0, p{fid}_a1, ...
  candidate_a <- paste0("p", fid, "_a", 0:(n_max - 1))
  present_a <- intersect(candidate_a, names(dt))
  
  # Some fields may be instance-arrayed: p{fid}_i0_a0, p{fid}_i0_a1
  if (length(present_a) == 0) {
    candidate_ia <- paste0("p", fid, "_i0_a", 0:(n_max - 1))
    present_a <- intersect(candidate_ia, names(dt))
  }
  
  # Single-instance fallback (no array)
  if (length(present_a) == 0) {
    candidate_single <- c(paste0("p", fid), paste0("p", fid, "_i0"))
    present_a <- intersect(candidate_single, names(dt))
  }
  
  present_a
}

# ============================================================================
# 3. HELPER: Ever-exposed collapse
# ============================================================================
#
# Per-job ratings are coded (REPLACE mode) as one of:
#   "Rarely/never", "Sometimes", "Often", "Always",
#   "Prefer not to answer", "Do not know"
#
# Or numeric coding (1=Rarely/Never, 2=Sometimes, 3=Often, 4=Always).
#
# ever_exposed = TRUE if ANY array slot is "Sometimes" or higher,
#                FALSE if ALL array slots are "Rarely/never",
#                NA   if all slots are NA / missing / refused
# ============================================================================

EXPOSED_LEVELS <- c("Sometimes", "Often", "Always", "2", "3", "4")
UNEXPOSED_LEVELS <- c("Rarely/never", "Rarely / never", "Never", "1")
MISSING_LEVELS <- c("", "Prefer not to answer", "Do not know", "NA")

collapse_ever_exposed <- function(dt, cols) {
  if (length(cols) == 0) return(rep(NA_integer_, nrow(dt)))
  
  # Build logical matrices
  is_exposed_mat <- matrix(FALSE, nrow = nrow(dt), ncol = length(cols))
  is_unexposed_mat <- matrix(FALSE, nrow = nrow(dt), ncol = length(cols))
  
  for (j in seq_along(cols)) {
    val <- dt[[cols[j]]]
    if (is.numeric(val)) {
      is_exposed_mat[, j]   <- !is.na(val) & val >= 2
      is_unexposed_mat[, j] <- !is.na(val) & val == 1
    } else {
      val_chr <- as.character(val)
      val_chr[val_chr %in% MISSING_LEVELS] <- NA
      is_exposed_mat[, j]   <- !is.na(val_chr) & val_chr %in% EXPOSED_LEVELS
      is_unexposed_mat[, j] <- !is.na(val_chr) & val_chr %in% UNEXPOSED_LEVELS
    }
  }
  
  any_exposed   <- rowSums(is_exposed_mat)   > 0
  any_unexposed <- rowSums(is_unexposed_mat) > 0
  
  result <- fifelse(
    any_exposed, 1L,
    fifelse(any_unexposed, 0L, NA_integer_)
  )
  result
}

# ============================================================================
# 4. PROCESS EACH EXPOSURE
# ============================================================================

cat("\n--- Step 4: Per-exposure array collapse ---\n")

# Detect available exposures
work <- raw[, .(eid)]
exposures_processed <- character(0)
exposures_skipped <- character(0)

for (exp_name in names(work_exposures)) {
  fid <- work_exposures[[exp_name]]$fid
  array_cols <- find_array_cols(fid, raw)
  
  if (length(array_cols) == 0) {
    cat("  ", exp_name, " (Field ", fid, "): SKIPPED — no columns found\n",
        sep = "")
    exposures_skipped <- c(exposures_skipped, exp_name)
    next
  }
  
  var_name <- paste0("work_", exp_name)
  work[, (var_name) := collapse_ever_exposed(raw, array_cols)]
  
  n_yes <- sum(work[[var_name]] == 1L, na.rm = TRUE)
  n_no  <- sum(work[[var_name]] == 0L, na.rm = TRUE)
  n_na  <- sum(is.na(work[[var_name]]))
  
  cat("  ", var_name, " (Field ", fid, ", ", length(array_cols),
      " arrays): yes=", n_yes, " no=", n_no, " na=", n_na, "\n", sep = "")
  
  exposures_processed <- c(exposures_processed, exp_name)
}

# Adjust expected variable names per processed exposures
work_var_names <- paste0("work_", exposures_processed)

# ============================================================================
# 5. RENAME SHIFT WORK / NIGHT SHIFTS FOR CLARITY
# ============================================================================

# These are 'ever' indicators by construction; rename to make explicit
if ("work_shift_work" %in% names(work)) {
  setnames(work, "work_shift_work", "work_shift_work_ever")
  work_var_names[work_var_names == "work_shift_work"] <- "work_shift_work_ever"
}
if ("work_night_shifts" %in% names(work)) {
  setnames(work, "work_night_shifts", "work_night_shifts_ever")
  work_var_names[work_var_names == "work_night_shifts"] <- "work_night_shifts_ever"
}
if ("work_breathing_probs" %in% names(work)) {
  setnames(work, "work_breathing_probs", "work_breathing_problems")
  work_var_names[work_var_names == "work_breathing_probs"] <- "work_breathing_problems"
}

# ============================================================================
# 6. COMPOSITE CHEMICAL EXPOSURE SCORE
# ============================================================================
#
# Sum of binary chemical exposure indicators (range 0-4):
#   chemical_fumes + paints_thinners + pesticides + diesel
#
# Captures cumulative chemical exposure burden across taste-relevant
# pathways (chemosensory irritation, neurotoxicity).
# ============================================================================

cat("\n--- Step 5: Composite chemical exposure score ---\n")

chemical_components <- c("work_chemical_fumes", "work_paints_thinners",
                         "work_pesticides", "work_diesel")
chemical_present <- intersect(chemical_components, names(work))

if (length(chemical_present) >= 2) {
  # Treat NA as 0 for sum (conservative; participant didn't report exposure)
  # but require at least 1 non-NA value to assign a score
  work[, work_chemical_score := {
    mat <- as.matrix(.SD)
    n_valid <- rowSums(!is.na(mat))
    score <- rowSums(mat, na.rm = TRUE)
    fifelse(n_valid > 0, as.integer(score), NA_integer_)
  }, .SDcols = chemical_present]
  
  cat("  work_chemical_score (sum of ", length(chemical_present),
      " components): ", sep = "")
  cat("range [0,", max(work$work_chemical_score, na.rm = TRUE), "], ",
      "n_valid=", sum(!is.na(work$work_chemical_score)), "\n", sep = "")
  
  work_var_names <- c(work_var_names, "work_chemical_score")
} else {
  cat("  work_chemical_score: SKIPPED — insufficient chemical components\n")
}

# ============================================================================
# 7. ASSEMBLE + PER-GROUP OUTPUT
# ============================================================================

cat("\n--- Step 6: Assembly ---\n")

if (length(work_var_names) == 0) {
  cat("  WARNING: No work-environment variables successfully processed.\n")
  cat("  Cat 130 fields may not be present in this extract.\n")
  cat("  The work-environment family will be skipped — add a Methods note.\n")
} else {
  
  work_out <- work[, c("eid", work_var_names), with = FALSE]
  
  cat("  Output:", nrow(work_out), "rows x",
      length(work_var_names), "variables\n")
  
  work_respondents <- work_out[, .(has_data = any(!is.na(.SD))),
                           .SDcols = work_var_names, by = eid]
  cat("  Respondents:", format(sum(work_respondents$has_data), big.mark = ","), "\n")
  
  # Variable dictionary
  dict_rows <- list()
  for (vn in work_var_names) {
    if (vn == "work_chemical_score") {
      dict_rows[[vn]] <- dict_entry(
        vn, "derived", "Work environment", "Cat 130", "continuous (0-4)",
        paste0("Sum of binary chemical exposures: ",
               paste(chemical_present, collapse = " + ")),
        "Composite chemical exposure burden — chemosensory + neurotoxic taste pathways"
      )
    } else {
      # Recover field id from exposure registry
      exp_key <- sub("^work_", "", vn)
      exp_key <- sub("_ever$", "", exp_key)
      exp_key <- sub("breathing_problems", "breathing_probs", exp_key)
      
      if (exp_key %in% names(work_exposures)) {
        fid_val <- work_exposures[[exp_key]]$fid
        label   <- work_exposures[[exp_key]]$label
      } else {
        fid_val <- NA
        label   <- vn
      }
      dict_rows[[vn]] <- dict_entry(
        vn, fid_val, "Work environment", "Cat 130", "binary",
        paste0("Ever exposed (any job rated 'Sometimes' or higher) across ",
               "lifetime work history (up to 40 jobs)"),
        label
      )
    }
  }
  
  work_var_dict <- rbindlist(dict_rows)
  fwrite(work_var_dict, file.path(OUTPUT_DIR, "followup_exwas_work_environment_variable_dict.csv"))
  cat("  Saved variable dictionary:", nrow(work_var_dict), "variables\n")
  
  for (gname in names(pheno_list)) {
    cat("  --- Group:", toupper(gname), "---\n")
    gpheno <- pheno_list[[gname]]
    gdata <- merge(work_out, gpheno[, .(eid, taste_2w_strict)], by = "eid")
    gdata <- gdata[!is.na(taste_2w_strict)]
    n_cases <- sum(gdata$taste_2w_strict == 1, na.rm = TRUE)
    cat("    N=", nrow(gdata), " cases=", n_cases, "\n")
    
    out_file <- file.path(OUTPUT_DIR,
                          paste0("followup_exwas_work_environment_", gname, ".csv"))
    fwrite(gdata[, !"taste_2w_strict"], out_file)
    cat("    Saved:", out_file, "\n")
    
    miss_g <- missingness_report(gdata, work_var_names)
    fwrite(miss_g, file.path(OUTPUT_DIR,
                             paste0("followup_exwas_work_environment_", gname, "_missingness.csv")))
  }
}

cat("\n  Work-environment cleaning complete.\n")
if (length(exposures_skipped) > 0) {
  cat("  NOTE: Skipped exposures (not in data):",
      paste(exposures_skipped, collapse = ", "), "\n")
  cat("  → Extract from RAP if needed and re-run.\n")
}

# Last derivation step: close the loop on the export specification. Every field
# id in field_list.txt must have been claimed by one of 01 to 06, so an export
# carrying a field no step reads stops the run here rather than passing an
# unused column downstream.
check_field_list_complete()
