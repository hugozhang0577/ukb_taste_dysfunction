#!/usr/bin/env Rscript
# ============================================================================
# Follow-up ExWAS — assemble per-cohort analysis-ready matrices + dictionary
# ============================================================================
#
# No variable is dropped for statistical reasons here: inclusion thresholds
# (n_cases >= 100 and the rest) belong to the regression runner, which sees the
# case counts this script computes. The one flag written here is the
# item-level target-leakage exclusion, which is a definition of the exposure
# set rather than a filter on the data (Step 4).
#
# The script aggregates the per-family files, computes per-variable case
# counts, and writes a readable summary report backing the Methods counts.
#
# Run AFTER 01 to 06. This is a STANDALONE script: it reads everything from
# disk, no R session state required. Configure paths in CONFIGURATION below.
#
# Inputs (in OUTPUT_DIR1):
#   followup_exwas_<family>_{group1,group2,group3}.csv
#   followup_exwas_<family>_{group1,group2,group3}_missingness.csv
#   followup_exwas_<family>_variable_dict.csv
#
# Inputs (in PHENO_DIR):
#   phenotype_group1.csv
#   phenotype_group2.csv
#   phenotype_group3.csv
#
# Outputs (in OUTPUT_DIR2):
#   exwas_followup_{group}.csv               - merged exposure matrix
#   followup_exwas_variable_dictionary.csv        - unified dictionary
#   followup_exwas_combined_missingness.csv       - aggregated missingness
#   followup_exwas_variable_case_counts.csv       - per-var per-group N counts
#   followup_exwas_source_overlap.csv             - questionnaire completion patterns
#   followup_exwas_summary_report.txt             - human-readable report
#
# NOTE: analysis_ready CSVs contain ONLY exposures + eid (NO phenotype/
# covariates). The downstream regression script merges phenotype + exposures
# at runtime. This keeps phenotype as a single source of truth across the
# protein-, metabolite-, exposure-, and disease-wide scans.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================================
# UTILITY: Count non-NA per row (robust to any column count)
# ============================================================================
# Avoids rowSums/as.matrix pitfalls that fail on single-column data.tables
# in some R/data.table versions.

n_nonNA_per_row <- function(dt) {
  if (is.null(dt) || ncol(dt) == 0) return(integer(0))
  Reduce("+", lapply(dt, function(x) as.integer(!is.na(x))))
}

# ============================================================================
# CONFIGURATION
# ============================================================================

# Paths resolved under $PROJECT_DIR (override individual vars as needed).
#   OUTPUT_DIR1 : per-section files written by the cleaning scripts (READ)
#   OUTPUT_DIR2 : analysis-ready matrices + dictionary (engine exposure inputs;
#                 point the regression engine's EXPOSURE_DIR here)
#   PHENO_DIR   : per-cohort phenotype files (READ)
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
OUTPUT_DIR1 <- Sys.getenv("DERIVE_PER_SECTION_DIR",
                          unset = file.path(PROJECT_DIR, "output",
                                            "followup_exwas", "derive", "per_section"))
OUTPUT_DIR2 <- Sys.getenv("EXPOSURE_DIR",
                          unset = file.path(PROJECT_DIR, "output",
                                            "followup_exwas", "derive"))
PHENO_DIR   <- Sys.getenv("PHENO_DIR", unset = file.path(PROJECT_DIR, "input", "analysis_ready"))

# Phenotype file naming convention: {prefix}{group}.csv
# Files expected: phenotype_group1.csv, phenotype_group2.csv, phenotype_group3.csv
PHENO_FILE_PREFIX <- "phenotype_"

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR2)) {
  dir.create(OUTPUT_DIR2, recursive = TRUE)
  cat("Created output directory:", OUTPUT_DIR2, "\n")
}

# Exposure families, with the UK Biobank online-follow-up category each comes
# from. The family name is the value stored in the dictionary and carried into
# every results table, so nothing downstream needs a code-to-name lookup.
FAMILIES <- c(
  "Food preferences",       # Category 1039
  "Diet (24-hour recall)",  # Categories 100117 (nutrients) and 100112 (supplements)
  "Digestive health",       # Category 153
  "Mental health",          # Category 136, MHQ1
  "Experience of pain",     # Category 154
  "Cognitive function",     # Category 116
  "Work environment"        # Category 130
)

# Family names are display strings, so the per-family file names use a slug of
# the name rather than a second vocabulary that could drift away from it.
family_slug <- function(x) gsub("^_|_$", "", gsub("_+", "_", gsub("[^A-Za-z0-9]+", "_", tolower(x))))

GROUPS <- c("group1", "group2", "group3")
GROUP_LABELS <- c(
  group1 = "Pure White British",
  group2 = "Other White",
  group3 = "Non-White"
)

# Outcome variable in phenotype files
PHENO_VAR <- "taste_2w_strict"

# Meta columns to exclude from variable lists
META_COLS <- c("eid", "nutr_n_valid_instances", "nutr_single_instance")

# ============================================================================
# STEP 0: LOAD PHENOTYPE
# ============================================================================

cat("================================================================\n")
cat("  follow-up ExWAS Post-Processing\n")
cat("================================================================\n")
cat("Data input  (per-section CSVs):", OUTPUT_DIR1, "\n")
cat("Pheno input (phenotype CSVs): ", PHENO_DIR, "\n")
cat("Output      (merged datasets):", OUTPUT_DIR2, "\n")

# Always load phenotype from disk (standalone script — no R session dependency)
# Phenotype is loaded ONLY to compute case counts and respondent statistics.
# Covariates are NOT touched here — they belong to the downstream regression script.
cat("\nLoading phenotype files from:", PHENO_DIR, "\n")
pheno_list <- list()
for (g in GROUPS) {
  pf <- file.path(PHENO_DIR, paste0(PHENO_FILE_PREFIX, g, ".csv"))
  if (file.exists(pf)) {
    pheno_list[[g]] <- fread(pf, select = c("eid", PHENO_VAR))
    cat("  Loaded ", g, ": ", format(nrow(pheno_list[[g]]), big.mark = ","),
        " rows\n", sep = "")
  } else {
    cat("  ERROR: ", pf, " not found\n", sep = "")
  }
}

# Cohort sizes & case counts
cohort_stats <- list()
for (g in GROUPS) {
  if (is.null(pheno_list[[g]])) next
  p <- pheno_list[[g]]
  if (!PHENO_VAR %in% names(p)) {
    cat("  WARNING:", PHENO_VAR, "not in", g, "\n")
    next
  }
  p_with_outcome <- p[!is.na(get(PHENO_VAR))]
  n_total <- nrow(p_with_outcome)
  n_cases <- sum(p_with_outcome[[PHENO_VAR]] == 1)
  n_controls <- sum(p_with_outcome[[PHENO_VAR]] == 0)
  cohort_stats[[g]] <- list(
    n_total = n_total,
    n_cases = n_cases,
    n_controls = n_controls,
    case_rate = n_cases / n_total
  )
}

cat("\nCohort sizes (with phenotype):\n")
for (g in GROUPS) {
  cs <- cohort_stats[[g]]
  if (is.null(cs)) next
  cat(sprintf("  %s: total=%s, cases=%s (%.2f%%), controls=%s\n",
              g,
              format(cs$n_total, big.mark = ","),
              format(cs$n_cases, big.mark = ","),
              100 * cs$case_rate,
              format(cs$n_controls, big.mark = ",")))
}

# ============================================================================
# STEP 1: FILE AUDIT
# ============================================================================

cat("\n--- Step 1: File audit ---\n")

audit <- data.table()
for (src in FAMILIES) {
  for (g in GROUPS) {
    df <- file.path(OUTPUT_DIR1,
                    paste0("followup_exwas_", family_slug(src), "_", g, ".csv"))
    mf <- file.path(OUTPUT_DIR1,
                    paste0("followup_exwas_", family_slug(src), "_", g, "_missingness.csv"))
    audit <- rbind(audit, data.table(
      source = src, group = g,
      data_file = df,
      miss_file = mf,
      data_exists = file.exists(df),
      miss_exists = file.exists(mf)
    ))
  }
}

n_data_missing <- sum(!audit$data_exists)
if (n_data_missing > 0) {
  cat("  WARNING:", n_data_missing, "data files missing:\n")
  for (i in which(!audit$data_exists)) {
    cat("    -", basename(audit$data_file[i]), "\n")
  }
}
cat("  Files found:", sum(audit$data_exists), "/", nrow(audit), "\n")

# ============================================================================
# STEP 2: AGGREGATE PER-SECTION MISSINGNESS (NO RECOMPUTATION)
# ============================================================================

cat("\n--- Step 2: Aggregate per-section missingness ---\n")

miss_list <- list()
for (i in seq_len(nrow(audit))) {
  if (!audit$miss_exists[i]) next
  m <- fread(audit$miss_file[i])
  m[, source := audit$source[i]]
  m[, group  := audit$group[i]]
  miss_list[[paste(audit$source[i], audit$group[i], sep = "_")]] <- m
}

if (length(miss_list) > 0) {
  combined_miss <- rbindlist(miss_list, fill = TRUE)
  fwrite(combined_miss,
         file.path(OUTPUT_DIR2, "followup_exwas_combined_missingness.csv"))
  cat("  Aggregated missingness rows:", nrow(combined_miss), "\n")
  cat("  Saved: followup_exwas_combined_missingness.csv\n")
} else {
  combined_miss <- data.table()
  cat("  WARNING: no missingness files found\n")
}

# ============================================================================
# STEP 3: PER-VARIABLE CASE COUNTS (FOR DOWNSTREAM REGRESSION FILTER)
# ============================================================================
#
# For each variable * each group, compute:
#   - n_complete = participants with non-NA exposure AND non-NA outcome
#   - n_cases_complete = of those, how many are cases
#   - n_controls_complete = of those, how many are controls
#
# These counts drive the downstream regression filter
# (e.g., n_cases_complete >= 100 -> include in ExWAS).
# ============================================================================

cat("\n--- Step 3: Per-variable case counts ---\n")

case_count_rows <- list()
per_section_data <- list()  # cache for Step 5

for (src in FAMILIES) {
  for (g in GROUPS) {
    df <- audit[source == src & group == g, data_file]
    if (length(df) == 0 || !file.exists(df)) next
    if (is.null(pheno_list[[g]])) next
    
    dat <- fread(df)
    var_cols <- setdiff(names(dat), META_COLS)
    if (length(var_cols) == 0) next
    
    # Cache for Step 5 merging
    per_section_data[[paste(src, g, sep = "_")]] <- dat
    
    # Merge with outcome
    pheno <- pheno_list[[g]][!is.na(get(PHENO_VAR)),
                             .(eid, outcome = get(PHENO_VAR))]
    merged <- merge(dat, pheno, by = "eid")
    n_total_g <- nrow(merged)
    
    for (v in var_cols) {
      complete_idx <- !is.na(merged[[v]])
      n_complete <- sum(complete_idx)
      if (n_complete == 0) {
        n_cases_c <- 0L
        n_controls_c <- 0L
        v_var <- NA_real_
      } else {
        sub <- merged[complete_idx]
        n_cases_c <- sum(sub$outcome == 1)
        n_controls_c <- sum(sub$outcome == 0)
        v_var <- if (n_complete >= 2) {
          suppressWarnings(var(merged[[v]], na.rm = TRUE))
        } else NA_real_
      }
      
      case_count_rows[[length(case_count_rows) + 1L]] <- data.table(
        source = src,
        group = g,
        variable = v,
        n_total_group = n_total_g,
        n_complete = n_complete,
        complete_pct = round(100 * n_complete / n_total_g, 2),
        n_cases_complete = n_cases_c,
        n_controls_complete = n_controls_c,
        case_rate_complete = if (n_complete > 0) round(100 * n_cases_c / n_complete, 3) else NA_real_,
        variance = v_var
      )
    }
  }
}

case_counts <- rbindlist(case_count_rows)
fwrite(case_counts,
       file.path(OUTPUT_DIR2, "followup_exwas_variable_case_counts.csv"))
cat("  Per-variable case count rows:", nrow(case_counts), "\n")
cat("  Saved: followup_exwas_variable_case_counts.csv\n")

# ============================================================================
# STEP 4: UNIFIED VARIABLE DICTIONARY (NO FILTERING)
# ============================================================================
#
# No variable is dropped here. The one flag the dictionary does carry is
# exclude_primary, for items whose wording restates the outcome rather than
# predicting it (item-level target leakage). The regression runner drops those
# and reports them; everything else is decided downstream on case counts.
#
# The exclusion is declarative: after this step the exclude_primary state of
# the dictionary is fully determined by EXCLUDED_VARS below, whatever it was
# before.
#
#   gi_life_interference — a digestive-health item whose "life interference"
#     self-rating wording overlaps directly with the functional-impact
#     criterion of the taste outcome. It is also unavailable at clinical
#     deployment, so it is excluded from both the association scan and any
#     prediction model. The remaining items — including the PHQ-15 somatic
#     items, which describe symptom burden rather than its consequences — are
#     retained and join the comorbidity-correlate layer in interpretation.

EXCLUDED_VARS <- c(
  gi_life_interference = "item_level_target_leakage_life_interference_wording"
)

cat("\n--- Step 4: Unified variable dictionary ---\n")

dict_list <- list()
for (src in FAMILIES) {
  dict_file <- file.path(OUTPUT_DIR1,
                         paste0("followup_exwas_", family_slug(src), "_variable_dict.csv"))
  if (!file.exists(dict_file)) {
    cat("  ", src, ": no dictionary file\n", sep = "")
    next
  }
  d <- fread(dict_file)
  dict_list[[src]] <- d
}

if (length(dict_list) > 0) {
  unified_dict <- rbindlist(dict_list, fill = TRUE)

  missing_excl <- setdiff(names(EXCLUDED_VARS), unified_dict$var_name)
  if (length(missing_excl) > 0)
    stop("variables flagged for exclusion are not in the dictionary: ",
         paste(missing_excl, collapse = ", "))

  unified_dict[, exclude_primary := FALSE]
  unified_dict[, exclude_reason  := NA_character_]
  for (v in names(EXCLUDED_VARS))
    unified_dict[var_name == v, `:=`(exclude_primary = TRUE,
                                     exclude_reason  = EXCLUDED_VARS[[v]])]

  fwrite(unified_dict,
         file.path(OUTPUT_DIR2, "followup_exwas_variable_dictionary.csv"))
  cat("  Unified dictionary:", nrow(unified_dict), "variables\n")
  cat("  Flagged exclude_primary:", sum(unified_dict$exclude_primary), "-",
      paste(names(EXCLUDED_VARS), collapse = ", "), "\n")
  cat("  Available to the scan:", sum(!unified_dict$exclude_primary), "\n")
} else {
  unified_dict <- data.table()
}

# ============================================================================
# STEP 5: PER-GROUP MERGING (NO FILTERING - ALL VARIABLES)
# ============================================================================

cat("\n--- Step 5: Per-group merging ---\n")

merge_stats <- list()

for (g in GROUPS) {
  cat("\n  --- Group:", toupper(g), "---\n")
  
  merged <- NULL
  for (src in FAMILIES) {
    key <- paste(src, g, sep = "_")
    dat <- per_section_data[[key]]
    if (is.null(dat)) {
      cat("    ", src, ": no data\n", sep = "")
      next
    }
    
    keep_cols <- setdiff(names(dat), c("nutr_n_valid_instances",
                                       "nutr_single_instance"))
    dat <- dat[, ..keep_cols]
    
    var_cols <- setdiff(keep_cols, "eid")
    n_resp <- sum(n_nonNA_per_row(dat[, ..var_cols]) > 0)
    cat("    ", src, ": ", format(nrow(dat), big.mark = ","),
        " rows, ", length(var_cols), " vars, ",
        format(n_resp, big.mark = ","), " respondents\n", sep = "")
    
    if (is.null(merged)) {
      merged <- dat
    } else {
      merged <- merge(merged, dat, by = "eid", all = TRUE)
    }
  }
  
  if (is.null(merged)) {
    cat("    No data merged for", g, "\n")
    next
  }
  
  out_file <- file.path(OUTPUT_DIR2,
                        paste0("exwas_followup_", g, ".csv"))
  fwrite(merged, out_file)
  
  n_vars_total <- ncol(merged) - 1L
  any_data_idx <- n_nonNA_per_row(merged[, -"eid"]) > 0
  n_any <- sum(any_data_idx)
  
  cat("    MERGED: ", format(nrow(merged), big.mark = ","),
      " participants x ", n_vars_total, " variables\n", sep = "")
  cat("    Participants with >=1 follow-up ExWAS variable: ",
      format(n_any, big.mark = ","), "\n", sep = "")
  
  merge_stats[[g]] <- list(
    n_total_eids = nrow(merged),
    n_any_exposure_data = n_any,
    n_total_vars = n_vars_total
  )
}

# ============================================================================
# STEP 6: SOURCE OVERLAP DIAGNOSTIC
# ============================================================================

cat("\n--- Step 6: Source overlap ---\n")

overlap_list <- list()

for (g in GROUPS) {
  overlap_dt <- NULL
  for (src in FAMILIES) {
    key <- paste(src, g, sep = "_")
    dat <- per_section_data[[key]]
    if (is.null(dat)) next
    
    var_cols <- setdiff(names(dat), META_COLS)
    if (length(var_cols) == 0) next
    
    has_data_idx <- n_nonNA_per_row(dat[, ..var_cols]) > 0
    src_resp <- dat[has_data_idx, .(eid)]
    src_resp[, (src) := 1L]
    
    if (is.null(overlap_dt)) {
      overlap_dt <- src_resp
    } else {
      overlap_dt <- merge(overlap_dt, src_resp, by = "eid", all = TRUE)
    }
  }
  
  if (is.null(overlap_dt)) next
  
  src_cols <- intersect(FAMILIES, names(overlap_dt))
  if (length(src_cols) == 0) next
  
  # Replace NA with 0 for each source flag
  for (sc in src_cols) {
    overlap_dt[is.na(get(sc)), (sc) := 0L]
  }
  
  # Compute n_sources by iterative addition (avoids as.matrix on .SD)
  overlap_dt[, n_sources := 0L]
  for (sc in src_cols) {
    overlap_dt[, n_sources := n_sources + as.integer(get(sc))]
  }
  
  dist <- overlap_dt[, .N, by = n_sources][order(n_sources)]
  dist[, group := g]
  overlap_list[[g]] <- dist
}

if (length(overlap_list) > 0) {
  source_overlap <- rbindlist(overlap_list)
  fwrite(source_overlap,
         file.path(OUTPUT_DIR2, "followup_exwas_source_overlap.csv"))
  cat("  Saved: followup_exwas_source_overlap.csv\n")
} else {
  source_overlap <- data.table()
}

# ============================================================================
# STEP 7: HUMAN-READABLE SUMMARY REPORT
# ============================================================================

cat("\n--- Step 7: Generate summary report ---\n")

report_file <- file.path(OUTPUT_DIR2, "followup_exwas_summary_report.txt")
con <- file(report_file, "w")

wln <- function(...) cat(..., "\n", file = con, sep = "")

wln("================================================================")
wln("  Follow-up ExWAS Post-Processing Summary Report")
wln("  Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
wln("================================================================")
wln("")
wln("Philosophy: NO variable filtering applied at this stage.")
wln("Downstream regression script applies case-count threshold")
wln("(e.g., n_cases_complete >= 100) per variable x group.")
wln("")

# ----- Section 1: Cohort -----
wln("[1] COHORT SIZES (with valid taste_2w_strict)")
wln("--------------------------------------------------------")
wln(sprintf("%-10s %-22s %12s %12s %10s",
            "Group", "Label", "Total", "Cases", "Case%"))
for (g in GROUPS) {
  cs <- cohort_stats[[g]]
  if (is.null(cs)) next
  wln(sprintf("%-10s %-22s %12s %12s %9.2f%%",
              g, GROUP_LABELS[g],
              format(cs$n_total, big.mark = ","),
              format(cs$n_cases, big.mark = ","),
              100 * cs$case_rate))
}
wln("")

# ----- Section 2: Per-source x per-group -----
wln("[2] PER-SOURCE x PER-GROUP SUMMARY")
wln("--------------------------------------------------------")

for (src in FAMILIES) {
  wln("")
  wln("  ", src)
  wln("  ", strrep("-", 56))
  wln(sprintf("  %-10s %12s %12s %12s %8s",
              "Group", "Respondents", "Resp_w_pheno", "Cases_resp", "N_vars"))
  
  for (g in GROUPS) {
    key <- paste(src, g, sep = "_")
    dat <- per_section_data[[key]]
    if (is.null(dat)) {
      wln(sprintf("  %-10s %12s %12s %12s %8s",
                  g, "-", "-", "-", "-"))
      next
    }
    
    var_cols <- setdiff(names(dat), META_COLS)
    n_resp <- sum(n_nonNA_per_row(dat[, ..var_cols]) > 0)
    
    if (!is.null(pheno_list[[g]])) {
      pheno <- pheno_list[[g]][!is.na(get(PHENO_VAR)),
                               .(eid, outcome = get(PHENO_VAR))]
      has_data_idx2 <- n_nonNA_per_row(dat[, ..var_cols]) > 0
      resp_eids <- dat[has_data_idx2, .(eid)]
      m <- merge(resp_eids, pheno, by = "eid")
      n_resp_pheno <- nrow(m)
      n_cases_resp <- sum(m$outcome == 1)
    } else {
      n_resp_pheno <- NA_integer_
      n_cases_resp <- NA_integer_
    }
    
    wln(sprintf("  %-10s %12s %12s %12s %8d",
                g,
                format(n_resp, big.mark = ","),
                format(n_resp_pheno, big.mark = ","),
                format(n_cases_resp, big.mark = ","),
                length(var_cols)))
  }
}
wln("")

# ----- Section 3: Source overlap -----
wln("")
wln("[3] SOURCE OVERLAP (number of questionnaires completed)")
wln("--------------------------------------------------------")

if (nrow(source_overlap) > 0) {
  for (g in GROUPS) {
    sub <- source_overlap[group == g]
    if (nrow(sub) == 0) next
    wln("")
    wln("  ", g, " (", GROUP_LABELS[g], "):")
    
    denom <- if (!is.null(cohort_stats[[g]])) cohort_stats[[g]]$n_total else sum(sub$N)
    
    for (i in seq_len(nrow(sub))) {
      ns <- sub$n_sources[i]
      n  <- sub$N[i]
      pct <- 100 * n / denom
      wln(sprintf("    %d source(s): %12s  (%5.2f%% of cohort)",
                  ns, format(n, big.mark = ","), pct))
    }
    
    n_any_in_overlap <- sum(sub[n_sources > 0, N])
    wln(sprintf("    Any source : %12s  (%5.2f%% of cohort)",
                format(n_any_in_overlap, big.mark = ","),
                100 * n_any_in_overlap / denom))
  }
}
wln("")

# ----- Section 4: Final merged datasets -----
wln("")
wln("[4] FINAL MERGED ANALYSIS-READY DATASETS")
wln("--------------------------------------------------------")
wln(sprintf("  %-10s %15s %15s %12s",
            "Group", "Total_eids", "With_any_exp", "N_variables"))
for (g in GROUPS) {
  ms <- merge_stats[[g]]
  if (is.null(ms)) next
  wln(sprintf("  %-10s %15s %15s %12d",
              g,
              format(ms$n_total_eids, big.mark = ","),
              format(ms$n_any_exposure_data, big.mark = ","),
              ms$n_total_vars))
}
wln("")

# ----- Section 5: Variable case-count distribution -----
wln("")
wln("[5] VARIABLE CASE-COUNT DISTRIBUTION (group1)")
wln("--------------------------------------------------------")
wln("    How many variables have >=N cases with non-missing exposure?")
wln("    Useful for choosing the regression filter threshold.")
wln("")

if (nrow(case_counts) > 0) {
  cc_g1 <- case_counts[group == "group1"]
  thresholds <- c(20, 50, 100, 200, 500, 1000)
  
  wln(sprintf("  %-12s %12s %12s",
              "Threshold", "N_vars", "% of total"))
  total_vars_g1 <- nrow(cc_g1)
  for (th in thresholds) {
    n_pass <- sum(cc_g1$n_cases_complete >= th)
    wln(sprintf("  cases>=%-5d %12d %11.1f%%",
                th, n_pass, 100 * n_pass / total_vars_g1))
  }
  
  wln("")
  wln("  Per-source breakdown at cases>=100 threshold (group1):")
  for (src in FAMILIES) {
    sub <- cc_g1[source == src]
    if (nrow(sub) == 0) next
    n_pass <- sum(sub$n_cases_complete >= 100)
    wln(sprintf("    %s: %d / %d variables pass",
                src, n_pass, nrow(sub)))
  }
}
wln("")

# ----- Section 6: Top variables by case count -----
wln("")
wln("[6] TOP 10 VARIABLES BY CASE COUNT (group1)")
wln("--------------------------------------------------------")
if (nrow(case_counts) > 0) {
  top10 <- case_counts[group == "group1"][order(-n_cases_complete)][1:10]
  wln(sprintf("  %-30s %-5s %12s %12s %10s",
              "Variable", "Src", "N_complete", "N_cases", "Case%"))
  for (i in seq_len(nrow(top10))) {
    wln(sprintf("  %-30s %-5s %12s %12s %9.3f%%",
                substr(top10$variable[i], 1, 30),
                top10$source[i],
                format(top10$n_complete[i], big.mark = ","),
                format(top10$n_cases_complete[i], big.mark = ","),
                top10$case_rate_complete[i]))
  }
}
wln("")

# ----- Section 7: Bottom variables by case count -----
wln("")
wln("[7] BOTTOM 10 VARIABLES BY CASE COUNT (group1, n_complete>0)")
wln("--------------------------------------------------------")
wln("    Candidates for downstream filtering")
if (nrow(case_counts) > 0) {
  bot10 <- case_counts[group == "group1" & n_complete > 0][order(n_cases_complete)][1:10]
  wln(sprintf("  %-30s %-5s %12s %12s %10s",
              "Variable", "Src", "N_complete", "N_cases", "Case%"))
  for (i in seq_len(nrow(bot10))) {
    wln(sprintf("  %-30s %-5s %12s %12s %9.3f%%",
                substr(bot10$variable[i], 1, 30),
                bot10$source[i],
                format(bot10$n_complete[i], big.mark = ","),
                format(bot10$n_cases_complete[i], big.mark = ","),
                bot10$case_rate_complete[i]))
  }
}
wln("")

# ----- Section 8: Files generated -----
wln("")
wln("[8] FILES GENERATED")
wln("--------------------------------------------------------")
wln("  Merged datasets:")
for (g in GROUPS) {
  wln("    exwas_followup_", g, ".csv")
}
wln("  Variable dictionary:")
wln("    followup_exwas_variable_dictionary.csv")
wln("  Diagnostics:")
wln("    followup_exwas_combined_missingness.csv")
wln("    followup_exwas_variable_case_counts.csv")
wln("    followup_exwas_source_overlap.csv")
wln("    followup_exwas_summary_report.txt  (this file)")
wln("")
wln("================================================================")
wln("  END OF REPORT")
wln("================================================================")

close(con)

cat("  Saved:", basename(report_file), "\n")
# Export and dx upload to RAP  (the per-cohort exposure matrices, the unified
# dictionary, the diagnostics and the summary report)

# ============================================================================
# DONE
# ============================================================================

cat("\n================================================================\n")
cat("  follow-up ExWAS post-processing complete.\n")
cat("================================================================\n")
cat("\nAll outputs written to:", OUTPUT_DIR2, "\n")
cat("\nfollowup_exwas_summary_report.txt carries the per-source variable and\n")
cat("case counts that the Methods section reports.\n")

