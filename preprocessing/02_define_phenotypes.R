# =============================================================================
# Phenotype and chemosensory-feature definitions
# =============================================================================
# Derives, from the same seven questionnaire columns of the base table:
#   (A) the taste phenotypes, including the primary outcome taste_2w_strict
#   (B) the three machine-learning smell features
#
# Input
#   input/analysis_ready/base_table_full.csv, built by 01_build_base_table.R,
#   which renames the questionnaire fields to the canonical columns used here:
#     p28615 -> taste_change   p28612 -> smell_change
#     p28616 -> taste_time     p28613 -> smell_time
#     p28617 -> taste_extent   p28614 -> smell_extent
#
# Outputs
#   output/gwas_phenotypes/taste_gwas_phenotypes.csv         -> 04_gwas_sample_qc.R
#   output/gwas_phenotypes/taste_gwas_phenotype_summary.csv
#   output/smell_features/smell_features.csv                 -> the ML feature set
#   output/smell_features/smell_feature_summary.csv
#
# The two sections deliberately work on separate copies of the input, because
# they treat a negative smell_change code differently: the taste definitions
# read it as "no smell change" (so the participant can still serve as a control),
# whereas the smell features read it as missing (an unknown smell status must not
# become a smell_any of 0). Sharing one mutated column would corrupt one of them.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

project_dir <- Sys.getenv("PROJECT_DIR", unset = ".")
input_dir <- file.path(project_dir, "input", "analysis_ready")
pheno_dir <- file.path(project_dir, "output", "gwas_phenotypes")
smell_dir <- file.path(project_dir, "output", "smell_features")
dir.create(pheno_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(smell_dir, recursive = TRUE, showWarnings = FALSE)

normalise_negative_to_na <- function(x) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  x[x < 0] <- NA_real_
  x
}

base_cols <- c("eid", "taste_change", "taste_time", "taste_extent",
               "smell_change", "smell_time", "smell_extent")
base_tbl <- fread(file.path(input_dir, "base_table_full.csv"), select = base_cols)

missing <- setdiff(base_cols, names(base_tbl))
if (length(missing) > 0) {
  stop("Missing required columns in the base table: ", paste(missing, collapse = ", "))
}


# =============================================================================
# (A) Taste phenotypes
# =============================================================================

# Basic Cleaning ==========
# Definition logic:
# Negative duration or impact codes are treated as missing. A negative
# smell_change code is set to 0, so that a participant who reported no taste
# change but gave no usable smell answer still qualifies as a strict control.

taste <- copy(base_tbl)
taste[, `:=`(
  taste_change = suppressWarnings(as.numeric(as.character(taste_change))),
  smell_change = suppressWarnings(as.numeric(as.character(smell_change))),
  taste_time_clean = normalise_negative_to_na(taste_time),
  taste_extent_clean = normalise_negative_to_na(taste_extent)
)]

taste[smell_change < 0, smell_change := 0]
taste[, only_smell := as.integer(taste_change == 0 & smell_change == 1)]


# Basic Taste Phenotype ==========
# Definition logic:
# taste_basic_strict:
#   case = any reported taste change
#   control = no taste change and no smell change
#   isolated smell change = excluded
# taste_basic_relaxed:
#   case = any reported taste change
#   control = no taste change, regardless of smell-change status

taste[, taste_basic_strict := fcase(
  taste_change == 1, 1,
  taste_change == 0 & smell_change == 0, 0,
  default = NA_real_
)]

taste[, taste_basic_relaxed := fcase(
  taste_change == 1, 1,
  taste_change == 0, 0,
  default = NA_real_
)]


# Two-week or Impact Taste Phenotype ==========
# Definition logic:
# taste_2w_strict is the primary taste phenotype used in the manuscript:
#   case = taste change lasting at least 2 weeks, or taste change affecting
#          daily life — either criterion suffices, not both
#   control = no taste change and no smell change
#   taste change not meeting the duration/impact rule = excluded
#   isolated smell change = excluded
#
# Source coding of the questionnaire items:
#   taste_time_clean >= 1 corresponds to at least 2 weeks
#   taste_extent_clean == 1 corresponds to daily-life impact

taste[, has_time_extent_info := fcase(
  taste_change == 0, TRUE,
  taste_change == 1, !is.na(taste_time_clean) | !is.na(taste_extent_clean),
  default = FALSE
)]

taste[, meets_2w_extent := taste_change == 1 &
        ((!is.na(taste_time_clean) & taste_time_clean >= 1) |
           (!is.na(taste_extent_clean) & taste_extent_clean == 1))]

taste[, taste_2w_strict := fcase(
  taste_change == 0 & smell_change == 0, 0,
  meets_2w_extent == TRUE, 1,
  taste_change == 1 & has_time_extent_info == TRUE, NA_real_,
  default = NA_real_
)]

taste[, taste_2w_relaxed := fcase(
  taste_change == 0, 0,
  meets_2w_extent == TRUE, 1,
  taste_change == 1 & has_time_extent_info == TRUE, NA_real_,
  default = NA_real_
)]


# Four-week or Impact Sensitivity Phenotype ==========
# Definition logic:
# A sensitivity phenotype, not the primary taste definition. It uses:
#   taste_time_clean >= 2, corresponding to at least 4 weeks
#   or taste_extent_clean == 1, corresponding to daily-life impact.

taste[, meets_4w_extent := taste_change == 1 &
        ((!is.na(taste_time_clean) & taste_time_clean >= 2) |
           (!is.na(taste_extent_clean) & taste_extent_clean == 1))]

taste[, taste_4w_strict := fcase(
  taste_change == 0 & smell_change == 0, 0,
  meets_4w_extent == TRUE, 1,
  taste_change == 1 & has_time_extent_info == TRUE, NA_real_,
  default = NA_real_
)]


# Taste Output Tables ==========

phenotype_cols <- c(
  "eid",
  "taste_change",
  "smell_change",
  "taste_time_clean",
  "taste_extent_clean",
  "only_smell",
  "has_time_extent_info",
  "meets_2w_extent",
  "meets_4w_extent",
  "taste_basic_strict",
  "taste_basic_relaxed",
  "taste_2w_strict",
  "taste_2w_relaxed",
  "taste_4w_strict"
)

taste_phenotypes <- taste[, ..phenotype_cols]
fwrite(taste_phenotypes, file.path(pheno_dir, "taste_gwas_phenotypes.csv"))
# Export and dx upload to RAP  (taste_gwas_phenotypes.csv feeds the sample-QC step)

pheno_names <- c(
  "taste_basic_strict",
  "taste_basic_relaxed",
  "taste_2w_strict",
  "taste_2w_relaxed",
  "taste_4w_strict"
)

summary_table <- rbindlist(lapply(pheno_names, function(v) {
  data.table(
    phenotype = v,
    n_total_non_missing = sum(!is.na(taste_phenotypes[[v]])),
    n_control = sum(taste_phenotypes[[v]] == 0, na.rm = TRUE),
    n_case = sum(taste_phenotypes[[v]] == 1, na.rm = TRUE),
    n_missing_or_excluded = sum(is.na(taste_phenotypes[[v]]))
  )
}))

summary_table[, case_fraction := n_case / n_total_non_missing]
fwrite(summary_table, file.path(pheno_dir, "taste_gwas_phenotype_summary.csv"))
print(summary_table)


# =============================================================================
# (B) Smell features for the machine-learning feature set
# =============================================================================

# Smell Value Cleaning ==========
# Definition logic:
# UKB negative response codes -1 and -3 are treated as missing for all three
# smell fields. Duration coding is retained as:
#   0 = <2 weeks
#   1 = 2-4 weeks
#   2 = 4-12 weeks
#   3 = >12 weeks
# For smell_extent:
#   0 = no daily-life impact
#   1 = impacts daily life

smell <- base_tbl[, .(eid, smell_change, smell_time, smell_extent)]
for (v in c("smell_change", "smell_time", "smell_extent")) {
  smell[, (v) := suppressWarnings(as.integer(as.character(get(v))))]
  smell[get(v) %in% c(-1L, -3L), (v) := NA_integer_]
}


# Smell Feature Construction ==========
# Definition logic:
# Final ML features:
#   smell_any: binary indicator derived from p28612 / smell_change.
#   smell_time: duration band derived from p28613.
#   smell_extent: daily-life impact derived from p28614.
# smell_time and smell_extent are conditional follow-up features and are set
# to missing when smell_any is 0 or missing.

smell[, smell_any := fcase(
  smell_change == 1L, 1L,
  smell_change == 0L, 0L,
  default = NA_integer_
)]

smell[is.na(smell_any) | smell_any == 0L,
      `:=`(smell_time = NA_integer_, smell_extent = NA_integer_)]

smell_features <- smell[, .(eid, smell_any, smell_time, smell_extent)]
fwrite(smell_features, file.path(smell_dir, "smell_features.csv"))
# Export and dx upload to RAP  (smell_features.csv feeds the ML feature set and
# the smell sensitivity analyses)


# Smell Feature Summary ==========
# Definition logic:
# Export aggregate distributions for review. No subject-level values are
# printed to the console.

smell_summary <- rbindlist(list(
  smell_features[, .(feature = "smell_any",
                     level = as.character(smell_any),
                     n = .N), by = smell_any][, smell_any := NULL],
  smell_features[, .(feature = "smell_time",
                     level = as.character(smell_time),
                     n = .N), by = smell_time][, smell_time := NULL],
  smell_features[, .(feature = "smell_extent",
                     level = as.character(smell_extent),
                     n = .N), by = smell_extent][, smell_extent := NULL]
), fill = TRUE)

fwrite(smell_summary, file.path(smell_dir, "smell_feature_summary.csv"))
print(smell_summary)
