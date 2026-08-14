#!/usr/bin/env Rscript
# =============================================================================
# Disease-wide association scan — per-cohort PheCode case/control matrices
# =============================================================================
#
# fo_long -> temporal filter -> per-cohort PheCode case/control matrix.
# No phenotype or covariate is merged in here: the matrices stay separate and
# the regression engine joins them to the phenotype at run time, so no
# pre-merged subject-level file is ever written.
#
# PheCode mapping uses PheCode Map v1.2 via the PheWAS package createPhenotypes:
#   - min.code.count = 1        single-occurrence inclusion threshold
#   - add.phecode.exclusions    PheCode-specific control exclusion ranges
#   - full.population.ids       defines the control population per cohort
# All resulting PheCodes, including rollup parent codes, are retained: each
# code's rollup or leaf status comes from the PheCode 1.2 definitions, so
# parent-level signals can be collapsed or kept at the reporting stage rather
# than being decided here.
#
# Inputs (under $PROJECT_DIR):
#   output/dwas/derive/fo_long.rds         (from 01)
#   input/analysis_ready/tastetime_assesstime.csv    eid + p53_i0 (baseline
#                                                    assessment date) + p28755
#                                                    (taste-assessment date)
#   input/analysis_ready/phenotype_group{1,2,3}.csv  cohort eid universe
# Outputs (under output/dwas/derive/):
#   phecode_matrix_group{1,2,3}.rds   eid + one logical column per PheCode
#   phecode_summary_group{1,2,3}.csv  per-PheCode case/control/excluded counts
#   time_anchors.csv
#
# Requires the PheWAS package (PheCode Map v1.2); see the README.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(PheWAS)
})

# =============================================================================
# Configuration
# =============================================================================
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
DERIVE_DIR  <- Sys.getenv("DERIVE_DIR",
                          unset = file.path(PROJECT_DIR, "output", "dwas", "derive"))
PHENO_DIR   <- Sys.getenv("PHENO_DIR", unset = file.path(PROJECT_DIR, "input", "analysis_ready"))

FO_LONG_FILE <- file.path(DERIVE_DIR, "fo_long.rds")
TIME_FILE    <- file.path(PHENO_DIR, "tastetime_assesstime.csv")
PHENO_FILES  <- list(
  group1 = file.path(PHENO_DIR, "phenotype_group1.csv"),
  group2 = file.path(PHENO_DIR, "phenotype_group2.csv"),
  group3 = file.path(PHENO_DIR, "phenotype_group3.csv")
)

MIN_CASES      <- 100   # primary: reported by dwas_regression.R downstream
MIN_CASES_SENS <- 50    # sensitivity / held-out validation

dir.create(DERIVE_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Step 1: load fo_long + time anchors
# =============================================================================
cat("=== Step 1: load ===\n")
fo_long <- readRDS(FO_LONG_FILE)
cat(sprintf("  fo_long: %s records, %s people\n",
            format(nrow(fo_long), big.mark = ","), format(uniqueN(fo_long$eid), big.mark = ",")))

time_dt <- fread(TIME_FILE)
setnames(time_dt, "p53_i0", "assess_date", skip_absent = TRUE)
setnames(time_dt, "p28755", "taste_date",  skip_absent = TRUE)
time_dt <- time_dt[, .(eid, assess_date, taste_date)]
time_dt[, assess_date := as.IDate(assess_date)]
time_dt[, taste_date  := as.IDate(taste_date)]

# =============================================================================
# Step 2: collect per-cohort eids
# =============================================================================
cat("\n=== Step 2: cohort eids ===\n")
group_eids <- list()
for (gname in names(PHENO_FILES)) {
  group_eids[[gname]] <- fread(PHENO_FILES[[gname]], select = "eid")$eid
  cat(sprintf("  %s: %s people\n", gname, format(length(group_eids[[gname]]), big.mark = ",")))
}
all_taste_eids <- unlist(group_eids)

# =============================================================================
# Step 3: temporal filter — keep diagnoses before the taste-assessment date
#   (records with an unknown date are retained as "time unknown")
# =============================================================================
cat("\n=== Step 3: temporal filter ===\n")
fo_long <- fo_long[eid %in% all_taste_eids]
fo_long <- merge(fo_long, time_dt, by = "eid", all.x = TRUE)
n_before <- nrow(fo_long)
fo_primary <- fo_long[is.na(first_date) | first_date < taste_date]
cat(sprintf("  %s -> %s records after filter\n",
            format(n_before, big.mark = ","), format(nrow(fo_primary), big.mark = ",")))
rm(fo_long); gc()

# =============================================================================
# Step 4: build the PheCode matrix per cohort
# =============================================================================
cat("\n=== Step 4: build PheCode matrices ===\n")
for (gname in names(group_eids)) {
  cat(sprintf("\n--- %s ---\n", gname))
  eids <- group_eids[[gname]]

  group_input <- data.frame(
    id            = fo_primary$eid[fo_primary$eid %in% eids],
    vocabulary_id = "ICD10CM",
    code          = fo_primary$icd10[fo_primary$eid %in% eids],
    count         = 1L
  )
  cat(sprintf("  diagnosis records: %s\n", format(nrow(group_input), big.mark = ",")))

  pheno_matrix <- createPhenotypes(
    group_input,
    min.code.count        = 1,
    add.phecode.exclusions = TRUE,
    full.population.ids   = eids
  )

  phecode_cols <- setdiff(names(pheno_matrix), "id")
  case_counts  <- sapply(phecode_cols, function(pc) sum(pheno_matrix[[pc]] == TRUE, na.rm = TRUE))
  cat(sprintf("  matrix: %d people x %d PheCodes\n", nrow(pheno_matrix), length(phecode_cols)))
  cat(sprintf("  case >= %d: %d PheCodes | case >= %d: %d (sensitivity)\n",
              MIN_CASES, sum(case_counts >= MIN_CASES),
              MIN_CASES_SENS, sum(case_counts >= MIN_CASES_SENS)))

  saveRDS(pheno_matrix, file.path(DERIVE_DIR, sprintf("phecode_matrix_%s.rds", gname)))

  summary_dt <- data.table(
    phecode    = names(case_counts),
    n_case     = as.integer(case_counts),
    n_control  = sapply(names(case_counts), function(pc) sum(pheno_matrix[[pc]] == FALSE, na.rm = TRUE)),
    n_excluded = sapply(names(case_counts), function(pc) sum(is.na(pheno_matrix[[pc]])))
  )
  if ("pheinfo" %in% ls("package:PheWAS")) {
    desc <- as.data.table(PheWAS::pheinfo)
    desc_cols <- intersect(c("phecode", "description", "group", "groupnum"), names(desc))
    if (length(desc_cols) > 1) summary_dt <- merge(summary_dt, desc[, ..desc_cols], by = "phecode", all.x = TRUE)
  }
  setorder(summary_dt, -n_case)
  fwrite(summary_dt, file.path(DERIVE_DIR, sprintf("phecode_summary_%s.csv", gname)))
  cat(sprintf("  saved phecode_matrix_%s.rds + phecode_summary_%s.csv\n", gname, gname))
}

fwrite(time_dt[eid %in% all_taste_eids], file.path(DERIVE_DIR, "time_anchors.csv"))
# Export and dx upload to RAP  (the three PheCode matrices, their summaries and
# the time anchors)
cat("\n=== Done ===\n")
cat("Copy phecode_matrix_group{1,2,3}.rds into input/analysis_ready/, or point\n")
cat("PHECODE_DIR at this directory, before running 04_run_dwas.R.\n")
