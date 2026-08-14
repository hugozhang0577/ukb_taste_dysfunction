#!/usr/bin/env Rscript
# =============================================================================
# Disease-wide association scan — the diagnostic-count covariate
# =============================================================================
#
# Computes each participant's total count of distinct ICD-coded diagnoses
# (total_dx_count) as a health-care-utilisation proxy, and merges it into the
# per-cohort phenotype files. This term is a primary-model covariate in the
# disease scan, included to mitigate ascertainment bias from differential
# health-care contact: participants with more recorded contact accumulate more
# PheCodes for reasons unrelated to taste.
#
# Idempotent: any existing total_dx_count column is dropped and recomputed, so
# re-running leaves the phenotype files in the same state.
#
# Inputs (under $PROJECT_DIR):
#   output/dwas/derive/fo_long.rds          (from 01)
#   input/analysis_ready/phenotype_group{1,2,3}.csv
# Output: total_dx_count written back into each phenotype file.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
DERIVE_DIR  <- Sys.getenv("DERIVE_DIR",
                          unset = file.path(PROJECT_DIR, "output", "dwas", "derive"))
PHENO_DIR   <- Sys.getenv("PHENO_DIR", unset = file.path(PROJECT_DIR, "input", "analysis_ready"))

fo_long  <- readRDS(file.path(DERIVE_DIR, "fo_long.rds"))
dx_count <- fo_long[, .(total_dx_count = uniqueN(icd10)), by = eid]
cat(sprintf("total_dx_count computed for %s people\n",
            format(nrow(dx_count), big.mark = ",")))

GROUPS <- c("group1", "group2", "group3")
pheno_files <- file.path(PHENO_DIR, sprintf("phenotype_%s.csv", GROUPS))
missing <- pheno_files[!file.exists(pheno_files)]
if (length(missing) == length(pheno_files))
  stop("no phenotype files found under ", PHENO_DIR, "; nothing to add the covariate to")
if (length(missing) > 0)
  cat(sprintf("  WARNING: %d phenotype file(s) not found: %s\n",
              length(missing), paste(basename(missing), collapse = ", ")))

for (i in seq_along(GROUPS)) {
  pf <- pheno_files[i]
  if (!file.exists(pf)) next

  pheno <- fread(pf)
  if ("total_dx_count" %in% names(pheno)) pheno[, total_dx_count := NULL]
  n0 <- nrow(pheno)
  pheno <- merge(pheno, dx_count, by = "eid", all.x = TRUE)
  n_no_fo <- sum(is.na(pheno$total_dx_count))
  pheno[is.na(total_dx_count), total_dx_count := 0L]   # no FO records -> 0 diagnoses
  stopifnot(nrow(pheno) == n0)                          # merge must not change N

  cat(sprintf("  %s: %s people, %s with no First-Occurrences record (set to 0)\n",
              GROUPS[i], format(n0, big.mark = ","), format(n_no_fo, big.mark = ",")))
  print(summary(pheno$total_dx_count))
  fwrite(pheno, pf)
  # Export and dx upload to RAP  (the phenotype file, now carrying total_dx_count)
}
cat("Done.\n")
