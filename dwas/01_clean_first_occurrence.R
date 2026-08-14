#!/usr/bin/env Rscript
# =============================================================================
# Disease-wide association scan — clean the First Occurrences export into a long table
# =============================================================================
#
# FIRST step of the PheCode derivation:
#   merge export batches -> melt to long -> clean dates -> decode source
#   -> fo_long (eid | icd10 | first_date | source_code | source_label | source_objective)
#
# The First Occurrences resource records, per participant, the earliest date of
# each ICD-10 three-character code, integrated across hospital inpatient,
# primary-care, cancer- and death-registry, and baseline self-report sources.
# Each code has a "Date <code> first reported" field and a matching
# "Source of report of <code>" field; the two are paired by ICD-10 code here.
#
# Inputs (under $PROJECT_DIR):
#   input/raw/fo_batches/*.csv   Table-Exporter batches (eid + p<id> date/source fields)
#   input/raw/FO_dic.csv         First-Occurrences field dictionary (field_id, Description)
#   field_list.txt               the export specification, checked against the
#                                dictionary pairing below
# Output (under $PROJECT_DIR/output/dwas/derive/):
#   fo_long.rds (+ .csv), fo_lookup.rds, fo_icd10_frequency.csv
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# =============================================================================
# Configuration ($PROJECT_DIR-based; override as needed)
# =============================================================================
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
CODE_DIR    <- Sys.getenv("CODE_DIR", unset = ".")
BATCH_DIR   <- Sys.getenv("FO_BATCH_DIR", unset = file.path(PROJECT_DIR, "input", "raw", "fo_batches"))
FO_META_FILE <- Sys.getenv("FO_DIC", unset = file.path(PROJECT_DIR, "input", "raw", "FO_dic.csv"))
OUTPUT_DIR  <- Sys.getenv("DERIVE_DIR",
                          unset = file.path(PROJECT_DIR, "output", "dwas", "derive"))
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Step 0: build the field -> ICD-10 lookup from the FO dictionary
#   Each ICD-10 code has a "Date <code> first reported" field and a
#   "Source of report of <code>" field; pair them by ICD-10 code.
# =============================================================================
cat("=== Step 0: build FO lookup ===\n")
fo_meta <- fread(FO_META_FILE)

id_col   <- grep("field", names(fo_meta), ignore.case = TRUE, value = TRUE)[1]
desc_col <- grep("desc",  names(fo_meta), ignore.case = TRUE, value = TRUE)[1]
if (is.na(id_col))   id_col   <- names(fo_meta)[1]
if (is.na(desc_col)) desc_col <- names(fo_meta)[2]

fo_meta[, icd10 := sub(".*(\\b[A-Z][0-9]{2}\\b).*", "\\1", get(desc_col))]
fo_meta[, type  := fifelse(grepl("^Date", get(desc_col)), "date", "source")]

fo_date   <- fo_meta[type == "date",   .(field_date = get(id_col), icd10)]
fo_source <- fo_meta[type == "source", .(field_source = get(id_col), icd10)]
fo_lookup <- merge(fo_date, fo_source, by = "icd10")
cat(sprintf("  lookup: %d ICD-10 codes\n", nrow(fo_lookup)))

# -----------------------------------------------------------------------------
# field_list.txt is the export specification: the field names this scan needs,
# one per line, in the form the RAP Table Exporter takes (an `eid` line followed
# by field names). It is checked against the pairing just built from the
# dictionary, so the list a reproducer exports and the fields this step expects
# cannot drift apart. Both directions stop the run and name the fields.
# -----------------------------------------------------------------------------
FIELD_LIST <- file.path(CODE_DIR, "field_list.txt")
requested_fields <- setdiff(trimws(readLines(FIELD_LIST)), c("", "eid"))
paired_fields    <- sort(unique(c(fo_lookup$field_date, fo_lookup$field_source)))
only_list <- setdiff(requested_fields, paired_fields)
only_dict <- setdiff(paired_fields, requested_fields)
if (length(only_list) > 0 || length(only_dict) > 0) {
  stop("field_list.txt and the First-Occurrences dictionary disagree.",
       if (length(only_list)) paste0("\n  in field_list.txt only: ",
                                     paste(head(sort(only_list), 20), collapse = ", "),
                                     if (length(only_list) > 20)
                                       sprintf(" ... (%d total)", length(only_list))),
       if (length(only_dict)) paste0("\n  in the dictionary only: ",
                                     paste(head(sort(only_dict), 20), collapse = ", "),
                                     if (length(only_dict) > 20)
                                       sprintf(" ... (%d total)", length(only_dict))))
}
cat(sprintf("  field_list.txt agrees with the dictionary: %d fields (%d ICD-10 codes)\n",
            length(requested_fields), nrow(fo_lookup)))

# =============================================================================
# Step 1: merge all export batches by eid
# =============================================================================
cat("\n=== Step 1: merge batches ===\n")
csv_files <- list.files(BATCH_DIR, pattern = "\\.csv$", full.names = TRUE)
cat(sprintf("  %d batch files\n", length(csv_files)))

merged <- NULL
for (f in csv_files) {
  dt <- fread(f)
  if (is.null(merged)) {
    merged <- dt
  } else {
    new_cols <- setdiff(names(dt), names(merged))
    if (length(new_cols) > 0)
      merged <- merge(merged, dt[, c("eid", new_cols), with = FALSE], by = "eid", all = TRUE)
  }
}
cat(sprintf("  merged: %d people x %d cols\n", nrow(merged), ncol(merged)))

# Keep only lookup entries whose fields are actually present. field_list.txt
# already agrees with the dictionary, so anything dropped here means the export
# is short of its own specification — report it rather than letting codes
# disappear silently.
n_codes_before <- nrow(fo_lookup)
fo_lookup <- fo_lookup[field_date %in% names(merged) & field_source %in% names(merged)]
cat(sprintf("  usable mappings: %d -> %d ICD-10 codes\n", n_codes_before, nrow(fo_lookup)))
if (nrow(fo_lookup) < n_codes_before) {
  dropped <- setdiff(fo_date$icd10, fo_lookup$icd10)
  cat(sprintf("  WARNING: %d code(s) in field_list.txt are absent from the export: %s\n",
              length(dropped), paste(head(sort(dropped), 20), collapse = ", ")))
}

# =============================================================================
# Step 2-4: melt date and source columns to long, then merge
# =============================================================================
cat("\n=== Step 2-4: melt to long ===\n")
date_to_icd <- setNames(fo_lookup$icd10, fo_lookup$field_date)
fo_dates <- melt(merged[, c("eid", fo_lookup$field_date), with = FALSE],
                 id.vars = "eid", variable.name = "field",
                 value.name = "first_date", variable.factor = FALSE)
fo_dates[, icd10 := date_to_icd[field]][, field := NULL]
fo_dates <- fo_dates[!is.na(first_date)]

source_to_icd <- setNames(fo_lookup$icd10, fo_lookup$field_source)
fo_sources <- melt(merged[, c("eid", fo_lookup$field_source), with = FALSE],
                   id.vars = "eid", variable.name = "field",
                   value.name = "source_code", variable.factor = FALSE)
fo_sources[, icd10 := source_to_icd[field]][, field := NULL]
fo_sources <- fo_sources[!is.na(source_code)]

fo_long <- merge(fo_dates, fo_sources, by = c("eid", "icd10"), all.x = TRUE)
rm(fo_dates, fo_sources); gc()
cat(sprintf("  fo_long: %s records\n", format(nrow(fo_long), big.mark = ",")))

# =============================================================================
# Step 5: clean dates — handle UK Biobank sentinel codes
#   1909-09-09 / 2037-07-07  -> placeholder, drop the record
#   1900/1901/1902/1903      -> date unreliable, keep record but set date NA
# =============================================================================
cat("\n=== Step 5: clean dates ===\n")
fo_long[, first_date := as.Date(first_date)]
n_before <- nrow(fo_long)
fo_long[, date_flag := fcase(
  first_date %in% as.Date(c("1909-09-09", "2037-07-07")), "remove",
  first_date %in% as.Date(c("1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03")), "date_unknown",
  default = "valid"
)]
fo_long <- fo_long[date_flag != "remove"]
fo_long[date_flag == "date_unknown", first_date := as.IDate(NA)]
fo_long[, date_flag := NULL]
cat(sprintf("  %s -> %s (removed %s placeholder records)\n",
            format(n_before, big.mark = ","), format(nrow(fo_long), big.mark = ","),
            format(n_before - nrow(fo_long), big.mark = ",")))

# =============================================================================
# Step 6: decode source codes (objective = hospital/GP/death; subjective = self-report)
# =============================================================================
cat("\n=== Step 6: decode source ===\n")
SOURCE_LABELS <- data.table(
  source_code  = c(20L, 21L, 30L, 31L, 40L, 41L, 50L, 51L),
  source_label = c("death_only", "death_plus", "gp_only", "gp_plus",
                   "hes_only", "hes_plus", "selfreport_only", "selfreport_plus"),
  source_objective = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)
)
fo_long[, source_code := as.integer(source_code)]
fo_long <- merge(fo_long, SOURCE_LABELS, by = "source_code", all.x = TRUE)

# =============================================================================
# Step 7: QC summary
# =============================================================================
cat("\n=== Step 7: QC ===\n")
cat(sprintf("  records: %s | participants: %s | ICD-10 codes: %d\n",
            format(nrow(fo_long), big.mark = ","),
            format(uniqueN(fo_long$eid), big.mark = ","),
            uniqueN(fo_long$icd10)))
cat(sprintf("  records with a valid date: %.1f%%\n", mean(!is.na(fo_long$first_date)) * 100))

# =============================================================================
# Step 8: save
# =============================================================================
cat("\n=== Step 8: save ===\n")
saveRDS(fo_long, file.path(OUTPUT_DIR, "fo_long.rds"))
fwrite(fo_long,  file.path(OUTPUT_DIR, "fo_long.csv"))
saveRDS(fo_lookup, file.path(OUTPUT_DIR, "fo_lookup.rds"))

icd_freq <- fo_long[, .(
  n_persons    = uniqueN(eid),
  n_records    = .N,
  pct_has_date = round(mean(!is.na(first_date)) * 100, 1)
), by = icd10][order(-n_persons)]
fwrite(icd_freq, file.path(OUTPUT_DIR, "fo_icd10_frequency.csv"))
# Export and dx upload to RAP  (fo_long.rds/.csv, fo_lookup.rds, the frequency table)
cat(sprintf("  saved fo_long.rds + companions to %s\n", OUTPUT_DIR))
