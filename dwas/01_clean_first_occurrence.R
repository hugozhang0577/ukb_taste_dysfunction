# Long table of prior diagnoses from the UK Biobank First Occurrence fields:
# one row per participant and three-character ICD-10 code, with first date and source.
library(data.table)
setwd(Sys.getenv("PROJECT_DIR", "."))

OUT_DIR <- "output/dwas/derive"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Each ICD-10 code has a "Date <code> first reported" field and a "Source of report" field
fo_meta <- fread("input/raw/FO_dic.csv")
fo_meta[, icd10 := sub(".*(\\b[A-Z][0-9]{2}\\b).*", "\\1", Description)]
fo_meta[, type  := fifelse(grepl("^Date", Description), "date", "source")]
fo_lookup <- merge(fo_meta[type == "date",   .(field_date   = field_id, icd10)],
                   fo_meta[type == "source", .(field_source = field_id, icd10)], by = "icd10")

# Merge the export batches by eid
merged <- NULL
for (f in list.files("input/raw/fo_batches", pattern = "\\.csv$", full.names = TRUE)) {
  dt <- fread(f)
  if (is.null(merged)) { merged <- dt; next }
  new_cols <- setdiff(names(dt), names(merged))
  if (length(new_cols) > 0)
    merged <- merge(merged, dt[, c("eid", new_cols), with = FALSE], by = "eid", all = TRUE)
}
fo_lookup <- fo_lookup[field_date %in% names(merged) & field_source %in% names(merged)]

to_long <- function(fields, value_name) {
  long <- melt(merged[, c("eid", fields), with = FALSE], id.vars = "eid",
               variable.name = "field", value.name = value_name, variable.factor = FALSE)
  long[, icd10 := setNames(fo_lookup$icd10, fields)[field]][, field := NULL]
  long[!is.na(get(value_name))]
}
fo_long <- merge(to_long(fo_lookup$field_date, "first_date"),
                 to_long(fo_lookup$field_source, "source_code"),
                 by = c("eid", "icd10"), all.x = TRUE)

# UK Biobank placeholder dates: 1909-09-09 / 2037-07-07 drop the record,
# 1900-01-01 to 1903-03-03 mean the date is unknown (record kept, date set to NA)
fo_long[, first_date := as.Date(first_date)]
fo_long <- fo_long[!first_date %in% as.Date(c("1909-09-09", "2037-07-07"))]
fo_long[first_date %in% as.Date(c("1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03")),
        first_date := as.IDate(NA)]

# Source of report (death register, primary care, hospital, self-report)
SOURCE_LABELS <- data.table(
  source_code      = c(20L, 21L, 30L, 31L, 40L, 41L, 50L, 51L),
  source_label     = c("death_only", "death_plus", "gp_only", "gp_plus",
                       "hes_only", "hes_plus", "selfreport_only", "selfreport_plus"),
  source_objective = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)
)
fo_long[, source_code := as.integer(source_code)]
fo_long <- merge(fo_long, SOURCE_LABELS, by = "source_code", all.x = TRUE)

saveRDS(fo_long, file.path(OUT_DIR, "fo_long.rds"))
