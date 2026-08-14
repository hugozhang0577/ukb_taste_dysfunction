#!/usr/bin/env Rscript
# =============================================================================
# 05_covid_timing_flags.R  (SARS-CoV-2 test history, timed against the outcome)
# =============================================================================
#
# Per-participant SARS-CoV-2 positivity flags, each defined relative to the date
# the participant completed the health and well-being questionnaire. Figure 5E
# and the subtype analysis in ../subtype_characterisation/05_covid_association.R
# read the output.
#
# This step is separate from 01_build_base_table.R because the test results are
# not participant fields. They are record-level tables, one row per test, held as
# their own entities and exported separately (see covid_field_list.txt). A
# participant appears once per test, or not at all.
#
# The whole point is the timing. The outcome is a questionnaire answered on a
# known date, so an infection recorded AFTER that date cannot have caused the
# reported symptom; only `before` flags are usable as exposures, and the `post`
# flag exists to make that asymmetry visible rather than to be modelled.
#
# What these flags are not: they record a documented positive test, not
# infection. Testing was neither universal nor random, and its intensity changed
# over the pandemic, so `n_tests_before_taste` is carried alongside as a crude
# handle on how much a participant was tested at all.
#
# Usage: Rscript 05_covid_timing_flags.R
# Input:  input/raw/covid19_result_{england,scotland,wales}.csv
#           per covid_field_list.txt: eid, specdate, result
#         input/analysis_ready/tastetime_assesstime.csv   (eid, p53_i0, p28755)
# Output: output/covid/covid_temporal_flags.csv
# =============================================================================

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
suppressPackageStartupMessages(library(data.table))

REGIONS   <- c("england", "scotland", "wales")
RAW_DIR   <- "input/raw"
TIME_FILE <- "input/analysis_ready/tastetime_assesstime.csv"
OUT_DIR   <- "output/covid"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

DISTANT_DAYS <- 30L   # "distant" means at least this long before the questionnaire

# ---- [1] the three regional test tables -------------------------------------
read_region <- function(region) {
  f <- file.path(RAW_DIR, sprintf("covid19_result_%s.csv", region))
  if (!file.exists(f)) stop("test-result export not found: ", f,
                            "\n  export the covid19_result_", region,
                            " entity with covid_field_list.txt")
  d <- fread(f)
  miss <- setdiff(c("eid", "specdate", "result"), names(d))
  if (length(miss))
    stop(f, " is missing: ", paste(miss, collapse = ", "),
         "\n  the export must carry all three columns in one file; two separate ",
         "exports joined by row order would attach dates to the wrong tests")
  d <- d[, .(eid = as.integer(eid), specdate = as.IDate(specdate),
             result = as.integer(result), region = region)]
  cat(sprintf("  %-9s %8d tests, %7d positive, %s to %s\n", region, nrow(d),
              sum(d$result == 1, na.rm = TRUE),
              min(d$specdate, na.rm = TRUE), max(d$specdate, na.rm = TRUE)))
  d
}
cat("[1] test records\n")
tests <- rbindlist(lapply(REGIONS, read_region))

n0 <- nrow(tests)
tests <- tests[!is.na(eid) & !is.na(specdate) & !is.na(result)]
cat(sprintf("  pooled %d tests, %d dropped for a missing id, date or result\n",
            n0, n0 - nrow(tests)))
bad <- setdiff(unique(tests$result), c(0L, 1L))
if (length(bad)) stop("result is not 0/1; found: ", paste(bad, collapse = ", "))
cat(sprintf("  %d participants tested at least once\n", uniqueN(tests$eid)))

# ---- [2] the questionnaire date ---------------------------------------------
if (!file.exists(TIME_FILE)) stop("date anchors not found: ", TIME_FILE)
anchors <- fread(TIME_FILE)
setnames(anchors, "p28755", "taste_date", skip_absent = TRUE)
if (!"taste_date" %in% names(anchors))
  stop(TIME_FILE, " has no p28755 column (questionnaire completion date)")
anchors <- anchors[, .(eid = as.integer(eid), taste_date = as.IDate(taste_date))]
anchors <- anchors[!is.na(taste_date)]
cat(sprintf("[2] %d participants with a questionnaire date (%s to %s)\n",
            nrow(anchors), min(anchors$taste_date), max(anchors$taste_date)))

# ---- [3] the flags -----------------------------------------------------------
# Every participant with a questionnaire date gets a row, including those who
# were never tested: absence of a test is informative here (it is a zero, not a
# missing value), and dropping them would silently restrict the denominator to
# the tested.
d <- merge(anchors, tests, by = "eid", all.x = TRUE, allow.cartesian = TRUE)
cat(sprintf("[3] %d participant-test rows; %d participants never tested\n",
            nrow(d), anchors[!eid %in% tests$eid, .N]))

pos <- d[result == 1]
flags <- d[, .(
  phe_ever_tested        = as.integer(any(!is.na(specdate))),
  n_tests_before_taste   = sum(!is.na(specdate) & specdate < taste_date)
), by = .(eid, taste_date)]

pos_flags <- pos[, .(
  phe_ever_positive          = 1L,
  phe_pos_before_taste       = as.integer(any(specdate <  taste_date)),
  phe_pos_distant_pre_taste  = as.integer(any(specdate <  taste_date - DISTANT_DAYS)),
  phe_pos_post_taste         = as.integer(any(specdate >  taste_date)),
  n_positive_before_taste    = sum(specdate < taste_date),
  first_positive_date        = min(specdate)
), by = eid]

out <- merge(flags, pos_flags, by = "eid", all.x = TRUE)
for (v in c("phe_ever_positive", "phe_pos_before_taste", "phe_pos_distant_pre_taste",
            "phe_pos_post_taste", "n_positive_before_taste"))
  set(out, which(is.na(out[[v]])), v, 0L)
out[, first_positive_days_to_taste :=
      as.integer(taste_date - first_positive_date)]   # positive = before the questionnaire
out[, first_positive_date := NULL]

cat("\n=== flags ===\n")
print(out[, .(n = .N,
              tested        = sum(phe_ever_tested),
              ever_positive = sum(phe_ever_positive),
              pos_before    = sum(phe_pos_before_taste),
              pos_after     = sum(phe_pos_post_taste))])
cat(sprintf("positivity among the tested: %.1f%%\n",
            100 * out[phe_ever_tested == 1, mean(phe_ever_positive)]))

# A participant cannot be positive before the questionnaire without being
# positive at all; catching that here beats finding it in an odds ratio.
stopifnot(out[phe_pos_before_taste == 1 & phe_ever_positive == 0, .N] == 0,
          out[n_positive_before_taste > n_tests_before_taste, .N] == 0)

fwrite(out, file.path(OUT_DIR, "covid_temporal_flags.csv"))
# Export and dx upload to RAP  (the timing flags the subtype analysis reads)
cat("\nwrote", file.path(OUT_DIR, "covid_temporal_flags.csv"), "\n")
cat("copy it to input/analysis_ready/ for ../subtype_characterisation/05\n")
