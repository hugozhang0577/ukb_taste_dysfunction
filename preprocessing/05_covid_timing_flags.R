# SARS-CoV-2 test timing relative to the taste questionnaire date: ever tested, ever positive,
# and positive tests before, more than 30 days before, and after the questionnaire.
library(data.table)
setwd(Sys.getenv("PROJECT_DIR", "."))
dir.create("output/covid", recursive = TRUE, showWarnings = FALSE)

DISTANT_DAYS <- 30L

tests <- rbindlist(lapply(c("england", "scotland", "wales"), function(region) {
  d <- fread(sprintf("input/raw/covid19_result_%s.csv", region))
  d[, .(eid = as.integer(eid), specdate = as.IDate(specdate), result = as.integer(result), region = region)]
}))
tests <- tests[!is.na(eid) & !is.na(specdate) & !is.na(result)]

# Questionnaire completion date (field 28755)
anchors <- fread("input/analysis_ready/tastetime_assesstime.csv")
setnames(anchors, "p28755", "taste_date")
anchors <- anchors[, .(eid = as.integer(eid), taste_date = as.IDate(taste_date))]
anchors <- anchors[!is.na(taste_date)]

d <- merge(anchors, tests, by = "eid", all.x = TRUE, allow.cartesian = TRUE)
flags <- d[, .(
  phe_ever_tested      = as.integer(any(!is.na(specdate))),
  n_tests_before_taste = sum(!is.na(specdate) & specdate < taste_date)
), by = .(eid, taste_date)]
pos_flags <- d[result == 1, .(
  phe_ever_positive         = 1L,
  phe_pos_before_taste      = as.integer(any(specdate <  taste_date)),
  phe_pos_distant_pre_taste = as.integer(any(specdate <  taste_date - DISTANT_DAYS)),
  phe_pos_post_taste        = as.integer(any(specdate >  taste_date)),
  n_positive_before_taste   = sum(specdate < taste_date),
  first_positive_date       = min(specdate)
), by = eid]

out <- merge(flags, pos_flags, by = "eid", all.x = TRUE)
for (v in c("phe_ever_positive", "phe_pos_before_taste", "phe_pos_distant_pre_taste",
            "phe_pos_post_taste", "n_positive_before_taste"))
  set(out, which(is.na(out[[v]])), v, 0L)
out[, first_positive_days_to_taste := as.integer(taste_date - first_positive_date)]   # > 0: before the questionnaire
out[, first_positive_date := NULL]
fwrite(out, "output/covid/covid_temporal_flags.csv")
