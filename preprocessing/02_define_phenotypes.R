# Taste outcome definitions (primary taste_2w_strict and four variants) and smell features
# from the taste and smell questionnaire fields.
library(data.table)
setwd(Sys.getenv("PROJECT_DIR", "."))
dir.create("output/gwas_phenotypes", recursive = TRUE, showWarnings = FALSE)
dir.create("output/smell_features", recursive = TRUE, showWarnings = FALSE)

negative_to_na <- function(x) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  x[x < 0] <- NA_real_
  x
}

base_tbl <- fread("input/analysis_ready/base_table_full.csv",
                  select = c("eid", "taste_change", "taste_time", "taste_extent",
                             "smell_change", "smell_time", "smell_extent"))

taste <- copy(base_tbl)
taste[, `:=`(
  taste_change       = suppressWarnings(as.numeric(as.character(taste_change))),
  smell_change       = suppressWarnings(as.numeric(as.character(smell_change))),
  taste_time_clean   = negative_to_na(taste_time),
  taste_extent_clean = negative_to_na(taste_extent)
)]

# For the taste outcome, a do-not-know / prefer-not-to-answer smell response counts as no smell change
taste[smell_change < 0, smell_change := 0]
taste[, only_smell := as.integer(taste_change == 0 & smell_change == 1)]

# Controls (strict): neither taste nor smell change; smell-only participants are excluded (NA)
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
taste[, has_time_extent_info := fcase(
  taste_change == 0, TRUE,
  taste_change == 1, !is.na(taste_time_clean) | !is.na(taste_extent_clean),
  default = FALSE
)]

# Primary case: taste change lasting >= 2 weeks (time code >= 1) OR affecting daily life (extent code 1)
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

# Variant: lasting >= 4 weeks (time code >= 2) OR affecting daily life
taste[, meets_4w_extent := taste_change == 1 &
        ((!is.na(taste_time_clean) & taste_time_clean >= 2) |
           (!is.na(taste_extent_clean) & taste_extent_clean == 1))]
taste[, taste_4w_strict := fcase(
  taste_change == 0 & smell_change == 0, 0,
  meets_4w_extent == TRUE, 1,
  taste_change == 1 & has_time_extent_info == TRUE, NA_real_,
  default = NA_real_
)]

fwrite(taste[, .(eid, taste_change, smell_change, taste_time_clean, taste_extent_clean, only_smell,
                 has_time_extent_info, meets_2w_extent, meets_4w_extent, taste_basic_strict,
                 taste_basic_relaxed, taste_2w_strict, taste_2w_relaxed, taste_4w_strict)],
       "output/gwas_phenotypes/taste_gwas_phenotypes.csv")

# Smell features: any change, and its duration / extent among those reporting a change
smell <- base_tbl[, .(eid, smell_change, smell_time, smell_extent)]
for (v in c("smell_change", "smell_time", "smell_extent")) {
  smell[, (v) := suppressWarnings(as.integer(as.character(get(v))))]
  smell[get(v) %in% c(-1L, -3L), (v) := NA_integer_]
}
smell[, smell_any := fcase(smell_change == 1L, 1L, smell_change == 0L, 0L, default = NA_integer_)]
smell[is.na(smell_any) | smell_any == 0L, `:=`(smell_time = NA_integer_, smell_extent = NA_integer_)]
fwrite(smell[, .(eid, smell_any, smell_time, smell_extent)], "output/smell_features/smell_features.csv")
