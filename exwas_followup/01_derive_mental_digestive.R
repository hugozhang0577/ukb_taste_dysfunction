#!/usr/bin/env Rscript
# Follow-up questionnaire exposures, block 1: load the follow-up extract, set non-answers to
# missing, and derive the digestive-health (Category 153) and mental-health (Category 136) exposures.
# 01-06 share one R session (02-06 use `raw`, `pheno_list` and save_section() defined here):
#   Rscript -e 'for (f in normalizePath(sort(Sys.glob("exwas_followup/0[1-6]_*.R")))) source(f)'
library(data.table)
library(stringr)

setwd(Sys.getenv("PROJECT_DIR", "."))
OUT_DIR <- "output/followup_exwas/derive/per_section"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

raw <- fread("input/raw/followup_exwas_merged.csv", na.strings = "NA")
for (cn in names(raw)[sapply(raw, is.character)])
  raw[get(cn) %in% c("", "Do not know", "Prefer not to answer", "do not wish to answer"),
      (cn) := NA_character_]

pheno_list <- sapply(c("group1", "group2", "group3"), function(g)
  fread(sprintf("input/analysis_ready/phenotype_%s.csv", g), na.strings = c("", "NA")),
  simplify = FALSE)

# Write one questionnaire block: its variable list (name, family, type) and, per cohort,
# the exposures of participants with a non-missing outcome.
save_section <- function(dt, family, var_type) {
  stem <- file.path(OUT_DIR, paste0("followup_exwas_",
                                    gsub("^_|_$", "", gsub("[^a-z0-9]+", "_", tolower(family)))))
  fwrite(data.table(var_name = names(var_type), domain = family, var_type = unname(var_type)),
         paste0(stem, "_variable_dict.csv"))
  for (g in names(pheno_list))
    fwrite(merge(dt[, c("eid", names(var_type)), with = FALSE],
                 pheno_list[[g]][!is.na(taste_2w_strict), .(eid)], by = "eid"),
           paste0(stem, "_", g, ".csv"))
}

recode_yesno <- function(x) fcase(x == "Yes", 1L, x == "No", 0L)
lead_number  <- function(x) as.numeric(str_extract(x, "^\\d+\\.?\\d*"))   # "0 - Not at all" -> 0
# sum of item scores; missing if any item is missing
item_sum <- function(cols, map) rowSums(vapply(cols, function(cn) unname(map[raw[[cn]]]),
                                               integer(nrow(raw))))

# Digestive health
dig <- raw[, .(eid)]
yes_no <- c(ibs_ever = "p21024", abdom_pain_6mo = "p21027", abdom_distension = "p21038",
            sensitive_stomach = "p21064", ibs_family_history = "p21065", coeliac_gluten = "p21068")
for (v in names(yes_no)) dig[, (v) := recode_yesno(raw[[yes_no[[v]]]])]
dig[, abdom_discomfort_freq := c("Never" = 0L, "Less than one day a month" = 1L,
                                 "One day a month" = 2L, "Two to three days a month" = 3L,
                                 "One day a week" = 4L, "More than one day a week" = 5L,
                                 "Every day" = 6L)[raw$p21025]]
numeric_items <- c(abdom_pain_severity = "p21036", bowel_satisfaction = "p21040",   # 0-10 ratings
                   gi_life_interference = "p21041", bowel_movements_day = "p21044")
for (v in names(numeric_items)) dig[, (v) := lead_number(raw[[numeric_items[[v]]]])]
phq12 <- c(phq12_headache = "p21051", phq12_dizziness = "p21053",       # PHQ-12 somatic items
           phq12_nausea = "p21059", phq12_fatigue = "p21060")
phq12_map <- c("Not bothered at all" = 0L, "Bothered a little" = 1L, "Bothered a lot" = 2L)
for (v in names(phq12)) dig[, (v) := phq12_map[raw[[phq12[[v]]]]]]
dig[, phq12_somatic_total := rowSums(.SD), .SDcols = names(phq12)]     # missing if any item is

save_section(dig, "Digestive health", c(
  ibs_ever = "binary", abdom_discomfort_freq = "ordinal", abdom_pain_6mo = "binary",
  abdom_pain_severity = "continuous", abdom_distension = "binary",
  bowel_satisfaction = "continuous", gi_life_interference = "continuous",
  bowel_movements_day = "continuous", sensitive_stomach = "binary",
  ibs_family_history = "binary", coeliac_gluten = "binary",
  phq12_headache = "ordinal", phq12_dizziness = "ordinal", phq12_nausea = "ordinal",
  phq12_fatigue = "ordinal", phq12_somatic_total = "derived_continuous"))

# Mental health: PHQ-9 and GAD-7 (items scored 0-3, caseness >= 10), AUDIT (caseness >= 8), mania
answer_504 <- c("Not at all" = 0L, "Several days" = 1L, "More than half the days" = 2L,
                "Nearly every day" = 3L)
audit_freq <- c("Never" = 0L, "Less than monthly" = 1L, "Monthly" = 2L, "Weekly" = 3L,
                "Daily or almost daily" = 4L)

mental <- raw[, .(eid)]
mental[, phq9_total := item_sum(c("p20514", "p20510", "p20517", "p20519", "p20511",
                                  "p20507", "p20508", "p20518", "p20513"), answer_504)]
mental[, gad7_total := item_sum(c("p20506", "p20509", "p20520", "p20515", "p20516",
                                  "p20505", "p20512"), answer_504)]
mental[, audit_total :=
  item_sum("p20414", c("Never" = 0L, "Monthly or less" = 1L, "2 to 4 times a month" = 2L,
                       "2 to 3 times a week" = 3L, "4 or more times a week" = 4L)) +
  item_sum("p20403", c("1 or 2" = 0L, "3 or 4" = 1L, "5 or 6" = 2L, "7, 8 or 9" = 3L,
                       "7, 8, or 9" = 3L, "10 or more" = 4L)) +
  item_sum(c("p20416", "p20413", "p20407", "p20412", "p20409", "p20408"), audit_freq) +
  item_sum(c("p20411", "p20405"), c("No" = 0L, "Yes, but not in the last year" = 2L,
                                    "Yes, during the last year" = 4L))]
mental[, `:=`(phq9_binary  = as.integer(phq9_total >= 10),
              gad7_binary  = as.integer(gad7_total >= 10),
              audit_binary = as.integer(audit_total >= 8),
              mania_excitable = recode_yesno(raw$p20501),
              mania_irritable = recode_yesno(raw$p20502),
              mania_duration  = c("Less than 24 hours" = 1L,
                                  "At least a day, but less than a week" = 2L,
                                  "A week or more" = 3L)[raw$p20492],
              mania_severity  = c("No problems" = 0L,
                                  "Needed treatment or caused problems with work, relationships, finances, the law or other aspects of life." = 1L)[raw$p20493],
              # number of manifestations ticked in the multi-select field
              mania_manifestation_count = fifelse(is.na(raw$p20548), NA_integer_,
                                                  as.integer(str_count(raw$p20548, "\\|") + 1L)))]

save_section(mental, "Mental health", c(
  phq9_total = "derived_continuous", phq9_binary = "derived_binary",
  gad7_total = "derived_continuous", gad7_binary = "derived_binary",
  audit_total = "derived_continuous", audit_binary = "derived_binary",
  mania_excitable = "binary", mania_irritable = "binary", mania_duration = "ordinal",
  mania_severity = "ordinal", mania_manifestation_count = "derived_continuous"))
