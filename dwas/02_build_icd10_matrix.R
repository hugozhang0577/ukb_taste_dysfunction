# Per-cohort binary matrix of prior diagnoses: three-character ICD-10 codes first
# recorded before the taste assessment date (or undated), plus a code dictionary.
library(data.table)
setwd(Sys.getenv("PROJECT_DIR", "."))

OUT_DIR <- "output/dwas/derive"
GROUPS  <- c("group1", "group2", "group3")

# ICD-10 chapters, and the reporting domains of the scan (by chapter letter)
CHAPTERS <- data.table(
  chapter_num = c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
                  "XI", "XII", "XIII", "XIV", "XV", "XVI", "XVII"),
  from = c("A00", "C00", "D50", "E00", "F00", "G00", "H00", "H60", "I00", "J00",
           "K00", "L00", "M00", "N00", "O00", "P00", "Q00"),
  to   = c("B99", "D48", "D89", "E90", "F99", "G99", "H59", "H95", "I99", "J99",
           "K93", "L99", "M99", "N99", "O99", "P96", "Q99"),
  chapter = c("Infectious and parasitic", "Neoplasms", "Blood and immune",
              "Endocrine, nutritional and metabolic", "Mental and behavioural",
              "Nervous system", "Eye and adnexa", "Ear and mastoid", "Circulatory",
              "Respiratory", "Digestive", "Skin and subcutaneous tissue",
              "Musculoskeletal and connective tissue", "Genitourinary",
              "Pregnancy and childbirth", "Perinatal", "Congenital")
)
code_chapter <- function(cd) vapply(cd, function(x) {
  k <- which(CHAPTERS$from <= x & x <= CHAPTERS$to)
  if (length(k) == 1) k else NA_integer_
}, integer(1))
domain_of <- function(cd) {
  l <- substr(cd, 1, 1)
  fcase(l == "F", "C-Mental", l == "G", "C-Neuro", l == "E", "C-Endocrine",
        l == "I", "C-Circulatory", l == "J", "C-Respiratory", l == "K", "C-Digestive",
        default = "C-Other")
}

fo_long <- readRDS(file.path(OUT_DIR, "fo_long.rds"))
time_dt <- fread("input/analysis_ready/tastetime_assesstime.csv", select = c("eid", "p28755"))
setnames(time_dt, "p28755", "taste_date")
time_dt[, taste_date := as.IDate(taste_date)]

lab <- fread("input/reference/coding19.tsv", select = c("coding", "meaning"))[nchar(coding) == 3]
lab[, description := sub("^[A-Z][0-9]{2} ", "", meaning)]

group_eids <- lapply(setNames(GROUPS, GROUPS), function(g)
  unique(fread(sprintf("input/analysis_ready/phenotype_%s.csv", g), select = "eid")$eid))

# Keep diagnoses dated before the taste assessment, or undated
fo <- fo_long[eid %in% unlist(group_eids), .(eid, icd10, first_date = as.IDate(first_date))]
fo <- merge(fo, time_dt, by = "eid", all.x = TRUE)
fo_prior <- unique(fo[is.na(first_date) | first_date < taste_date, .(eid, icd10)])

# Exposure matrix per cohort: one logical column per code, no frequency floor
all_codes <- character(0)
for (g in GROUPS) {
  ids   <- group_eids[[g]]
  f     <- fo_prior[eid %in% ids]
  codes <- sort(unique(f$icd10))
  all_codes <- union(all_codes, codes)
  rows <- split(match(f$eid, ids), f$icd10)
  mat  <- data.table(eid = ids)
  for (cd in codes) {
    v <- logical(length(ids)); v[rows[[cd]]] <- TRUE
    set(mat, j = cd, value = v)
  }
  saveRDS(mat, file.path(OUT_DIR, sprintf("icd10_matrix_%s.rds", g)))
}

# Code dictionary: data-coding 19 name, ICD-10 chapter, scan domain, oral (K00-K14)
dict <- merge(data.table(icd10 = sort(all_codes)), lab[, .(icd10 = coding, description)],
              by = "icd10", all.x = TRUE)
k <- code_chapter(dict$icd10)
dict[, `:=`(chapter_num = CHAPTERS$chapter_num[k], chapter = CHAPTERS$chapter[k],
            domain = domain_of(icd10), oral = grepl("^K(0[0-9]|1[0-4])$", icd10))]
fwrite(dict, file.path(OUT_DIR, "icd10_dictionary.csv"))
# The matrices and the dictionary are copied to input/analysis_ready/ for 04_run_dwas.R
