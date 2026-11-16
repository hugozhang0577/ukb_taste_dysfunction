# Base table: taste/smell questionnaire fields, demographics, physical measurements and blood
# assays, plus derived covariates (smoking, alcohol, hypertension, diabetes, cancer, oral health,
# oral surgery, taste-affecting medication), for participants with a yes/no answer on taste change.
library(data.table)
requested <- trimws(readLines("preprocessing/field_list.txt"))
requested <- requested[nzchar(requested)]
setwd(Sys.getenv("PROJECT_DIR", "."))
dir.create("output/base_table", recursive = TRUE, showWarnings = FALSE)

data <- as.data.frame(fread("input/raw/ukb_base_fields.csv"))
data <- data[, intersect(requested, names(data)), drop = FALSE]

recode_labelled <- function(x, labels) {
  xc  <- trimws(as.character(x))
  out <- unname(labels[xc])
  num <- suppressWarnings(as.numeric(xc))
  out[is.na(out)] <- num[is.na(out)]
  out
}
# 1 yes, 0 no; do not know / prefer not to answer -> NA
yes_no <- function(x) {
  xc <- trimws(as.character(x))
  ifelse(xc %in% c("1", "Yes"), 1, ifelse(xc %in% c("0", "No"), 0, NA_real_))
}
drop_codes <- function(x, codes) { x[!is.na(x) & x %in% codes] <- NA; x }
# Append derived columns (same row order) and drop their source columns
attach_derived <- function(base, derived, source_cols) {
  base <- cbind(base, as.data.frame(derived)[, setdiff(names(derived), "eid"), drop = FALSE])
  base[, setdiff(names(base), source_cols), drop = FALSE]
}

# Definite yes/no answer on loss or change in taste (field 28615)
data <- data[data$p28615 %in% c(0, 1), , drop = FALSE]

rename_map <- c(
  p28615 = "taste_change", p28616 = "taste_time", p28617 = "taste_extent",
  p28612 = "smell_change", p28613 = "smell_time", p28614 = "smell_extent",
  p31 = "sex", p21003_i0 = "age_baseline", p22000 = "batch_number",
  p48_i0 = "Waist_circumference", p49_i0 = "Hip_circumference", p50_i0 = "Standing_height",
  p74_i0 = "Fasting_time", p20150_i0 = "FEV1", p20151_i0 = "FVC", p20258 = "FEV1_FVC_ratio",
  p21001_i0 = "BMI",
  # blood counts
  p30000_i0 = "WBC_count", p30010_i0 = "RBC_count", p30020_i0 = "Hb", p30030_i0 = "Hct",
  p30040_i0 = "MCV", p30050_i0 = "MCH", p30060_i0 = "MCHC", p30070_i0 = "RDW",
  p30080_i0 = "Palatelet_count", p30100_i0 = "MPV", p30120_i0 = "Lymphocyte_count",
  p30130_i0 = "Monocyte_count", p30140_i0 = "Neutrophill_count", p30150_i0 = "Eosinophill_count",
  p30160_i0 = "Basophill_count", p30170_i0 = "Nucleated_RBC_count", p30250_i0 = "Reticulocyte_count",
  p30260_i0 = "MRV", p30270_i0 = "MSCV", p30280_i0 = "IRF", p30300_i0 = "HLR",
  # urine
  p30500_i0 = "Microalbumin_urine", p30510_i0 = "Creatinine_enzy_urine",
  p30520_i0 = "Potassium_urine", p30530_i0 = "Sodium_urine",
  # biochemistry
  p30600_i0 = "Albumin", p30610_i0 = "ALP", p30620_i0 = "ALT", p30630_i0 = "Apol_A",
  p30640_i0 = "Apol_B", p30650_i0 = "AST", p30660_i0 = "bilirubin", p30670_i0 = "Urea",
  p30680_i0 = "Calcium", p30690_i0 = "Cholesterol", p30700_i0 = "Creatinine",
  p30710_i0 = "C_reactive_protein", p30720_i0 = "Cystatin_C", p30730_i0 = "GGT",
  p30740_i0 = "Glucose", p30750_i0 = "HbA1c", p30760_i0 = "HDL_cholesterol",
  p30770_i0 = "IGF_1", p30780_i0 = "LDL_direct", p30790_i0 = "Lipo_A", p30800_i0 = "Oestradiol",
  p30810_i0 = "Phosphate", p30820_i0 = "Rheumatoid_factor", p30830_i0 = "SHBG",
  p30840_i0 = "Total_bilirubin", p30850_i0 = "Testosterone", p30860_i0 = "Total_protein",
  p30870_i0 = "Triglycerides", p30880_i0 = "Urate", p30890_i0 = "VD")
names(data)[match(names(rename_map), names(data))] <- unname(rename_map)

# Age in 2023 (taste questionnaire) from year of birth (field 34)
data$age <- 2023L - data$p34
data$batch_number <- as.factor(data$batch_number)
data <- data[, setdiff(names(data), "p34"), drop = FALSE]

# Pulse and blood pressure: mean of the two automated readings
data$paulse <- rowMeans(cbind(data$p102_i0_a0,  data$p102_i0_a1),  na.rm = TRUE)
data$DBP    <- rowMeans(cbind(data$p4079_i0_a0, data$p4079_i0_a1), na.rm = TRUE)
data$SBP    <- rowMeans(cbind(data$p4080_i0_a0, data$p4080_i0_a1), na.rm = TRUE)
data <- data[, setdiff(names(data), c("p102_i0_a0", "p102_i0_a1", "p4079_i0_a0", "p4079_i0_a1",
                                      "p4080_i0_a0", "p4080_i0_a1")), drop = FALSE]

# Smoking (fields 1239 current, 1249 past, 2644 occasional): 1 non-current, 0 current
cur  <- drop_codes(data$p1239_i0, -3)
past <- drop_codes(data$p1249_i0, -3)
occ  <- drop_codes(data$p2644_i0, c(-3, -1))
smoking_bin <- rep(NA_real_, nrow(data))
smoking_bin[!is.na(cur) & cur == 0 & (
  (!is.na(past) & past == 4) |
  (past %in% c(2, 3)    & !is.na(occ) & occ == 0) |
  (past %in% c(1, 2, 3) & !is.na(occ) & occ == 1))] <- 1
smoking_bin[cur %in% c(1, 2)] <- 0
data$smoking_bin     <- smoking_bin
data$smoking_unknown <- as.integer(is.na(smoking_bin))
data$smoking      <- ifelse(data$smoking_unknown == 1, 0, smoking_bin)  # primary: unclassifiable -> current
data$smoking_sens <- ifelse(data$smoking_unknown == 1, 1, smoking_bin)  # sensitivity: unclassifiable -> non-current

# Alcohol: average daily drinks from weekly intake by drink type; 1 = above 2 (men) or 1 (women)
alcohol_vars <- c("p1568_i0", "p1578_i0", "p1588_i0", "p1598_i0", "p1608_i0", "p5364_i0")
data[alcohol_vars] <- lapply(data[alcohol_vars], function(x) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  ifelse(is.na(x) | x %in% c(-1, -3), 0, x)
})
data$alcohol_unit <- (data$p1568_i0 * 2 + data$p1578_i0 * 2 + data$p1588_i0 * 2 +
                      data$p1598_i0 * 1 + data$p1608_i0 * 2 + data$p5364_i0 * 2) / 7
data$drink <- ifelse(
  (data$sex == 1 & data$alcohol_unit <= 2) | (data$sex == 0 & data$alcohol_unit <= 1), 0, 1)

# Hypertension: measured BP >= 140/90, or self-reported diagnosis (6150) or medication (6153/6177)
labels_6150 <- c("Heart attack" = 1, "Angina" = 2, "Stroke" = 3, "High blood pressure" = 4,
                 "None of the above" = -7, "Prefer not to answer" = -3)
labels_med  <- c("Cholesterol lowering medication" = 1, "Blood pressure medication" = 2,
                 "Insulin" = 3, "Hormone replacement therapy" = 4,
                 "Oral contraceptive pill or minipill" = 5, "None of the above" = -7,
                 "Do not know" = -1, "Prefer not to answer" = -3)
for (src in c("p6150_i0", "p6150_i1", "p6150_i2"))
  data[[sub("^p", "n", src)]] <- recode_labelled(data[[src]], labels_6150)
for (src in c("p6153_i0", "p6177_i0"))
  data[[sub("^p", "n", src)]] <- recode_labelled(data[[src]], labels_med)
data$hypertension <- with(data, ifelse(
  SBP >= 140 | DBP >= 90 | n6150_i0 == 4 | n6150_i1 == 4 | n6150_i2 == 4, 1, NA))
data$hypertension <- with(data, ifelse(
  (SBP < 140 & DBP > 0) & (DBP < 90 & DBP > 0), 0, hypertension))
data$hypertension <- with(data, ifelse(
  (sex == 0 & n6153_i0 == 2) | (sex == 1 & n6177_i0 == 2), 1, hypertension))

# Diabetes: self-reported doctor diagnosis at any instance (field 2443)
diab_vars <- paste0("p2443_i", 0:3)
data[diab_vars] <- lapply(data[diab_vars], yes_no)
data$diabetes <- ifelse(rowSums(data[diab_vars], na.rm = TRUE) > 0, 1, 0)

# Prevalent cancer: any cancer diagnosed at or before the baseline age (field 40008)
cancer_age_cols <- grep("^p40008_i", names(data), value = TRUE)
data$cancer <- as.integer(Reduce(`|`, lapply(cancer_age_cols, function(cc)
  !is.na(data[[cc]]) & data[[cc]] <= data$age_baseline)))

# Oral health (field 6149, mouth/teeth problems): each problem at baseline and at any instance
has_code <- function(x, code) as.integer(vapply(x, function(s)
  !(is.na(s) || s == "" || s == "-7") && code %in% strsplit(as.character(s), "\\|")[[1]],
  logical(1)))
oral_items <- c(bleeding_gums = 1, painful_gums = 2, loose_teeth = 3, toothache = 4,
                mouth_ulcers = 5, dentures = 6)
oral_cols <- paste0("p6149_i", 0:3)
oral <- data.table(eid = data$eid)
for (nm in names(oral_items)) {
  code <- as.character(oral_items[[nm]])
  oral[, paste0("oral_", nm, "_baseline") := has_code(data$p6149_i0, code)]
  oral[, paste0("oral_", nm) := as.integer(
    Reduce(`+`, lapply(oral_cols, function(v) has_code(data[[v]], code))) > 0)]
}
oral[, periodontal_indicator_baseline := as.integer(
  oral_bleeding_gums_baseline == 1 | oral_painful_gums_baseline == 1 | oral_loose_teeth_baseline == 1)]
oral[, periodontal_indicator := as.integer(
  oral_bleeding_gums == 1 | oral_painful_gums == 1 | oral_loose_teeth == 1)]
oral[, oral_problem_count_baseline := rowSums(.SD), .SDcols = paste0("oral_", names(oral_items), "_baseline")]
oral[, oral_problem_count := rowSums(.SD), .SDcols = paste0("oral_", names(oral_items))]
oral[, any_oral_problem_baseline := as.integer(oral_problem_count_baseline > 0)]
oral[, any_oral_problem := as.integer(oral_problem_count > 0)]
data <- attach_derived(data, oral, oral_cols)

# Oral and dental surgery (OPCS-4 procedures, field 41272): lifetime and before the baseline visit
taste_affecting <- list(
  lip = c("F01", "F02", "F03", "F04", "F05", "F06"), tongue = c("F22", "F23", "F24", "F26"),
  palate = c("F28", "F29", "F30", "F32"), jaw = "F18", tonsil = c("F34", "F36"),
  salivary = c("F44", "F45", "F46", "F48", "F50", "F51", "F52", "F53", "F55", "F56", "F58"))
non_taste_affecting <- list(
  tooth_extraction = c("F08", "F09", "F10"), preprosthetic = "F11", apex = "F12",
  restoration = "F13", orthodontic = c("F14", "F15"), tooth_other = c("F16", "F17"),
  gingiva = "F20", mouth_other = c("F38", "F39", "F40", "F42", "F43"), denture = "F63")

has_opcs <- function(x, codes) {
  pattern <- paste0("(^|\\|)(", paste(codes, collapse = "|"), ")")
  as.integer(vapply(x, function(s) !(is.na(s) || s == "") && grepl(pattern, s), logical(1)))
}
# Keep the procedures whose date (field 41282, same order) is before the baseline date
codes_before <- function(codes, dates, baseline_date) {
  if (is.na(codes) || codes == "" || is.na(baseline_date)) return("")
  code_vec <- strsplit(as.character(codes), "\\|")[[1]]
  n <- min(length(code_vec), length(dates))
  keep <- !is.na(dates[seq_len(n)]) & dates[seq_len(n)] < baseline_date
  paste(code_vec[seq_len(n)][keep], collapse = "|")
}

surgery_cols <- c("p41272", "p53_i0", grep("^p41282_a", names(data), value = TRUE))
surgery <- as.data.table(data[, c("eid", surgery_cols)])
date_cols <- grep("^p41282_a", names(surgery), value = TRUE)
surgery[, baseline_date := as.IDate(p53_i0)]
date_matrix <- apply(as.matrix(surgery[, ..date_cols]), 2, as.IDate)
codes_baseline <- vapply(seq_len(nrow(surgery)), function(i)
  codes_before(surgery$p41272[i], date_matrix[i, ], surgery$baseline_date[i]), character(1))

surg <- data.table(eid = surgery$eid)
for (grp in list(taste_affecting, non_taste_affecting)) for (nm in names(grp)) {
  surg[, paste0("surg_", nm, "_full")     := has_opcs(surgery$p41272, grp[[nm]])]
  surg[, paste0("surg_", nm, "_baseline") := has_opcs(codes_baseline, grp[[nm]])]
}
any_of <- function(sd) as.integer(rowSums(sd) > 0)
surg[, surg_taste_affecting_full := any_of(.SD),
     .SDcols = paste0("surg_", names(taste_affecting), "_full")]
surg[, surg_taste_affecting_baseline := any_of(.SD),
     .SDcols = paste0("surg_", names(taste_affecting), "_baseline")]
surg[, surg_non_taste_affecting_full := any_of(.SD),
     .SDcols = paste0("surg_", names(non_taste_affecting), "_full")]
surg[, surg_non_taste_affecting_baseline := any_of(.SD),
     .SDcols = paste0("surg_", names(non_taste_affecting), "_baseline")]
data <- attach_derived(data, surg, surgery_cols)

# Taste-affecting medication (field 20003): high-risk drugs from the literature review
drug_codes  <- fread("input/reference/taste_drugs_full_list.csv")
drug_review <- fread("input/reference/taste_affecting_drugs_literature.csv")
med_codes <- unique(as.character(
  drug_codes[drug_name %in% drug_review[high_risk == TRUE, drug_name], ukb_code]))
med_cols <- grep("^p20003_", names(data), value = TRUE)
med <- as.data.table(data[, c("eid", med_cols)])
any_code <- function(sd)
  as.integer(apply(sd, 1, function(x) any(as.character(x) %in% med_codes, na.rm = TRUE)))
med[, taste_drug_user_numeric := any_code(.SD), .SDcols = med_cols]
med[, taste_drug_user_baseline_numeric := any_code(.SD),
    .SDcols = grep("^p20003_i0_", med_cols, value = TRUE)]
data <- attach_derived(data, med[, .(eid, taste_drug_user_numeric, taste_drug_user_baseline_numeric)],
                       med_cols)

fwrite(data, "output/base_table/base_table_full.csv")
drop_cols <- c(cancer_age_cols, "p1239_i0", "p1249_i0", "p2644_i0", alcohol_vars, diab_vars,
               paste0("p6150_i", 0:2), "p6153_i0", "p6177_i0",
               grep("^n6150|^n6153|^n6177", names(data), value = TRUE),
               "alcohol_unit", "smoking_unknown", "Rheumatoid_factor")
fwrite(data[, setdiff(names(data), drop_cols), drop = FALSE], "output/base_table/base_table_covar.csv")
