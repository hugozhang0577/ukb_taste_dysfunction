#> in   input/raw/ukb_base_fields.csv
#> in   input/reference/taste_affecting_drugs_literature.csv
#> in   input/reference/taste_drugs_full_list.csv
#> ship output/base_table/base_table_full.csv -> input/analysis_ready/
#> ship output/base_table/base_table_covar.csv -> input/analysis_ready/
suppressPackageStartupMessages(library(data.table))

project_dir <- Sys.getenv("PROJECT_DIR", unset = ".")
CODE_DIR    <- Sys.getenv("CODE_DIR",    unset = ".")
EXTRACT     <- Sys.getenv("UKB_EXTRACT", unset = "ukb_base_fields.csv")
INPUT_DIR   <- file.path(project_dir, "input", "raw")
REF_DIR     <- file.path(project_dir, "input", "reference")
OUTPUT_DIR  <- file.path(project_dir, "output", "base_table")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Step 0. Load the extract and check it against the requested field list
requested <- trimws(readLines(file.path(CODE_DIR, "field_list.txt")))
requested <- requested[nzchar(requested)]

data <- as.data.frame(fread(file.path(INPUT_DIR, EXTRACT), showProgress = FALSE))

names(data) <- sub("^participant\\.", "", names(data))

absent <- setdiff(requested, names(data))
# Age at cancer diagnosis carries as many instances as the release holds.
absent <- grep("^p40008_i", absent, value = TRUE, invert = TRUE)
if (length(absent)) {
  stop("Fields absent from the extract: ", paste(absent, collapse = ", "))
}

data <- data[, intersect(requested, names(data)), drop = FALSE]
cat(sprintf("extract: %d rows x %d cols\n", nrow(data), ncol(data)))

recode_labelled <- function(x, labels) {
  xc  <- trimws(as.character(x))
  out <- unname(labels[xc])
  num <- suppressWarnings(as.numeric(xc))
  out[is.na(out)] <- num[is.na(out)]
  out
}

# Yes/no fields: 1 yes, 0 no, -1 do not know, -3 prefer not to answer -> NA.
yes_no <- function(x) {
  xc <- trimws(as.character(x))
  ifelse(xc %in% c("1", "Yes"), 1, ifelse(xc %in% c("0", "No"), 0, NA_real_))
}

drop_codes <- function(x, codes) {
  x[!is.na(x) & x %in% codes] <- NA
  x
}

# Attach derived columns keyed on eid, then release the source columns.
attach_derived <- function(base, derived, source_cols) {
  stopifnot(identical(as.character(base$eid), as.character(derived$eid)))
  new_cols <- setdiff(names(derived), "eid")
  base <- cbind(base, as.data.frame(derived)[, new_cols, drop = FALSE])
  base[, setdiff(names(base), source_cols), drop = FALSE]
}

# Step 1.1 Taste / smell severity outcome fields
n_before <- nrow(data)
data <- data[data$p28615 %in% c(0, 1), , drop = FALSE]
cat(sprintf("definite yes/no on taste_change: %d -> %d rows (%d dropped)\n",
            n_before, nrow(data), n_before - nrow(data)))

rename_map <- c(
  p28615 = "taste_change", p28616 = "taste_time",  p28617 = "taste_extent",
  p28612 = "smell_change", p28613 = "smell_time",  p28614 = "smell_extent"
)
for (old in names(rename_map)) names(data)[names(data) == old] <- rename_map[[old]]

# Step 1.2 Demographics
names(data)[names(data) == "p31"] <- "sex"
data$age <- 2023L - data$p34
names(data)[names(data) == "p21003_i0"] <- "age_baseline"
names(data)[names(data) == "p22000"] <- "batch_number"
data$batch_number <- as.factor(data$batch_number)
cat(sprintf("genotyping batches: %d levels\n", nlevels(data$batch_number)))
data <- data[, setdiff(names(data), "p34"), drop = FALSE]

# Step 1.3 Anthropometry, blood pressure, spirometry, BMI
names(data)[names(data) == "p48_i0"] <- "Waist_circumference"
names(data)[names(data) == "p49_i0"] <- "Hip_circumference"
names(data)[names(data) == "p50_i0"] <- "Standing_height"
names(data)[names(data) == "p74_i0"] <- "Fasting_time"

# Pulse and blood pressure: mean of the two automated readings.
data$paulse <- rowMeans(cbind(data$p102_i0_a0,  data$p102_i0_a1),  na.rm = TRUE)
data$DBP    <- rowMeans(cbind(data$p4079_i0_a0, data$p4079_i0_a1), na.rm = TRUE)
data$SBP    <- rowMeans(cbind(data$p4080_i0_a0, data$p4080_i0_a1), na.rm = TRUE)
data <- data[, setdiff(names(data), c(
  "p102_i0_a0",  "p102_i0_a1",
  "p4079_i0_a0", "p4079_i0_a1",
  "p4080_i0_a0", "p4080_i0_a1")), drop = FALSE]

names(data)[names(data) == "p20150_i0"] <- "FEV1"
names(data)[names(data) == "p20151_i0"] <- "FVC"
names(data)[names(data) == "p20258"]    <- "FEV1_FVC_ratio"
names(data)[names(data) == "p21001_i0"] <- "BMI"

# Step 1.4 Blood counts, urine assay, blood biochemistry (baseline instance).
biomarker_map <- c(
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
  p30870_i0 = "Triglycerides", p30880_i0 = "Urate", p30890_i0 = "VD"
)
for (fid in names(biomarker_map)) {
  names(data)[names(data) == fid] <- biomarker_map[[fid]]
}

# Step 1.5 Smoking status (binary: 1 non-current, 0 current; missing kept separate).
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
# The two variants differ only in how the unclassifiable participants are handled.
data$smoking      <- ifelse(data$smoking_unknown == 1, 0, smoking_bin)  # primary: missing -> current
data$smoking_sens <- ifelse(data$smoking_unknown == 1, 1, smoking_bin)  # sensitivity: missing -> non-current

# Step 1.6 Alcohol (weekly units -> sex-specific high/low-risk drinking).
alcohol_vars <- c("p1568_i0", "p1578_i0", "p1588_i0", "p1598_i0", "p1608_i0", "p5364_i0")
data[alcohol_vars] <- lapply(data[alcohol_vars], function(x) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  ifelse(is.na(x) | x %in% c(-1, -3), 0, x)
})
data$alcohol_unit <- (data$p1568_i0 * 2 + data$p1578_i0 * 2 + data$p1588_i0 * 2 +
                      data$p1598_i0 * 1 + data$p1608_i0 * 2 + data$p5364_i0 * 2) / 7
# 1 = high-risk drinking, 0 = low-risk (sex-specific weekly-unit threshold).
data$drink <- ifelse(
  (data$sex == 1 & data$alcohol_unit <= 2) | (data$sex == 0 & data$alcohol_unit <= 1), 0, 1)

# Step 1.7 Hypertension (measured BP OR self-reported diagnosis/medication).
labels_6150 <- c("Heart attack" = 1, "Angina" = 2, "Stroke" = 3,
                 "High blood pressure" = 4, "None of the above" = -7,
                 "Prefer not to answer" = -3)
labels_med  <- c("Cholesterol lowering medication" = 1, "Blood pressure medication" = 2,
                 "Insulin" = 3, "Hormone replacement therapy" = 4,
                 "Oral contraceptive pill or minipill" = 5, "None of the above" = -7,
                 "Do not know" = -1, "Prefer not to answer" = -3)

for (src in c("p6150_i0", "p6150_i1", "p6150_i2")) {
  data[[sub("^p", "n", src)]] <- recode_labelled(data[[src]], labels_6150)
}
for (src in c("p6153_i0", "p6177_i0")) {
  data[[sub("^p", "n", src)]] <- recode_labelled(data[[src]], labels_med)
}

data$hypertension <- with(data, ifelse(
  SBP >= 140 | DBP >= 90 | n6150_i0 == 4 | n6150_i1 == 4 | n6150_i2 == 4, 1, NA))
data$hypertension <- with(data, ifelse(
  (SBP < 140 & DBP > 0) & (DBP < 90 & DBP > 0), 0, hypertension))
data$hypertension <- with(data, ifelse(
  (sex == 0 & n6153_i0 == 2) | (sex == 1 & n6177_i0 == 2), 1, hypertension))

# Step 1.8 Diabetes (self-reported doctor diagnosis, any instance of p2443).
diab_vars <- paste0("p2443_i", 0:3)
data[diab_vars] <- lapply(data[diab_vars], yes_no)
data$diabetes <- ifelse(rowSums(data[diab_vars], na.rm = TRUE) > 0, 1, 0)

# Step 1.9 Prevalent cancer (age at any cancer diagnosis <= baseline age; p40008).
cancer_age_cols <- grep("^p40008_i", names(data), value = TRUE)
data$cancer <- as.integer(Reduce(`|`, lapply(cancer_age_cols, function(cc)
  !is.na(data[[cc]]) & data[[cc]] <= data$age_baseline)))

cat(sprintf("derived: hypertension %d / diabetes %d / prevalent cancer %d cases\n",
            sum(data$hypertension == 1, na.rm = TRUE),
            sum(data$diabetes == 1,     na.rm = TRUE),
            sum(data$cancer == 1,       na.rm = TRUE)))

# Step 2. Oral-health covariates (field 6149, mouth/teeth dental problems)
detect_multiresponse_code <- function(x, code) {
  if (is.na(x) || x == "" || x == "-7") return(FALSE)
  code %in% strsplit(as.character(x), "\\|")[[1]]
}

detect_multiresponse_code_vec <- function(x, code) {
  as.integer(vapply(x, detect_multiresponse_code, logical(1), code = as.character(code)))
}

derive_oral_health_covariates <- function(oral) {
  instance_cols <- paste0("p6149_i", 0:3)

  labels <- list(
    bleeding_gums = 1,
    painful_gums = 2,
    loose_teeth = 3,
    toothache = 4,
    mouth_ulcers = 5,
    dentures = 6
  )

  out <- oral[, .(eid)]
  for (nm in names(labels)) {
    code <- labels[[nm]]
    out[, paste0("oral_", nm, "_baseline") := detect_multiresponse_code_vec(oral$p6149_i0, code)]
    inst <- lapply(instance_cols, function(v) detect_multiresponse_code_vec(oral[[v]], code))
    out[, paste0("oral_", nm) := as.integer(Reduce(`+`, inst) > 0)]
  }

  base_vars <- paste0("oral_", names(labels), "_baseline")
  full_vars <- paste0("oral_", names(labels))

  out[, periodontal_indicator_baseline := as.integer(
    oral_bleeding_gums_baseline == 1 |
      oral_painful_gums_baseline == 1 |
      oral_loose_teeth_baseline == 1
  )]
  out[, periodontal_indicator := as.integer(
    oral_bleeding_gums == 1 |
      oral_painful_gums == 1 |
      oral_loose_teeth == 1
  )]

  out[, oral_problem_count_baseline := rowSums(.SD), .SDcols = base_vars]
  out[, oral_problem_count := rowSums(.SD), .SDcols = full_vars]
  out[, any_oral_problem_baseline := as.integer(oral_problem_count_baseline > 0)]
  out[, any_oral_problem := as.integer(oral_problem_count > 0)]

  out[]
}

oral_cols <- paste0("p6149_i", 0:3)
data <- attach_derived(
  data,
  derive_oral_health_covariates(as.data.table(data[, c("eid", oral_cols)])),
  oral_cols)
cat(sprintf("oral health: %d with any problem at baseline\n",
            sum(data$any_oral_problem_baseline == 1, na.rm = TRUE)))

# Step 3. Oral and dental surgery covariates (OPCS-4 operative procedures)
taste_affecting_codes <- list(
  lip = c("F01", "F02", "F03", "F04", "F05", "F06"),
  tongue = c("F22", "F23", "F24", "F26"),
  palate = c("F28", "F29", "F30", "F32"),
  jaw = c("F18"),
  tonsil = c("F34", "F36"),
  salivary = c("F44", "F45", "F46", "F48", "F50", "F51", "F52", "F53", "F55", "F56", "F58")
)

non_taste_affecting_codes <- list(
  tooth_extraction = c("F08", "F09", "F10"),
  preprosthetic = c("F11"),
  apex = c("F12"),
  restoration = c("F13"),
  orthodontic = c("F14", "F15"),
  tooth_other = c("F16", "F17"),
  gingiva = c("F20"),
  mouth_other = c("F38", "F39", "F40", "F42", "F43"),
  denture = c("F63")
)

has_any_opcs <- function(x, codes) {
  if (is.na(x) || x == "") return(FALSE)
  pattern <- paste0("(^|\\|)(", paste(codes, collapse = "|"), ")")
  grepl(pattern, x)
}

has_any_opcs_vec <- function(x, codes) {
  as.integer(vapply(x, has_any_opcs, logical(1), codes = codes))
}

filter_codes_before_baseline <- function(codes, dates, baseline_date) {
  if (is.na(codes) || codes == "" || is.na(baseline_date)) return("")
  code_vec <- strsplit(as.character(codes), "\\|")[[1]]
  n <- min(length(code_vec), length(dates))
  if (n == 0) return("")
  keep <- !is.na(dates[seq_len(n)]) & dates[seq_len(n)] < baseline_date
  paste(code_vec[seq_len(n)][keep], collapse = "|")
}

derive_oral_surgery_covariates <- function(surgery) {
  date_cols <- grep("^p41282_a", names(surgery), value = TRUE)
  if (length(date_cols) == 0) stop("No procedure-date columns matching ^p41282_a were found")

  surgery[, baseline_date := as.IDate(p53_i0)]
  date_matrix <- as.matrix(surgery[, ..date_cols])
  date_matrix <- apply(date_matrix, 2, as.IDate)

  surgery[, p41272_baseline := {
    out <- character(.N)
    for (i in seq_len(.N)) {
      out[i] <- filter_codes_before_baseline(p41272[i], date_matrix[i, ], baseline_date[i])
    }
    out
  }]

  out <- surgery[, .(eid)]
  for (nm in names(taste_affecting_codes)) {
    out[, paste0("surg_", nm, "_full") := has_any_opcs_vec(surgery$p41272, taste_affecting_codes[[nm]])]
    out[, paste0("surg_", nm, "_baseline") := has_any_opcs_vec(surgery$p41272_baseline, taste_affecting_codes[[nm]])]
  }
  for (nm in names(non_taste_affecting_codes)) {
    out[, paste0("surg_", nm, "_full") := has_any_opcs_vec(surgery$p41272, non_taste_affecting_codes[[nm]])]
    out[, paste0("surg_", nm, "_baseline") := has_any_opcs_vec(surgery$p41272_baseline, non_taste_affecting_codes[[nm]])]
  }

  taste_full <- paste0("surg_", names(taste_affecting_codes), "_full")
  taste_base <- paste0("surg_", names(taste_affecting_codes), "_baseline")
  non_full <- paste0("surg_", names(non_taste_affecting_codes), "_full")
  non_base <- paste0("surg_", names(non_taste_affecting_codes), "_baseline")

  out[, surg_taste_affecting_full := as.integer(rowSums(.SD) > 0), .SDcols = taste_full]
  out[, surg_taste_affecting_baseline := as.integer(rowSums(.SD) > 0), .SDcols = taste_base]
  out[, surg_non_taste_affecting_full := as.integer(rowSums(.SD) > 0), .SDcols = non_full]
  out[, surg_non_taste_affecting_baseline := as.integer(rowSums(.SD) > 0), .SDcols = non_base]

  out[]
}

surgery_cols <- c("p41272", "p53_i0", grep("^p41282_a", names(data), value = TRUE))
data <- attach_derived(
  data,
  derive_oral_surgery_covariates(as.data.table(data[, c("eid", surgery_cols)])),
  surgery_cols)
cat(sprintf("oral surgery: %d taste-affecting before baseline / %d lifetime\n",
            sum(data$surg_taste_affecting_baseline == 1, na.rm = TRUE),
            sum(data$surg_taste_affecting_full == 1,     na.rm = TRUE)))

# Step 4. Taste-affecting medication covariates (field 20003, self-reported)
drug_codes <- fread(file.path(REF_DIR, "taste_drugs_full_list.csv"))
drug_review <- fread(file.path(REF_DIR, "taste_affecting_drugs_literature.csv"))
high_risk_names <- drug_review[high_risk == TRUE, drug_name]
med_codes <- unique(as.character(drug_codes[drug_name %in% high_risk_names, ukb_code]))
cat(sprintf("medication codes: %d high-risk of %d published\n",
            length(med_codes), uniqueN(drug_codes$ukb_code)))

derive_medication_covariates <- function(medication, codes) {
  med_cols_all <- grep("^p20003_", names(medication), value = TRUE)
  if (length(med_cols_all) == 0) stop("No medication columns matching ^p20003_")
  med_cols_baseline <- grep("^p20003_i0_", med_cols_all, value = TRUE)
  if (length(med_cols_baseline) == 0) stop("No baseline medication columns matching ^p20003_i0_")

  any_code <- function(sd) {
    as.integer(apply(sd, 1, function(x) any(as.character(x) %in% codes, na.rm = TRUE)))
  }

  medication[, taste_drug_user_numeric := any_code(.SD), .SDcols = med_cols_all]
  medication[, taste_drug_user_baseline_numeric := any_code(.SD), .SDcols = med_cols_baseline]

  medication[, .(eid, taste_drug_user_numeric, taste_drug_user_baseline_numeric)]
}

med_cols <- grep("^p20003_", names(data), value = TRUE)
data <- attach_derived(
  data,
  derive_medication_covariates(as.data.table(data[, c("eid", med_cols)]), med_codes),
  med_cols)
cat(sprintf("medication: %d users at baseline / %d across instances\n",
            sum(data$taste_drug_user_baseline_numeric == 1, na.rm = TRUE),
            sum(data$taste_drug_user_numeric == 1,          na.rm = TRUE)))

# Step 5. Outputs
fwrite(data, file.path(OUTPUT_DIR, "base_table_full.csv"))

drop_cols <- c(
  cancer_age_cols,
  "p1239_i0", "p1249_i0", "p2644_i0",
  alcohol_vars,
  paste0("p2443_i", 0:3),
  paste0("p6150_i", 0:2), "p6153_i0", "p6177_i0",
  grep("^n6150|^n6153|^n6177", names(data), value = TRUE),
  "alcohol_unit", "smoking_unknown", "Rheumatoid_factor"
)
data_covar <- data[, setdiff(names(data), drop_cols), drop = FALSE]
fwrite(data_covar, file.path(OUTPUT_DIR, "base_table_covar.csv"))

cat(sprintf("base_table_full.csv : %d rows x %d cols\n", nrow(data), ncol(data)))
cat(sprintf("base_table_covar.csv: %d rows x %d cols\n", nrow(data_covar), ncol(data_covar)))
