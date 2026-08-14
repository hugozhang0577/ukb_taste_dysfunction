#!/usr/bin/env Rscript
# ============================================================================
# Baseline ExWAS — exposure derivation: cleaning, recoding, dictionary, and
#                  per-cohort analysis-ready matrices
# ============================================================================
#
# Turns the baseline questionnaire and assessment-centre export into the
# analysis-ready exposure matrices the regression engine consumes, together
# with the variable dictionary that tells the engine each variable's type and
# its false-discovery-rate family.
#
# Sections 1-13 clean and recode the raw fields and write one file per exposure
# family per cohort. Sections 14-19 read those files back, merge them into one
# matrix per cohort, attach case counts, and write the summary report. The
# per-family files are kept rather than passed in memory: they are inspectable
# outputs in their own right, and section 15 audits them before the merge, so a
# family that failed to write is caught rather than silently dropped.
#
# No variable is excluded here for missingness or rarity. The case-count
# threshold (n_cases >= 100) is applied downstream in the regression engine, so
# that the dictionary stays a complete record of what was derived.
#
# Exposure families (the false-discovery-rate groups; Benjamini-Hochberg is
# applied within each, never pooled across them, because the families differ
# enormously in size and pooling would let one large family set the threshold
# for all the others):
#   Demographics & SES                 Lifestyle & sleep
#   Dietary intake & preferences       Anthropometric & physiological
#   Oral health                        General & mental health
#   Sensory function & pain            Clinical screening & treatment
#   Blood biochemistry & haematology   Environmental exposures
#   Handedness & laterality            Reproductive (female-only)
# The family name is the value stored in the dictionary's `source` column, so
# results tables and figure legends read directly from it with no lookup.
# Reproductive items are female-only and form their own family, refitted on
# female participants without a sex term.
#
# Inputs
#   input/raw/baseline_exwas_table_export.csv   Table Exporter output
#   input/analysis_ready/phenotype_group{1,2,3}.csv
#     already restricted to the taste cohort (taste_2w_strict in {0,1});
#     the phenotype file stays the single source of truth for outcome and
#     covariates across all four scans, so the exposure matrices below carry
#     exposures and eid only
#
# Outputs
#   output/baseline_exwas/derive/per_source/    one file per family per cohort,
#                                              plus per-family missingness and
#                                              dictionaries
#   output/baseline_exwas/derive/
#     exwas_baseline_{group}.csv                 merged exposure matrix
#     baseline_exwas_variable_dictionary.csv     unified dictionary
#     baseline_exwas_combined_missingness.csv
#     baseline_exwas_variable_case_counts.csv
#     baseline_exwas_source_overlap.csv
#     baseline_exwas_summary_report.txt
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

# ============================================================================
# 0. CONFIGURATION — paths resolved under $PROJECT_DIR (override as needed)
# ============================================================================
# Expected layout:
#   $PROJECT_DIR/input/raw/        Table Exporter output + Field 1920 export
#   $PROJECT_DIR/input/analysis_ready/  per-cohort phenotype files
#   $PROJECT_DIR/output/baseline_exwas/derive/per_source/   (written here)
# The per-cohort phenotype files are already restricted to the taste cohort
# (taste_2w_strict in {0,1}).

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
RAW_DIR     <- Sys.getenv("RAW_DIR",     unset = file.path(PROJECT_DIR, "input", "raw"))
PHENO_DIR   <- Sys.getenv("PHENO_DIR",   unset = file.path(PROJECT_DIR, "input", "analysis_ready"))
PER_SOURCE_DIR  <- Sys.getenv("DERIVE_PER_SOURCE_DIR",
                          unset = file.path(PROJECT_DIR, "output",
                                            "baseline_exwas", "derive", "per_source"))

EXPORT_CSV <- file.path(RAW_DIR, "baseline_exwas_table_export.csv")

PHENO_FILES <- list(
  group1 = file.path(PHENO_DIR, "phenotype_group1.csv"),
  group2 = file.path(PHENO_DIR, "phenotype_group2.csv"),
  group3 = file.path(PHENO_DIR, "phenotype_group3.csv")
)


# Exposure families: the false-discovery-rate groups. Benjamini-Hochberg is
# applied within each family and never pooled across them, because the families
# differ enormously in size and pooling would let one large family set the
# threshold for all the others. The family name is the value stored in the
# dictionary, so results tables and figure legends need no lookup.
# Reproductive items are female-only and form their own family; they are
# excluded from the list below because that sweep is run separately.
FAMILIES <- c(
  "Demographics & SES", "Lifestyle & sleep", "Dietary intake & preferences",
  "Anthropometric & physiological", "Oral health", "General & mental health",
  "Sensory function & pain", "Clinical screening & treatment",
  "Blood biochemistry & haematology", "Environmental exposures",
  "Handedness & laterality"
)

# Family names are display strings, so the per-family file names use a slug of
# the name rather than a second vocabulary that could drift away from it.
family_slug <- function(x) gsub("_+", "_", gsub("[^A-Za-z0-9]+", "_", tolower(x)))

DERIVE_DIR <- Sys.getenv("EXPOSURE_DIR",
                         unset = file.path(PROJECT_DIR, "output",
                                           "baseline_exwas", "derive"))

dir.create(PER_SOURCE_DIR, showWarnings = FALSE, recursive = TRUE)

cat("================================================================\n")
cat("  ExWAS baseline ExWAS — Data Cleaning Pipeline\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

# ============================================================================
# 1. DATA LOADING
# ============================================================================

cat("--- Step 1: Loading data ---\n")

raw <- fread(EXPORT_CSV, na.strings = c("", "NA"))
cat("  baseline ExWAS raw:", nrow(raw), "rows x", ncol(raw), "columns\n")

# Load all group phenotype files
pheno_list <- list()
for (gname in names(PHENO_FILES)) {
  fp <- PHENO_FILES[[gname]]
  if (file.exists(fp)) {
    pheno_list[[gname]] <- fread(fp, na.strings = c("", "NA"))
    cat("  Phenotype", gname, ":", nrow(pheno_list[[gname]]), "rows x",
        ncol(pheno_list[[gname]]), "cols\n")
  } else {
    cat("  WARNING:", gname, "file not found:", fp, "- skipping\n")
  }
}

pheno <- pheno_list[["group1"]]

# ============================================================================
# 2. FIELD SELECTION — Apply REVIEW decisions + exclude non-EXPORT
# ============================================================================

cat("\n--- Step 2: Field selection ---\n")

# --- Define EXPORT field IDs by domain ---
# Based on audit report + REVIEW decisions

field_registry <- list(
  
  demographics = list(
    fields = c(845, 6138, 21000, 738, 6142),
    labels = c("age_education", "qualifications", "ethnic_background",
               "household_income", "employment_status")
  ),
  
  lifestyle = list(
    fields = c(1160, 1170, 1180, 1190, 1200, 1210, 1220, 20161, 22040),
    labels = c("sleep_duration", "getting_up", "chronotype", "nap_daytime",
               "insomnia", "snoring", "daytime_dozing", "pack_years",
               "physical_activity_met")
  ),
  
  dietary = list(
    fields = c(1289, 1299, 1309, 1319, 1329, 1339, 1349, 1359, 1369, 1379,
               1389, 1408, 1418, 1428, 1438, 1448, 1458, 1468, 1478, 1488,
               1498, 1508, 1518, 1528, 1538, 1548, 2654, 3680, 6144),
    labels = c("cooked_veg", "salad_veg", "fresh_fruit", "dried_fruit",
               "oily_fish", "nonoily_fish", "processed_meat", "poultry",
               "beef", "lamb", "pork", "cheese", "milk_type", "spread_type",
               "bread_intake", "bread_type", "cereal_intake", "cereal_type",
               "salt_added", "tea_intake", "coffee_intake", "coffee_type",
               "hot_drink_temp", "water_intake", "diet_change_5yr",
               "diet_variation", "spread_details", "age_last_meat",
               "never_eat")
  ),
  
  anthropometric = list(
    fields = c(46, 47, 48, 49, 50, 95, 1687, 4079, 4080, 21002,
               23099, 23100, 23101, 23127),
    labels = c("grip_left", "grip_right", "waist_circ", "hip_circ",
               "height", "pulse_rate", "body_size_age10",
               "dbp", "sbp", "weight",
               "body_fat_pct", "whole_body_fat", "fat_free_mass",
               "trunk_fat_pct")
  ),
  
  reproductive_female = list(
    fields = c(2714, 2724, 3581),
    labels = c("age_menarche", "had_menopause", "age_menopause"),
    female_only = TRUE
  ),
  
  general_health = list(
    fields = c(137, 2178, 2188, 2296, 2306),
    labels = c("n_medications", "overall_health", "longstanding_illness",
               "falls_last_year", "weight_change_1yr")
  ),
  
  mental_health = list(
    fields = c(1031, 1920, 2020, 2050, 2060, 2110, 4598, 4631, 6145, 20127),
    labels = c("friend_visits", "mood_swings", "loneliness",
               "depressed_mood_2wk", "unenthusiasm_2wk", "able_confide",
               "ever_depressed_2wk", "ever_unenthusiastic_2wk",
               "life_stressors", "neuroticism_score")
  ),
  
  sensory_pain = list(
    fields = c(2247, 2956, 3404, 3414, 3571, 3741, 3773, 3799,
               4067, 4803, 6159, 20019, 20021),
    labels = c("hearing_difficulty", "general_pain_3mo", "neck_pain_3mo",
               "hip_pain_3mo", "back_pain_3mo", "stomach_pain_3mo",
               "knee_pain_3mo", "headache_3mo", "facial_pain_3mo",
               "tinnitus", "pain_types_last_month",
               "srt_left", "srt_right")
  ),
  
  clinical_screen = list(
    fields = c(2443, 2453, 4501, 6150, 20016, 20018, 20023),
    labels = c("diabetes_sr", "cancer_sr", "family_nonaccidental_death",
               "vascular_problems", "fluid_intelligence",
               "prospective_memory", "reaction_time")
  ),
  
  biochemistry = list(
    fields = c(30600, 30610, 30620, 30630, 30640, 30650, 30660, 30670,
               30680, 30690, 30700, 30710, 30720, 30730, 30740, 30750,
               30760, 30770, 30780, 30790, 30800, 30810, 30820, 30830,
               30840, 30850, 30860, 30870, 30880, 30890),
    labels = c("albumin", "alp", "alt", "apoa", "apob", "ast",
               "direct_bilirubin", "urea", "calcium", "cholesterol",
               "creatinine", "crp", "cystatin_c", "ggt", "glucose",
               "hba1c", "hdl", "igf1", "ldl", "lpa",
               "oestradiol", "phosphate", "rheumatoid_factor", "shbg",
               "total_bilirubin", "testosterone", "total_protein",
               "triglycerides", "urate", "vitamin_d")
  ),
  
  haematology = list(
    fields = c(30000, 30010, 30020, 30030, 30040, 30050, 30060, 30070,
               30080, 30090, 30100, 30110, 30120, 30130, 30140, 30150,
               30160, 30170, 30180, 30190, 30200, 30210, 30220, 30230,
               30240, 30250, 30260, 30270, 30280, 30290, 30300),
    labels = c("wbc", "rbc", "haemoglobin", "haematocrit", "mcv", "mch",
               "mchc", "rdw", "platelet", "platelet_crit", "mpv",
               "platelet_dw", "lymphocyte", "monocyte", "neutrophil",
               "eosinophil", "basophil", "nrbc", "lymphocyte_pct",
               "monocyte_pct", "neutrophil_pct", "eosinophil_pct",
               "basophil_pct", "nrbc_pct", "reticulocyte_pct",
               "reticulocyte", "mean_retic_vol", "mean_sphered_vol",
               "immature_retic", "hls_retic_pct", "hls_retic")
  ),
  
  environment = list(
    fields = c(24006),
    labels = c("pm25_2010")
  ),
  
  handedness = list(
    fields = c(1707),
    labels = c("handedness")
  )
)

# Build flat lookup: field_id → (domain, label, female_only)
field_lookup <- rbindlist(lapply(names(field_registry), function(dom) {
  reg <- field_registry[[dom]]
  data.table(
    field_id    = reg$fields,
    domain      = dom,
    var_label   = reg$labels,
    female_only = isTRUE(reg$female_only)
  )
}))

cat("  Registered fields:", nrow(field_lookup), "\n")

# field_list.txt is the export specification: the RAP column names this scan
# needs, one per line, spelled out in full because the Table Exporter matches
# names literally. The field ids behind those names are checked against the
# registry above, so the list a reproducer exports and the fields this script
# expects cannot drift apart.
FIELD_LIST <- file.path(Sys.getenv("CODE_DIR", unset = "."), "field_list.txt")
requested_cols   <- setdiff(trimws(readLines(FIELD_LIST)), c("", "eid"))
requested_fields <- unique(as.integer(sub("^p(\\d+).*$", "\\1", requested_cols)))
requested_fields <- sort(requested_fields[!is.na(requested_fields)])
only_list <- setdiff(requested_fields, field_lookup$field_id)
only_reg  <- setdiff(field_lookup$field_id, requested_fields)
if (length(only_list) || length(only_reg)) {
  stop("field_list.txt and the field registry disagree.",
       if (length(only_list)) paste0("\n  in field_list.txt only: ",
                                     paste(sort(only_list), collapse = ", ")),
       if (length(only_reg))  paste0("\n  in the registry only:   ",
                                     paste(sort(only_reg), collapse = ", ")))
}
cat("  field_list.txt agrees with the registry:", length(requested_fields),
    "field ids across", length(requested_cols), "columns\n")
cat("  Domains:", paste(unique(field_lookup$domain), collapse = ", "), "\n")

# ============================================================================
# 3. COLUMN PARSING — Map CSV columns to field IDs
# ============================================================================

cat("\n--- Step 3: Column name parsing ---\n")

# Parse column name pattern: p{field_id}_i{instance}_a{array}  or  p{field_id}
parse_colname <- function(cn) {
  m <- str_match(cn, "^p(\\d+)(?:_i(\\d+))?(?:_a(\\d+))?$")
  if (is.na(m[1, 1])) return(NULL)
  data.table(
    colname  = cn,
    field_id = as.integer(m[1, 2]),
    instance = ifelse(is.na(m[1, 3]), 0L, as.integer(m[1, 3])),
    array_ix = ifelse(is.na(m[1, 4]), NA_integer_, as.integer(m[1, 4]))
  )
}

col_map <- rbindlist(lapply(setdiff(names(raw), "eid"), parse_colname))
col_map <- merge(col_map, field_lookup, by = "field_id", all.x = TRUE)

# Keep only EXPORT fields (those in our registry) at instance 0
export_cols <- col_map[!is.na(domain) & instance == 0]
cat("  Total CSV columns:", ncol(raw) - 1, "\n")
cat("  Matched to EXPORT fields:", nrow(export_cols), "\n")
cat("  Unique field IDs:", length(unique(export_cols$field_id)), "\n")

# Fields in registry but NOT found in data
missing_in_data <- field_lookup[!field_id %in% col_map$field_id]
if (nrow(missing_in_data) > 0) {
  cat("  WARNING — Registered fields not found in data:\n")
  print(missing_in_data[, .(field_id, domain, var_label)])
}

# ============================================================================
# 4. UNIVERSAL MISSING VALUE CLEANING
# ============================================================================

cat("\n--- Step 4: Universal missing value cleaning ---\n")

# REPLACE mode text labels for missing/refusal → NA
MISSING_LABELS <- c(
  "Do not know",
  "Prefer not to answer",
  "Not sure - had a hysterectomy"  # Field 2724 specific
)

# Track cleaning stats
clean_stats <- data.table(
  field_id = integer(), var_label = character(), domain = character(),
  n_dk = integer(), n_pnta = integer(), n_empty = integer(),
  n_valid = integer(), pct_missing_final = numeric()
)

# Apply to all character columns that are EXPORT fields
char_export <- export_cols[col_map[, .(colname, field_id)], 
                           on = "field_id", nomatch = 0]

for (cn in intersect(export_cols$colname, names(raw))) {
  if (is.character(raw[[cn]])) {
    # Count before cleaning
    n_dk   <- sum(raw[[cn]] %in% c("Do not know"), na.rm = TRUE)
    n_pnta <- sum(raw[[cn]] %in% c("Prefer not to answer",
                                   "Not sure - had a hysterectomy"), na.rm = TRUE)
    n_empty <- sum(raw[[cn]] == "" | is.na(raw[[cn]]), na.rm = TRUE)
    
    # Clean
    raw[get(cn) %in% MISSING_LABELS, (cn) := NA_character_]
    raw[get(cn) == "", (cn) := NA_character_]
    
    n_valid <- sum(!is.na(raw[[cn]]))
    
    info <- export_cols[colname == cn][1]
    clean_stats <- rbind(clean_stats, data.table(
      field_id = info$field_id, var_label = info$var_label,
      domain = info$domain,
      n_dk = n_dk, n_pnta = n_pnta, n_empty = n_empty,
      n_valid = n_valid,
      pct_missing_final = round(100 * (1 - n_valid / nrow(raw)), 1)
    ))
  }
}

cat("  Character columns cleaned:", nrow(clean_stats), "\n")
cat("  Median final missing rate:", 
    round(median(clean_stats$pct_missing_final, na.rm = TRUE), 1), "%\n")

# ============================================================================
# 5. MULTI-ARRAY AGGREGATION
# ============================================================================

cat("\n--- Step 5: Multi-array field aggregation ---\n")

# Blood pressure: mean of two automated readings (a0, a1)
aggregate_mean <- function(dt, field_id, var_name) {
  a0 <- paste0("p", field_id, "_i0_a0")
  a1 <- paste0("p", field_id, "_i0_a1")
  if (all(c(a0, a1) %in% names(dt))) {
    dt[, (var_name) := rowMeans(.SD, na.rm = TRUE), .SDcols = c(a0, a1)]
    # If both NA, result should be NA (not NaN)
    dt[is.nan(get(var_name)), (var_name) := NA_real_]
    cat("  Aggregated", var_name, "= mean(", a0, ",", a1, ")\n")
  } else {
    # Fallback: use whichever exists
    if (a0 %in% names(dt)) {
      dt[, (var_name) := get(a0)]
      cat("  Fallback", var_name, "=", a0, "(a1 not found)\n")
    }
  }
  invisible(dt)
}

raw <- aggregate_mean(raw, 4079, "dbp_mean")
raw <- aggregate_mean(raw, 4080, "sbp_mean")
raw <- aggregate_mean(raw, 95,   "pulse_mean")

# Lung function: best of 3 blows (max of a0, a1, a2)
aggregate_best <- function(dt, field_id, var_name) {
  cols <- paste0("p", field_id, "_i0_a", 0:2)
  cols <- intersect(cols, names(dt))
  if (length(cols) > 0) {
    dt[, (var_name) := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = cols]
    dt[is.infinite(get(var_name)), (var_name) := NA_real_]
    cat("  Aggregated", var_name, "= max(", paste(cols, collapse = ", "), ")\n")
  }
  invisible(dt)
}

raw <- aggregate_best(raw, 3062, "fev1_best")
raw <- aggregate_best(raw, 3063, "fvc_best")

# Derived: FEV1/FVC ratio
if (all(c("fev1_best", "fvc_best") %in% names(raw))) {
  raw[, fev1_fvc_ratio := fifelse(fvc_best > 0, fev1_best / fvc_best, NA_real_)]
  cat("  Derived: fev1_fvc_ratio = fev1_best / fvc_best\n")
}

# Grip strength: max of left and right
if (all(c("p46_i0", "p47_i0") %in% names(raw))) {
  raw[, grip_max := pmax(p46_i0, p47_i0, na.rm = TRUE)]
  raw[is.infinite(grip_max), grip_max := NA_real_]
  cat("  Derived: grip_max = max(left, right)\n")
}

# ============================================================================
# 6. MULTI-RESPONSE PIPE-SEPARATED FIELD DECOMPOSITION
# ============================================================================

cat("\n--- Step 6: Multi-response field decomposition ---\n")

# Generic pipe-separated decomposition into binary indicators
decompose_multiresponse <- function(dt, colname, prefix, 
                                    items_to_extract = NULL,
                                    none_label = "None of the above") {
  if (!colname %in% names(dt)) {
    cat("  SKIP", colname, "- not in data\n")
    return(dt)
  }
  
  vals <- dt[[colname]]
  vals[is.na(vals)] <- ""
  
  # Split all responses and find unique items
  all_items <- unique(unlist(strsplit(vals, "\\|")))
  all_items <- trimws(all_items)
  all_items <- all_items[all_items != "" & 
                           !all_items %in% c(none_label, "Do not know",
                                             "Prefer not to answer")]
  
  if (!is.null(items_to_extract)) {
    all_items <- intersect(items_to_extract, all_items)
  }
  
  # Create binary indicators
  for (item in all_items) {
    safe_name <- paste0(prefix, "_",
                        gsub("[^a-zA-Z0-9]", "_", tolower(item)))
    safe_name <- gsub("_+", "_", safe_name)
    safe_name <- gsub("_$", "", safe_name)
    dt[, (safe_name) := as.integer(grepl(item, get(colname), fixed = TRUE))]
  }
  
  # Count of endorsed items
  dt[, paste0(prefix, "_count") := 
       str_count(get(colname), "\\|") + 
       fifelse(get(colname) != "" & !get(colname) %in% 
                 c(none_label, "Do not know", "Prefer not to answer"), 1L, 0L)]
  dt[get(colname) == "" | is.na(get(colname)), 
     paste0(prefix, "_count") := NA_integer_]
  dt[get(colname) %in% c(none_label), paste0(prefix, "_count") := 0L]
  
  cat("  Decomposed", colname, "→", length(all_items), "binary +",
      paste0(prefix, "_count"), "\n")
  invisible(dt)
}

# ----------------------------------------------------------------------------
# Code → label substitution for multi-response fields
# ----------------------------------------------------------------------------
# UKB Table Exporter REPLACE mode substitutes labels for single-valued
# categorical fields, but NOT for pipe-delimited multi-value fields. Those
# remain as numeric codes (e.g., "1|3", "-7"). We substitute them here so
# that downstream decompose_multiresponse sees text labels and produces
# interpretable column names (e.g., stressor_death_of_a_close_relative
# instead of stressor_3), AND correctly filters out none-of-the-above
# codes (-7, "I eat all of the above" etc.).
#
# Code mappings taken from UK Biobank Data-Codings:
#   Field 6138: DC 100305  qualifications
#   Field 6142: DC 100306  employment status
#   Field 6144: DC 100516  never eat
#   Field 6145: DC 100697  life stressors
#   Field 6150: DC 100605  vascular/heart problems
#   Field 6159: DC 100692  pain types last month
# ----------------------------------------------------------------------------

multi_response_labels <- list(
  "p6138_i0" = c(
    "1"  = "College or University degree",
    "2"  = "A levels/AS levels or equivalent",
    "3"  = "O levels/GCSEs or equivalent",
    "4"  = "CSEs or equivalent",
    "5"  = "NVQ or HND or HNC or equivalent",
    "6"  = "Other professional qualifications eg: nursing, teaching",
    "-7" = "None of the above",
    "-3" = "Prefer not to answer"
  ),
  "p6142_i0" = c(
    "1"  = "In paid employment or self-employed",
    "2"  = "Retired",
    "3"  = "Looking after home and/or family",
    "4"  = "Unable to work because of sickness or disability",
    "5"  = "Unemployed",
    "6"  = "Doing unpaid or voluntary work",
    "7"  = "Full or part-time student",
    "-7" = "None of the above",
    "-3" = "Prefer not to answer"
  ),
  "p6144_i0" = c(
    "1"  = "Eggs or foods containing eggs",
    "2"  = "Dairy products",
    "3"  = "Wheat products",
    "4"  = "Sugar or foods/drinks containing sugar",
    "5"  = "I eat all of the above",
    "-3" = "Prefer not to answer"
  ),
  "p6145_i0" = c(
    "1"  = "Serious illness, injury or assault to yourself",
    "2"  = "Serious illness, injury or assault of a close relative",
    "3"  = "Death of a close relative",
    "4"  = "Death of a spouse or partner",
    "5"  = "Marital separation/divorce",
    "6"  = "Financial difficulties",
    "-7" = "None of the above",
    "-3" = "Prefer not to answer"
  ),
  "p6150_i0" = c(
    "1"  = "Heart attack",
    "2"  = "Angina",
    "3"  = "Stroke",
    "4"  = "High blood pressure",
    "-7" = "None of the above",
    "-3" = "Prefer not to answer"
  ),
  "p6159_i0" = c(
    "1"  = "Headache",
    "2"  = "Facial pain",
    "3"  = "Neck or shoulder pain",
    "4"  = "Back pain",
    "5"  = "Stomach or abdominal pain",
    "6"  = "Hip pain",
    "7"  = "Knee pain",
    "8"  = "Pain all over the body",
    "-7" = "None of the above",
    "-3" = "Prefer not to answer"
  )
)

# Helper: substitute numeric codes with text labels in a pipe-delimited column
# Skip silently if the column already contains text (idempotent).
substitute_multi_response_codes <- function(dt, colname, label_map) {
  if (!(colname %in% names(dt))) {
    cat("  SKIP", colname, "- not in data\n")
    return(invisible(dt))
  }
  raw_vals <- as.character(dt[[colname]])
  non_na <- raw_vals[!is.na(raw_vals) & raw_vals != ""]
  if (length(non_na) == 0) {
    cat("  SKIP", colname, "- all missing\n")
    return(invisible(dt))
  }
  # Heuristic: if first non-NA value has any letter, it's already text
  if (grepl("[A-Za-z]", non_na[1])) {
    cat("  ", colname, "already text-labeled, no substitution needed\n")
    return(invisible(dt))
  }
  # Vectorised split → map → rejoin
  split_lists <- strsplit(raw_vals, "|", fixed = TRUE)
  sub_lists <- lapply(split_lists, function(parts) {
    if (length(parts) == 0) return(character(0))
    mapped <- label_map[parts]
    mapped[is.na(mapped)] <- parts[is.na(mapped)]  # keep unknown codes as-is
    unname(mapped)
  })
  sub_vals <- vapply(sub_lists, function(x) paste(x, collapse = "|"),
                     character(1))
  sub_vals[is.na(raw_vals)] <- NA_character_
  dt[, (colname) := sub_vals]
  cat("  ", colname, ": substituted", length(label_map), "codes → labels\n")
  invisible(dt)
}

# Apply substitution to all multi-response fields
cat("\n  Substituting numeric codes → text labels for multi-response fields:\n")
for (cn in names(multi_response_labels)) {
  substitute_multi_response_codes(raw, cn, multi_response_labels[[cn]])
}
cat("\n")

# 6145: Life stressors (Illness, injury, bereavement, stress in last 2 years)
raw <- decompose_multiresponse(raw, "p6145_i0", "stressor")

# 6159: Pain types experienced in last month
raw <- decompose_multiresponse(raw, "p6159_i0", "pain_month")

# 6144: Never eat (eggs, dairy, wheat, sugar)
raw <- decompose_multiresponse(raw, "p6144_i0", "never_eat",
                               none_label = "I eat all of the above")

# 6150: Vascular/heart problems — extract specific conditions
raw <- decompose_multiresponse(raw, "p6150_i0", "vascular")

# 6138: Qualifications — binary indicators
raw <- decompose_multiresponse(raw, "p6138_i0", "qualification")

# 6142: Employment status — binary indicators  
raw <- decompose_multiresponse(raw, "p6142_i0", "employment")

# ============================================================================
# 7. DOMAIN-SPECIFIC RECODING
# ============================================================================

cat("\n--- Step 7: Domain-specific recoding ---\n")

# ---- 7a. Binary Yes/No fields ----
# Pattern: "Yes" → 1, "No" → 0, rest → NA (already cleaned in Step 4)

binary_yesno_fields <- c(
  "p1210_i0",   # snoring
  "p2020_i0",   # loneliness
  "p2188_i0",   # long-standing illness
  "p2247_i0",   # hearing difficulty
  "p2443_i0",   # diabetes self-report
  "p2956_i0",   # general pain 3+ months
  "p3404_i0",   # neck pain 3+ months
  "p3414_i0",   # hip pain 3+ months
  "p3571_i0",   # back pain 3+ months
  "p3741_i0",   # stomach pain 3+ months
  "p3773_i0",   # knee pain 3+ months
  "p3799_i0",   # headache 3+ months
  "p4067_i0",   # facial pain 3+ months
  "p4501_i0",   # family non-accidental death
  "p4598_i0",   # ever depressed ≥2 weeks
  "p4631_i0",   # ever unenthusiastic ≥2 weeks
  "p1920_i0"    # mood swings
)

for (cn in intersect(binary_yesno_fields, names(raw))) {
  new_name <- paste0(cn, "_bin")
  raw[, (new_name) := fcase(
    get(cn) == "Yes", 1L,
    get(cn) == "No",  0L,
    default = NA_integer_
  )]
}

# Cancer self-report (special wording)
if ("p2453_i0" %in% names(raw)) {
  raw[, p2453_i0_bin := fcase(
    grepl("^Yes", p2453_i0), 1L,
    p2453_i0 == "No",        0L,
    default = NA_integer_
  )]
}

cat("  Binary Yes/No recoded:", length(binary_yesno_fields), "fields\n")

# ---- 7b. Ordinal frequency fields ----

# Dietary frequency: "Never" to "Once or more daily" (9 levels)
diet_freq_fields <- c("p1329_i0", "p1339_i0", "p1349_i0", "p1359_i0",
                      "p1369_i0", "p1379_i0", "p1389_i0", "p1408_i0")

for (cn in intersect(diet_freq_fields, names(raw))) {
  new_name <- paste0(cn, "_ord")
  raw[, (new_name) := fcase(
    get(cn) == "Never",                 0L,
    get(cn) == "Less than once a week", 1L,
    get(cn) == "Once a week",           2L,
    get(cn) == "2-4 times a week",      3L,
    get(cn) == "5-6 times a week",      4L,
    get(cn) == "Once or more daily",    5L,
    default = NA_integer_
  )]
}
cat("  Dietary frequency ordinal:", length(diet_freq_fields), "fields\n")

# Salt added to food
if ("p1478_i0" %in% names(raw)) {
  raw[, salt_added_ord := fcase(
    p1478_i0 == "Never/rarely", 0L,
    p1478_i0 == "Sometimes",    1L,
    p1478_i0 == "Usually",      2L,
    p1478_i0 == "Always",       3L,
    default = NA_integer_
  )]
}

# ---- 7c. Numeric-as-text fields (portions/day, hours, etc.) ----
# These have numeric values stored as character due to "Less than one" option

numeric_text_fields <- c(
  "p1289_i0", "p1299_i0", "p1309_i0", "p1319_i0",  # vegetable/fruit portions
  "p1438_i0", "p1458_i0",                            # bread/cereal portions
  "p1488_i0", "p1498_i0", "p1528_i0",                # tea/coffee/water
  "p845_i0",                                          # age completed education
  "p1160_i0",                                         # sleep duration
  "p3680_i0",                                         # age last ate meat
  "p2714_i0",                                         # age menarche
  "p3581_i0"                                          # age menopause
)

for (cn in intersect(numeric_text_fields, names(raw))) {
  new_name <- paste0(cn, "_num")
  raw[, (new_name) := {
    x <- get(cn)
    x <- gsub("Less than one", "0.5", x, fixed = TRUE)
    suppressWarnings(as.numeric(x))
  }]
}
cat("  Numeric-as-text converted:", length(numeric_text_fields), "fields\n")

# ---- 7d. Ordinal categorical fields ----

# Overall health rating
if ("p2178_i0" %in% names(raw)) {
  raw[, health_rating_ord := fcase(
    p2178_i0 == "Excellent", 4L,
    p2178_i0 == "Good",      3L,
    p2178_i0 == "Fair",      2L,
    p2178_i0 == "Poor",      1L,
    default = NA_integer_
  )]
}

# Falls in last year
if ("p2296_i0" %in% names(raw)) {
  raw[, falls_ord := fcase(
    p2296_i0 == "No falls",            0L,
    p2296_i0 == "Only one fall",       1L,
    p2296_i0 == "More than one fall",  2L,
    default = NA_integer_
  )]
}

# Weight change compared with 1 year ago
if ("p2306_i0" %in% names(raw)) {
  raw[, weight_change_cat := fcase(
    p2306_i0 == "No - weigh about the same", 0L,
    p2306_i0 == "Yes - gained weight",       1L,
    p2306_i0 == "Yes - lost weight",         -1L,
    default = NA_integer_
  )]
}

# Getting up in morning
if ("p1170_i0" %in% names(raw)) {
  raw[, getting_up_ord := fcase(
    p1170_i0 == "Very easy",      4L,
    p1170_i0 == "Fairly easy",    3L,
    p1170_i0 == "Not very easy",  2L,
    p1170_i0 == "Not at all easy", 1L,
    default = NA_integer_
  )]
}

# Chronotype
if ("p1180_i0" %in% names(raw)) {
  raw[, chronotype_ord := fcase(
    p1180_i0 == "Definitely a 'morning' person",              2L,
    p1180_i0 == "More a 'morning' than 'evening' person",     1L,
    p1180_i0 == "More an 'evening' than a 'morning' person", -1L,
    p1180_i0 == "Definitely an 'evening' person",            -2L,
    default = NA_integer_
  )]
}

# Nap during day / Daytime dozing / Insomnia
ordinal_never_usually <- c("p1190_i0", "p1200_i0", "p1220_i0")
for (cn in intersect(ordinal_never_usually, names(raw))) {
  new_name <- paste0(cn, "_ord")
  raw[, (new_name) := fcase(
    get(cn) == "Never/rarely", 0L,
    get(cn) == "Sometimes",    1L,
    get(cn) == "Often",        2L,
    get(cn) == "Usually",      3L,
    default = NA_integer_
  )]
}

# Depressed mood / Unenthusiasm frequency (PHQ-2 style)
phq2_fields <- c("p2050_i0", "p2060_i0")
for (cn in intersect(phq2_fields, names(raw))) {
  new_name <- paste0(cn, "_ord")
  raw[, (new_name) := fcase(
    get(cn) == "Not at all",            0L,
    get(cn) == "Several days",          1L,
    get(cn) == "More than half the days", 2L,
    get(cn) == "Nearly every day",      3L,
    default = NA_integer_
  )]
}

# Able to confide / Friend visits frequency
freq_social_fields <- c("p1031_i0", "p2110_i0")
for (cn in intersect(freq_social_fields, names(raw))) {
  new_name <- paste0(cn, "_ord")
  raw[, (new_name) := fcase(
    get(cn) == "No friends/family",     0L,
    get(cn) == "Never or almost never", 0L,
    get(cn) == "Once every few months", 1L,
    get(cn) == "About once a month",    2L,
    get(cn) == "About once a week",     3L,
    get(cn) == "2-4 times a week",      4L,
    get(cn) == "Almost daily",          5L,
    default = NA_integer_
  )]
}

# Household income
if ("p738_i0" %in% names(raw)) {
  raw[, income_ord := fcase(
    p738_i0 == "Less than 18,000",    1L,
    p738_i0 == "18,000 to 30,999",    2L,
    p738_i0 == "31,000 to 51,999",    3L,
    p738_i0 == "52,000 to 100,000",   4L,
    p738_i0 == "Greater than 100,000", 5L,
    default = NA_integer_
  )]
}

# Diet variation
if ("p1548_i0" %in% names(raw)) {
  raw[, diet_variation_ord := fcase(
    p1548_i0 == "Never/rarely", 0L,
    p1548_i0 == "Sometimes",    1L,
    p1548_i0 == "Often",        2L,
    default = NA_integer_
  )]
}

# Hot drink temperature
if ("p1518_i0" %in% names(raw)) {
  raw[, hot_drink_temp_ord := fcase(
    p1518_i0 == "Do not drink hot drinks", 0L,
    p1518_i0 == "Warm",                    1L,
    p1518_i0 == "Hot",                     2L,
    p1518_i0 == "Very hot",                3L,
    default = NA_integer_
  )]
}

# Tinnitus (ordinal)
if ("p4803_i0" %in% names(raw)) {
  raw[, tinnitus_ord := fcase(
    p4803_i0 == "No, never",                               0L,
    p4803_i0 == "Yes, but not now, but have in the past",  1L,
    p4803_i0 == "Yes, now some of the time",               2L,
    p4803_i0 == "Yes, now most of the time",               3L,
    p4803_i0 == "Yes, now all of the time",                4L,
    default = NA_integer_
  )]
}

# Menopause (binary for female subsample)
if ("p2724_i0" %in% names(raw)) {
  raw[, menopause_bin := fcase(
    p2724_i0 == "Yes", 1L,
    p2724_i0 == "No",  0L,
    default = NA_integer_
  )]
}

# Prospective memory
if ("p20018_i0" %in% names(raw)) {
  raw[, prospective_memory_ord := fcase(
    p20018_i0 == "Correct recall on first attempt",  2L,
    p20018_i0 == "Correct recall on second attempt",  1L,
    p20018_i0 == "Instruction not recalled, either skipped or incorrect", 0L,
    default = NA_integer_
  )]
}

# Body size at age 10
if ("p1687_i0" %in% names(raw)) {
  raw[, body_size_10_ord := fcase(
    p1687_i0 == "Thinner",       -1L,
    p1687_i0 == "About average",  0L,
    p1687_i0 == "Plumper",        1L,
    default = NA_integer_
  )]
}

# Handedness
if ("p1707_i0" %in% names(raw)) {
  raw[, handedness_cat := fcase(
    p1707_i0 == "Right-handed",                          "right",
    p1707_i0 == "Left-handed",                           "left",
    p1707_i0 == "Use both right and left hands equally", "ambidextrous",
    default = NA_character_
  )]
  cat("  handedness_cat distribution:\n")
  print(table(raw$handedness_cat, useNA = "ifany"))
}

# NOTE: Categorical dietary/handedness variables are dummy-coded in Step 10
# via decompose_categorical(), not kept as single character columns.

# Diet change in last 5 years
if ("p1538_i0" %in% names(raw)) {
  raw[, diet_change_bin := fcase(
    p1538_i0 == "No",                             0L,
    p1538_i0 == "Yes, because of illness",        1L,
    p1538_i0 == "Yes, because of other reasons",  1L,
    default = NA_integer_
  )]
  raw[, diet_change_illness := as.integer(p1538_i0 == "Yes, because of illness")]
}

cat("  Domain-specific recoding completed\n")

# ============================================================================
# 8. CONTINUOUS VARIABLE QUALITY CONTROL
# ============================================================================

cat("\n--- Step 8: Continuous variable QC ---\n")

# Biologically plausible ranges for winsorization/exclusion
# Values outside these ranges → NA
bio_ranges <- list(
  # Anthropometric
  p50_i0     = c(130, 210),    # height (cm)
  p21002_i0  = c(30, 250),     # weight (kg)
  p48_i0     = c(50, 180),     # waist (cm)
  p49_i0     = c(50, 180),     # hip (cm)
  p46_i0     = c(0, 90),       # grip left (kg)
  p47_i0     = c(0, 90),       # grip right (kg)
  sbp_mean   = c(70, 250),     # SBP (mmHg)
  dbp_mean   = c(30, 150),     # DBP (mmHg)
  pulse_mean = c(25, 200),     # pulse (bpm)
  p23099_i0  = c(1, 75),       # body fat %
  
  # Biochemistry — standard clinical ranges (generous)
  p30600_i0  = c(15, 60),      # albumin (g/L)
  p30690_i0  = c(1, 20),       # total cholesterol (mmol/L)
  p30700_i0  = c(10, 700),     # creatinine (umol/L)
  p30710_i0  = c(0, 200),      # CRP (mg/L)
  p30740_i0  = c(1, 40),       # glucose (mmol/L)
  p30750_i0  = c(10, 200),     # HbA1c (mmol/mol)
  p30870_i0  = c(0, 30),       # triglycerides (mmol/L)
  p30880_i0  = c(50, 900),     # urate (umol/L)
  p30890_i0  = c(0, 400),      # vitamin D (nmol/L)
  
  # Lung function
  fev1_best  = c(0.3, 7),      # FEV1 (L)
  fvc_best   = c(0.5, 9),      # FVC (L)
  
  # Derived
  fev1_fvc_ratio = c(0.15, 1.1),
  
  # Other continuous
  p20161_i0  = c(0, 200),      # pack-years
  p22040_i0  = c(0, 20000),    # MET minutes/week
  p24006     = c(0, 100),      # PM2.5 (ug/m3)
  p20023_i0  = c(50, 5000),    # reaction time (ms)
  p20127_i0  = c(0, 12)        # neuroticism score
)

n_outliers_total <- 0
for (vname in names(bio_ranges)) {
  if (vname %in% names(raw)) {
    rng <- bio_ranges[[vname]]
    n_out <- sum(raw[[vname]] < rng[1] | raw[[vname]] > rng[2], na.rm = TRUE)
    if (n_out > 0) {
      raw[get(vname) < rng[1] | get(vname) > rng[2], (vname) := NA]
      n_outliers_total <- n_outliers_total + n_out
      cat("  ", vname, ": removed", n_out, "outliers outside [",
          rng[1], ",", rng[2], "]\n")
    }
  }
}
cat("  Total outlier values removed:", n_outliers_total, "\n")

# ============================================================================
# 9. DERIVED VARIABLES
# ============================================================================

cat("\n--- Step 9: Derived variables ---\n")

# WHR (waist-hip ratio)
if (all(c("p48_i0", "p49_i0") %in% names(raw))) {
  raw[, whr := fifelse(p49_i0 > 0, p48_i0 / p49_i0, NA_real_)]
  cat("  Derived: whr = waist / hip\n")
}

# PHQ-2 composite (sum of depressed mood + unenthusiasm)
if (all(c("p2050_i0_ord", "p2060_i0_ord") %in% names(raw))) {
  raw[, phq2_score := p2050_i0_ord + p2060_i0_ord]
  cat("  Derived: phq2_score = depressed_mood + unenthusiasm (0–6)\n")
}

# Pain count (number of pain sites with 3+ months duration)
pain_bin_cols <- paste0(c("p2956_i0", "p3404_i0", "p3414_i0", "p3571_i0",
                          "p3741_i0", "p3773_i0", "p3799_i0", "p4067_i0"),
                        "_bin")
pain_bin_exist <- intersect(pain_bin_cols, names(raw))
if (length(pain_bin_exist) > 0) {
  raw[, pain_site_count := rowSums(.SD, na.rm = TRUE), .SDcols = pain_bin_exist]
  # If ALL pain sites are NA, set count to NA (cannot use .SD in i-expression)
  all_pain_na <- rowSums(!is.na(raw[, ..pain_bin_exist])) == 0
  raw[all_pain_na, pain_site_count := NA_integer_]
  cat("  Derived: pain_site_count from", length(pain_bin_exist), "sites\n")
}

# SRT mean (average of left and right ear)
if (all(c("p20019_i0", "p20021_i0") %in% names(raw))) {
  raw[, srt_mean := rowMeans(.SD, na.rm = TRUE),
      .SDcols = c("p20019_i0", "p20021_i0")]
  raw[is.nan(srt_mean), srt_mean := NA_real_]
  cat("  Derived: srt_mean = mean(SRT_left, SRT_right)\n")
}

# ============================================================================
# 10. BUILD ANALYSIS-READY EXPOSURE MATRIX
# ============================================================================

cat("\n--- Step 10: Building exposure matrix ---\n")

# Collect all processed exposure columns
# Strategy: for each field, select the PROCESSED version preferentially

# Initialize variable dictionary
var_dict <- data.table(
  var_name     = character(),
  source_field = character(),
  domain       = character(),
  var_type     = character(),  # continuous, ordinal, binary, categorical
  female_only  = logical(),
  description  = character()
)

# Helper to register variables
register_var <- function(var_name, source_field, domain, var_type,
                         female_only = FALSE, description = "") {
  var_dict <<- rbind(var_dict, data.table(
    var_name = var_name, source_field = source_field,
    domain = domain, var_type = var_type,
    female_only = female_only, description = description
  ))
}

# ---- Continuous exposures (direct from raw, after QC) ----
continuous_direct <- list(
  # Anthropometric
  list("height",          "p50_i0",      "Anthropometric & physiological", "Standing height (cm)"),
  list("weight",          "p21002_i0",   "Anthropometric & physiological", "Weight (kg)"),
  list("waist_circ",      "p48_i0",      "Anthropometric & physiological", "Waist circumference (cm)"),
  list("hip_circ",        "p49_i0",      "Anthropometric & physiological", "Hip circumference (cm)"),
  list("grip_left",       "p46_i0",      "Anthropometric & physiological", "Hand grip strength left (kg)"),
  list("grip_right",      "p47_i0",      "Anthropometric & physiological", "Hand grip strength right (kg)"),
  list("body_fat_pct",    "p23099_i0",   "Anthropometric & physiological", "Body fat percentage"),
  list("whole_body_fat",  "p23100_i0",   "Anthropometric & physiological", "Whole body fat mass (kg)"),
  list("fat_free_mass",   "p23101_i0",   "Anthropometric & physiological", "Whole body fat-free mass (kg)"),
  list("trunk_fat_pct",   "p23127_i0",   "Anthropometric & physiological", "Trunk fat percentage"),
  
  # Derived anthropometric
  list("sbp",             "sbp_mean",    "Anthropometric & physiological", "Systolic BP mean of 2 readings (mmHg)"),
  list("dbp",             "dbp_mean",    "Anthropometric & physiological", "Diastolic BP mean of 2 readings (mmHg)"),
  list("pulse",           "pulse_mean",  "Anthropometric & physiological", "Pulse rate mean of 2 readings (bpm)"),
  list("grip_max",        "grip_max",    "Anthropometric & physiological", "Max grip strength L/R (kg)"),
  list("whr",             "whr",         "Anthropometric & physiological", "Waist-hip ratio"),
  list("fev1",            "fev1_best",   "Anthropometric & physiological", "FEV1 best of 3 blows (L)"),
  list("fvc",             "fvc_best",    "Anthropometric & physiological", "FVC best of 3 blows (L)"),
  list("fev1_fvc",        "fev1_fvc_ratio", "Anthropometric & physiological", "FEV1/FVC ratio"),
  
  # Lifestyle
  list("pack_years",      "p20161_i0",   "Lifestyle & sleep", "Pack-years of smoking"),
  list("physical_met",    "p22040_i0",   "Lifestyle & sleep", "Physical activity MET min/week"),
  
  # Mental health
  list("neuroticism",     "p20127_i0",   "General & mental health", "Neuroticism score (0–12)"),
  
  # Sensory
  list("srt_left",        "p20019_i0",   "Sensory function & pain", "Speech reception threshold left (dB)"),
  list("srt_right",       "p20021_i0",   "Sensory function & pain", "Speech reception threshold right (dB)"),
  list("srt_mean",        "srt_mean",    "Sensory function & pain", "Speech reception threshold mean (dB)"),
  
  # Cognitive
  list("fluid_intel",     "p20016_i0",   "Clinical screening & treatment", "Fluid intelligence score"),
  list("reaction_time",   "p20023_i0",   "Clinical screening & treatment", "Mean reaction time (ms)"),
  
  # General health
  list("n_medications",   "p137_i0",     "General & mental health", "Number of medications taken"),
  
  # Environment
  list("pm25",            "p24006",      "Environmental exposures", "PM2.5 air pollution 2010 (ug/m3)"),
  
  # Female-only
  list("age_menarche",    "p2714_i0_num", "Reproductive (female-only)", "Age at menarche (years)"),
  list("age_menopause",   "p3581_i0_num", "Reproductive (female-only)", "Age at menopause (years)")
)

# Biochemistry
biochem_vars <- list(
  list("albumin",      "p30600_i0"), list("alp",           "p30610_i0"),
  list("alt",          "p30620_i0"), list("apoa",          "p30630_i0"),
  list("apob",         "p30640_i0"), list("ast",           "p30650_i0"),
  list("direct_bili",  "p30660_i0"), list("urea",          "p30670_i0"),
  list("calcium",      "p30680_i0"), list("cholesterol",   "p30690_i0"),
  list("creatinine",   "p30700_i0"), list("crp",           "p30710_i0"),
  list("cystatin_c",   "p30720_i0"), list("ggt",           "p30730_i0"),
  list("glucose",      "p30740_i0"), list("hba1c",         "p30750_i0"),
  list("hdl",          "p30760_i0"), list("igf1",          "p30770_i0"),
  list("ldl",          "p30780_i0"), list("lpa",           "p30790_i0"),
  list("oestradiol",   "p30800_i0"), list("phosphate",     "p30810_i0"),
  list("rheumatoid_f", "p30820_i0"), list("shbg",          "p30830_i0"),
  list("total_bili",   "p30840_i0"), list("testosterone",  "p30850_i0"),
  list("total_protein","p30860_i0"), list("triglycerides", "p30870_i0"),
  list("urate",        "p30880_i0"), list("vitamin_d",     "p30890_i0")
)

# Haematology
haem_vars <- list(
  list("wbc",           "p30000_i0"), list("rbc",            "p30010_i0"),
  list("haemoglobin",   "p30020_i0"), list("haematocrit",    "p30030_i0"),
  list("mcv",           "p30040_i0"), list("mch",            "p30050_i0"),
  list("mchc",          "p30060_i0"), list("rdw",            "p30070_i0"),
  list("platelet",      "p30080_i0"), list("platelet_crit",  "p30090_i0"),
  list("mpv",           "p30100_i0"), list("platelet_dw",    "p30110_i0"),
  list("lymphocyte",    "p30120_i0"), list("monocyte",       "p30130_i0"),
  list("neutrophil",    "p30140_i0"), list("eosinophil",     "p30150_i0"),
  list("basophil",      "p30160_i0"), list("nrbc",           "p30170_i0"),
  list("lymphocyte_pct","p30180_i0"), list("monocyte_pct",   "p30190_i0"),
  list("neutrophil_pct","p30200_i0"), list("eosinophil_pct", "p30210_i0"),
  list("basophil_pct",  "p30220_i0"), list("nrbc_pct",       "p30230_i0"),
  list("retic_pct",     "p30240_i0"), list("retic_count",    "p30250_i0"),
  list("mean_retic_vol","p30260_i0"), list("mean_spher_vol", "p30270_i0"),
  list("immature_retic","p30280_i0"), list("hls_retic_pct",  "p30290_i0"),
  list("hls_retic",     "p30300_i0")
)

# ---- Collect all exposure columns ----
exposure_cols <- c("eid")
missing_cols  <- character()

# Helper: add column safely
# NOTE: uses raw[[source_col]] (base R extraction) instead of get(source_col)
# to avoid a data.table j-scope evaluation issue where get() with a variable
# argument can return incorrect values for character derived columns.
add_exposure <- function(new_name, source_col, domain, var_type,
                         female_only = FALSE, desc = "") {
  if (source_col %in% names(raw)) {
    raw[, (new_name) := raw[[source_col]]]
    register_var(new_name, source_col, domain, var_type, female_only, desc)
    exposure_cols <<- c(exposure_cols, new_name)
  } else {
    missing_cols <<- c(missing_cols, paste0(new_name, " (", source_col, ")"))
  }
}

# Register continuous direct
for (v in continuous_direct) {
  fo <- grepl("female|A5", v[[3]])
  add_exposure(v[[1]], v[[2]], v[[3]], "continuous", fo, v[[4]])
}

# Register biochemistry
for (v in biochem_vars) {
  add_exposure(v[[1]], v[[2]], "Blood biochemistry & haematology", "continuous", FALSE, 
               paste("Serum", v[[1]]))
}

# Register haematology
for (v in haem_vars) {
  add_exposure(v[[1]], v[[2]], "Blood biochemistry & haematology", "continuous", FALSE,
               paste("Blood", v[[1]]))
}

# Register ordinal/binary processed variables
ordinal_binary_vars <- list(
  # Binary
  list("snoring",         "p1210_i0_bin",   "Lifestyle & sleep",        "binary"),
  list("loneliness",      "p2020_i0_bin",   "General & mental health", "binary"),
  list("longstanding_ill","p2188_i0_bin",   "General & mental health","binary"),
  list("hearing_diff",    "p2247_i0_bin",   "Sensory function & pain",        "binary"),
  list("diabetes_sr",     "p2443_i0_bin",   "Clinical screening & treatment",        "binary"),
  list("cancer_sr",       "p2453_i0_bin",   "Clinical screening & treatment",        "binary"),
  list("general_pain_3mo","p2956_i0_bin",   "Sensory function & pain",        "binary"),
  list("neck_pain_3mo",   "p3404_i0_bin",   "Sensory function & pain",        "binary"),
  list("hip_pain_3mo",    "p3414_i0_bin",   "Sensory function & pain",        "binary"),
  list("back_pain_3mo",   "p3571_i0_bin",   "Sensory function & pain",        "binary"),
  list("stomach_pain_3mo","p3741_i0_bin",   "Sensory function & pain",        "binary"),
  list("knee_pain_3mo",   "p3773_i0_bin",   "Sensory function & pain",        "binary"),
  list("headache_3mo",    "p3799_i0_bin",   "Sensory function & pain",        "binary"),
  list("facial_pain_3mo", "p4067_i0_bin",   "Sensory function & pain",        "binary"),
  list("fam_nonacc_death","p4501_i0_bin",   "Clinical screening & treatment",        "binary"),
  list("ever_depressed",  "p4598_i0_bin",   "General & mental health", "binary"),
  list("ever_unenthu",    "p4631_i0_bin",   "General & mental health", "binary"),
  list("mood_swings",     "p1920_i0_bin",   "General & mental health", "binary"),
  list("menopause",       "menopause_bin",  "Reproductive (female-only)", "binary"),
  list("diet_change",     "diet_change_bin","Dietary intake & preferences",        "binary"),
  list("diet_change_ill", "diet_change_illness", "Dietary intake & preferences",   "binary"),
  
  # Ordinal
  list("health_rating",   "health_rating_ord",   "General & mental health", "ordinal"),
  list("falls",           "falls_ord",           "General & mental health", "ordinal"),
  list("weight_change",   "weight_change_cat",   "General & mental health", "ordinal"),
  list("getting_up",      "getting_up_ord",      "Lifestyle & sleep",         "ordinal"),
  list("chronotype",      "chronotype_ord",      "Lifestyle & sleep",         "ordinal"),
  list("nap_daytime",     "p1190_i0_ord",        "Lifestyle & sleep",         "ordinal"),
  list("insomnia",        "p1200_i0_ord",        "Lifestyle & sleep",         "ordinal"),
  list("daytime_dozing",  "p1220_i0_ord",        "Lifestyle & sleep",         "ordinal"),
  list("depressed_mood",  "p2050_i0_ord",        "General & mental health",  "ordinal"),
  list("unenthusiasm",    "p2060_i0_ord",        "General & mental health",  "ordinal"),
  list("phq2_score",      "phq2_score",          "General & mental health",  "ordinal"),
  list("friend_visits",   "p1031_i0_ord",        "General & mental health",  "ordinal"),
  list("able_confide",    "p2110_i0_ord",        "General & mental health",  "ordinal"),
  list("income",          "income_ord",          "Demographics & SES",         "ordinal"),
  list("salt_added",      "salt_added_ord",      "Dietary intake & preferences",         "ordinal"),
  list("diet_variation",  "diet_variation_ord",  "Dietary intake & preferences",         "ordinal"),
  list("hot_drink_temp",  "hot_drink_temp_ord",  "Dietary intake & preferences",         "ordinal"),
  list("tinnitus",        "tinnitus_ord",        "Sensory function & pain",         "ordinal"),
  list("prosp_memory",    "prospective_memory_ord", "Clinical screening & treatment",      "ordinal"),
  list("body_size_10",    "body_size_10_ord",    "Anthropometric & physiological",         "ordinal"),
  list("pain_site_count", "pain_site_count",     "Sensory function & pain",         "ordinal")
)

for (v in ordinal_binary_vars) {
  fo <- grepl("menopause", v[[1]])
  add_exposure(v[[1]], v[[2]], v[[3]], v[[4]], fo, "")
}

# Register dietary frequency ordinal
diet_freq_names <- c("oily_fish_freq", "nonoily_fish_freq", "processed_meat_freq",
                     "poultry_freq", "beef_freq", "lamb_freq", "pork_freq",
                     "cheese_freq")
diet_freq_sources <- paste0(c("p1329_i0", "p1339_i0", "p1349_i0", "p1359_i0",
                              "p1369_i0", "p1379_i0", "p1389_i0", "p1408_i0"),
                            "_ord")
for (i in seq_along(diet_freq_names)) {
  add_exposure(diet_freq_names[i], diet_freq_sources[i], "Dietary intake & preferences", "ordinal", FALSE, "")
}

# Register numeric-as-text converted dietary amounts
diet_amount_map <- list(
  c("cooked_veg",   "p1289_i0_num"), c("salad_veg",   "p1299_i0_num"),
  c("fresh_fruit",  "p1309_i0_num"), c("dried_fruit",  "p1319_i0_num"),
  c("bread_intake", "p1438_i0_num"), c("cereal_intake", "p1458_i0_num"),
  c("tea_intake",   "p1488_i0_num"), c("coffee_intake", "p1498_i0_num"),
  c("water_intake", "p1528_i0_num"), c("age_last_meat", "p3680_i0_num")
)
for (v in diet_amount_map) {
  add_exposure(v[1], v[2], "Dietary intake & preferences", "continuous", FALSE, "")
}

# Register non-dietary numeric-as-text (correct domains)
add_exposure("sleep_hours",   "p1160_i0_num", "Lifestyle & sleep", "continuous", FALSE,
             "Sleep duration (hours)")
add_exposure("age_education", "p845_i0_num",  "Demographics & SES", "continuous", FALSE,
             "Age completed full-time education (years)")

# ----------------------------------------------------------------------------
# Register categorical variables as DUMMY-CODED binary indicators
# ----------------------------------------------------------------------------
# ExWAS tests each exposure independently via univariable (+ covariates)
# logistic regression. Multi-level categorical variables (e.g., milk_type
# with 5 levels) cannot be passed as a single character column to glm().
# We decompose them into binary dummy indicators at cleaning time, one per
# level, so that each level is tested independently against the outcome.
# This is the standard ExWAS approach — equivalent to testing "Semi-skimmed
# milk drinker vs all others" for each category.
#
# NOTE: ethnic_bg (Field 21000) EXCLUDED from exposure pool
#       — used as group stratification variable, not an ExWAS exposure
# ----------------------------------------------------------------------------

decompose_categorical <- function(dt, colname, prefix, domain,
                                  exclude_levels = NULL) {
  if (!(colname %in% names(dt))) {
    cat("    SKIP", colname, "- not in data\n")
    return(character(0))
  }
  x <- dt[[colname]]
  # Get unique non-NA levels
  lvls <- sort(unique(x[!is.na(x)]))
  # Exclude specified levels (e.g., "Other type of milk")
  if (!is.null(exclude_levels)) {
    lvls <- setdiff(lvls, exclude_levels)
  }
  if (length(lvls) == 0) {
    cat("    SKIP", colname, "- no valid levels\n")
    return(character(0))
  }
  
  created <- character()
  for (lv in lvls) {
    safe_lv <- gsub("[^a-zA-Z0-9]", "_", tolower(lv))
    safe_lv <- gsub("_+", "_", safe_lv)
    safe_lv <- gsub("^_|_$", "", safe_lv)
    dummy_name <- paste0(prefix, "_", safe_lv)
    
    dt[, (dummy_name) := fifelse(!is.na(x) & x == lv, 1L, fifelse(is.na(x), NA_integer_, 0L))]
    register_var(dummy_name, colname, domain, "binary", FALSE,
                 paste0(prefix, ": ", lv))
    exposure_cols <<- c(exposure_cols, dummy_name)
    created <- c(created, dummy_name)
  }
  cat("    ", colname, "→", prefix, ":", length(created), "dummies from",
      length(lvls), "levels\n")
  invisible(created)
}

cat("  Dummy-coding categorical variables:\n")

# dietary categorical (6 variables)
decompose_categorical(raw, "p1418_i0", "milk_type",    "Dietary intake & preferences")
decompose_categorical(raw, "p1428_i0", "spread_type",  "Dietary intake & preferences")
decompose_categorical(raw, "p1448_i0", "bread_type",   "Dietary intake & preferences")
decompose_categorical(raw, "p1468_i0", "cereal_type",  "Dietary intake & preferences")
decompose_categorical(raw, "p1508_i0", "coffee_type",  "Dietary intake & preferences")
decompose_categorical(raw, "p2654_i0", "spread_detail","Dietary intake & preferences")

# handedness handedness (3 levels: right/left/ambidextrous)
decompose_categorical(raw, "handedness_cat", "handedness", "Handedness & laterality")

# Register decomposed multi-response variables (dynamic names)
decomposed_prefixes <- c("stressor", "pain_month", "never_eat", 
                         "vascular", "qualification", "employment")
for (pfx in decomposed_prefixes) {
  dcols <- grep(paste0("^", pfx, "_"), names(raw), value = TRUE)
  for (dc in dcols) {
    domain <- switch(pfx,
                     stressor = "General & mental health", pain_month = "Sensory function & pain", never_eat = "Dietary intake & preferences",
                     vascular = "Clinical screening & treatment", qualification = "Demographics & SES", employment = "Demographics & SES"
    )
    vtype <- ifelse(grepl("_count$", dc), "ordinal", "binary")
    register_var(dc, dc, domain, vtype, FALSE, paste("Decomposed from", pfx))
    exposure_cols <- c(exposure_cols, dc)
  }
}

cat("  Total exposure variables registered:", nrow(var_dict), "\n")
if (length(missing_cols) > 0) {
  cat("  Variables not found in data (skipped):", length(missing_cols), "\n")
  for (mc in missing_cols) cat("    -", mc, "\n")
}

# ============================================================================
# 11. SOURCE ASSIGNMENT + PHENOTYPE-SOURCED VARIABLE REGISTRATION
# ============================================================================
#
# Cleaning policy, shared with the follow-up and disease scans:
#   - NO variable exclusion by missingness threshold at cleaning stage
#   - All variables retained; case-count filtering (n_cases >= 100) is
#     applied downstream in the regression script
#
# This step:
#   (a) Adds a 'source' column to var_dict — the exposure family, which is the
#       false-discovery-rate grouping used downstream — and checks that every
#       variable carries a recognised family name.
#   (b) Registers phenotype-sourced exposure variables (oral health,
#       treatment history, derived renal function) that originate from
#       the existing PWAS/MWAS phenotype files rather than raw baseline ExWAS.
#       Rationale: ExWAS treats phenotypic exposures (oral health,
#       treatments) as candidate exposures, NOT covariates (Pipeline
#       Manual v1.1).
# ============================================================================

cat("\n--- Step 11: Source assignment & phenotype var registration ---\n")

# ---- Step 11a: Register phenotype-sourced variables FIRST ----
# These variables live in per-group phenotype files. They are registered
# here for dictionary completeness and loaded per-group in Step 12.
# source_field = "phenotype" marks them for phenotype-file lookup.
#
# NOTE: register_var rbinds 6-column rows, so this must run BEFORE adding
# the 'source' column to var_dict (otherwise rbind fails on column mismatch).

pheno_migration <- list(
  # Oral health (8) — baseline oral health indicators
  list("oral_bleeding_gums_baseline",     "phenotype", "Oral health",
       "binary",     FALSE, "Bleeding gums at baseline"),
  list("oral_painful_gums_baseline",      "phenotype", "Oral health",
       "binary",     FALSE, "Painful gums at baseline"),
  list("oral_loose_teeth_baseline",       "phenotype", "Oral health",
       "binary",     FALSE, "Loose teeth at baseline"),
  list("oral_toothache_baseline",         "phenotype", "Oral health",
       "binary",     FALSE, "Toothache at baseline"),
  list("oral_mouth_ulcers_baseline",      "phenotype", "Oral health",
       "binary",     FALSE, "Mouth ulcers at baseline"),
  list("oral_dentures_baseline",          "phenotype", "Oral health",
       "binary",     FALSE, "Dentures at baseline"),
  list("periodontal_indicator_baseline",  "phenotype", "Oral health",
       "binary",     FALSE, "Periodontal composite indicator at baseline"),
  list("any_oral_problem_baseline",       "phenotype", "Oral health",
       "binary",     FALSE, "Any oral problem composite at baseline"),
  
  # Clinical screening & treatment (4) — baseline treatment history
  list("surg_taste_affecting_baseline",     "phenotype", "Clinical screening & treatment",
       "binary",     FALSE, "Taste-affecting surgery before baseline"),
  list("surg_non_taste_affecting_baseline", "phenotype", "Clinical screening & treatment",
       "binary",     FALSE, "Non-taste-affecting surgery before baseline (negative control)"),
  list("surg_tooth_extraction_baseline",    "phenotype", "Clinical screening & treatment",
       "binary",     FALSE, "Tooth extraction surgery before baseline"),
  list("taste_drug_user_baseline_numeric",  "phenotype", "Clinical screening & treatment",
       "binary",     FALSE, "Taste-affecting medication use at baseline"),
  
  # Blood biochemistry & haematology (1) — derived renal function
  list("egfr",                              "phenotype", "Blood biochemistry & haematology",
       "continuous", FALSE, "Estimated GFR (derived from creatinine + cystatin C)")
)

for (v in pheno_migration) {
  register_var(
    var_name     = v[[1]],
    source_field = v[[2]],
    domain       = v[[3]],
    var_type     = v[[4]],
    female_only  = v[[5]],
    description  = v[[6]]
  )
}

cat("  Registered", length(pheno_migration), "phenotype-sourced variables\n")

# ---- Step 11b: Derive the source column ----
# Applied ONCE to the complete var_dict (now including phenotype migrations).
# Each variable is tagged with its exposure family at registration, so the
# false-discovery-rate grouping is the tag itself. `source` is kept as a separate
# column because that is the name the regression engine and the feature manifest
# read, while `domain` is what the dictionary reports.
var_dict[, source := domain]

unknown_family <- setdiff(unique(var_dict$source), c(FAMILIES, "Reproductive (female-only)"))
if (length(unknown_family) > 0) {
  stop("variables tagged with an unrecognised exposure family: ",
       paste(unknown_family, collapse = " | "))
}

cat("  Total variables in dictionary:", nrow(var_dict), "\n")

# Final per-source summary
source_summary <- var_dict[, .(
  n_vars        = .N,
  n_continuous  = sum(var_type == "continuous"),
  n_ordinal     = sum(var_type == "ordinal"),
  n_binary      = sum(var_type == "binary"),
  n_categorical = sum(var_type == "categorical"),
  n_from_raw    = sum(source_field != "phenotype"),
  n_from_pheno  = sum(source_field == "phenotype")
), by = source][order(source)]

cat("\n  Final variable counts by source:\n")
print(source_summary)
cat("  TOTAL:", sum(source_summary$n_vars), "exposure variables\n")
cat("\n  NOTE: No exclusion by missingness at cleaning stage.\n")
cat("        Case-count thresholds applied in regression script.\n")

# ============================================================================
# 12. PER-SOURCE PER-GROUP OUTPUT
# ============================================================================
#
# For each (source, group) pair, produce three files:
#   baseline_exwas_<source>_<group>.csv              — exposure data
#   baseline_exwas_<source>_<group>_missingness.csv  — per-variable missingness
#   baseline_exwas_<source>_variable_dict.csv        — source-specific dict (once)
#
# Design:
#   - Variables with source_field == "phenotype" are loaded from the
#     per-group phenotype file; all others come from the raw baseline ExWAS data.
#   - eid scope: restricted to participants with non-NA taste_2w_strict.
#   - Outputs contain EXPOSURES ONLY (no outcome, no covariates).
#   - Post-processing script merges per-source per-group files into
#     unified per-group analysis-ready datasets.
# ============================================================================

cat("\n--- Step 12: Per-source per-group output ---\n")

unique_sources <- sort(unique(var_dict$source))
cat("  Sources to process:", paste(unique_sources, collapse = ", "), "\n\n")

group_reports       <- list()
source_dict_written <- character()

for (gname in names(pheno_list)) {
  gpheno <- pheno_list[[gname]]
  if (!"taste_2w_strict" %in% names(gpheno)) {
    cat("  ", gname, ": SKIPPED (taste_2w_strict not found)\n")
    next
  }
  
  g_eids_ana <- gpheno[!is.na(taste_2w_strict), eid]
  g_cases    <- gpheno[!is.na(taste_2w_strict) & taste_2w_strict == 1, eid]
  g_controls <- gpheno[!is.na(taste_2w_strict) & taste_2w_strict == 0, eid]
  
  cat("  ── GROUP", toupper(gname), "──\n")
  cat("    Phenotype N:", nrow(gpheno),
      "| Analysis N:", length(g_eids_ana),
      "| Cases:", length(g_cases),
      "| Controls:", length(g_controls), "\n")
  
  source_stats <- list()
  
  for (src in unique_sources) {
    src_dict <- var_dict[source == src]
    if (nrow(src_dict) == 0) next
    
    # Split variables by data origin
    raw_vars   <- src_dict[source_field != "phenotype", var_name]
    pheno_vars <- src_dict[source_field == "phenotype", var_name]
    
    # Only variables actually present in each data source
    raw_available   <- intersect(raw_vars,   names(raw))
    pheno_available <- intersect(pheno_vars, names(gpheno))
    
    # Build raw portion (eid scope filter)
    if (length(raw_available) > 0) {
      raw_cols <- c("eid", raw_available)
      g_raw    <- raw[eid %in% g_eids_ana, ..raw_cols]
    } else {
      g_raw    <- data.table(eid = g_eids_ana)
    }
    
    # Build phenotype portion
    if (length(pheno_available) > 0) {
      pheno_cols <- c("eid", pheno_available)
      g_phe      <- gpheno[!is.na(taste_2w_strict), ..pheno_cols]
    } else {
      g_phe      <- NULL
    }
    
    # Merge raw + phenotype portions
    if (!is.null(g_phe)) {
      g_src <- merge(g_raw, g_phe, by = "eid", all.x = TRUE)
    } else {
      g_src <- g_raw
    }
    
    # Column order: eid + src_dict variable registry order
    keep_cols <- c("eid", intersect(src_dict$var_name, names(g_src)))
    g_src     <- g_src[, ..keep_cols]
    
    # Save exposure data
    out_data <- file.path(PER_SOURCE_DIR,
                          paste0("baseline_exwas_", family_slug(src), "_", gname, ".csv"))
    fwrite(g_src, out_data)
    
    # Per-variable missingness report (with source/group columns for later aggregation)
    var_cols <- setdiff(names(g_src), "eid")
    if (length(var_cols) > 0) {
      miss_g <- rbindlist(lapply(var_cols, function(cn) {
        n_total <- nrow(g_src)
        n_na    <- sum(is.na(g_src[[cn]]))
        data.table(
          variable    = cn,
          source      = src,
          group       = gname,
          n_total     = n_total,
          n_valid     = n_total - n_na,
          n_missing   = n_na,
          pct_missing = round(100 * n_na / n_total, 1)
        )
      }))
      out_miss <- file.path(PER_SOURCE_DIR,
                            paste0("baseline_exwas_", family_slug(src), "_", gname, "_missingness.csv"))
      fwrite(miss_g, out_miss)
    }
    
    # Source-specific variable dictionary (once per source, not per group)
    if (!(src %in% source_dict_written)) {
      out_dict <- file.path(PER_SOURCE_DIR,
                            paste0("baseline_exwas_", family_slug(src), "_variable_dict.csv"))
      fwrite(src_dict, out_dict)
      source_dict_written <- c(source_dict_written, src)
    }
    
    source_stats[[src]] <- list(
      n_rows       = nrow(g_src),
      n_vars       = length(var_cols),
      n_from_raw   = length(raw_available),
      n_from_pheno = length(pheno_available),
      median_miss  = if (length(var_cols) > 0)
        round(median(miss_g$pct_missing, na.rm = TRUE), 1)
      else NA_real_
    )
  }
  
  # Print per-source summary for this group
  cat(sprintf("    %-14s | %6s | %8s | %10s | %s\n",
              "Source", "N_vars", "from_raw", "from_pheno", "median_miss%"))
  for (src in names(source_stats)) {
    ss <- source_stats[[src]]
    cat(sprintf("    %-14s | %6d | %8d | %10d | %.1f\n",
                src, ss$n_vars, ss$n_from_raw, ss$n_from_pheno, ss$median_miss))
  }
  cat("\n")
  
  group_reports[[gname]] <- list(
    n_pheno      = nrow(gpheno),
    n_analysis   = length(g_eids_ana),
    n_cases      = length(g_cases),
    n_controls   = length(g_controls),
    source_stats = source_stats
  )
}

# Save unified variable dictionary (spans all sources)
fwrite(var_dict, file.path(PER_SOURCE_DIR, "baseline_exwas_variable_dict_all.csv"))
cat("  Saved: baseline_exwas_variable_dict_all.csv (unified across sources)\n")

# ============================================================================
# 13. QC REPORT
# ============================================================================

report_file <- file.path(PER_SOURCE_DIR, "baseline_exwas_cleaning_report.txt")
sink(report_file)

cat("================================================================\n")
cat("  ExWAS baseline ExWAS — Data Cleaning QC Report\n")
cat("  Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

cat("--- INPUT ---\n")
cat("  baseline ExWAS raw CSV:", nrow(raw), "participants x", ncol(raw), "columns\n")
cat("  Groups processed:", paste(names(group_reports), collapse = ", "), "\n\n")

cat("--- CLEANING OPERATIONS ---\n")
cat("  Missing value harmonization:\n")
cat("    'Do not know' → NA\n")
cat("    'Prefer not to answer' → NA\n")
cat("    Empty string → NA\n")
cat("    'Not sure - had a hysterectomy' → NA (Field 2724)\n")
cat("  Multi-array aggregation:\n")
cat("    Blood pressure (4079, 4080): mean of two automated readings\n")
cat("    Grip strength (46, 47): max of left/right\n")
cat("  Multi-response decomposition (pipe-separated fields):\n")
cat("    Fields 6145, 6159, 6144, 6150, 6138, 6142 → binary indicators + count\n")
cat("  Derived composite variables:\n")
cat("    WHR = waist / hip\n")
cat("    PHQ-2 score = sum of Field 2050 + Field 2060 (range 0–6)\n")
cat("    SRT mean = (left + right) / 2\n")
cat("  Biologically implausible values → NA:", n_outliers_total, "values\n\n")

cat("--- FIELD SELECTION DECISIONS ---\n")
cat("  INCLUDE: Female reproductive (2714, 2724, 3581)\n")
cat("  INCLUDE: Non-accidental family death (4501)\n")
cat("  INCLUDE: PM2.5 air pollution (24006)\n")
cat("  EXCLUDE: Hair/balding pattern (2395) — insufficient theoretical link\n")
cat("  EXCLUDE: Family history (20107, 20110, 20111) — complex coding, overlaps the disease scan\n")
cat("  EXCLUDE: Ethnic background (21000) — stratification variable, not exposure\n\n")

cat("--- PHENOTYPE-SOURCED VARIABLES (migrated to baseline ExWAS exposure pool) ---\n")
cat("  Origin: existing PWAS/MWAS phenotype files\n")
cat("  Rationale: ExWAS treats oral health and treatment history as\n")
cat("             candidate exposures, not covariates (Pipeline Manual v1.1)\n")
cat("  Migrated variables (13 total):\n")
cat("    Oral health (8):   oral_bleeding_gums/painful_gums/loose_teeth/\n")
cat("                       toothache/mouth_ulcers/dentures_baseline,\n")
cat("                       periodontal_indicator_baseline, any_oral_problem_baseline\n")
cat("    Treatment (4):     surg_taste_affecting_baseline,\n")
cat("                       surg_non_taste_affecting_baseline,\n")
cat("                       surg_tooth_extraction_baseline,\n")
cat("                       taste_drug_user_baseline_numeric\n")
cat("    Biochemistry (1):  egfr (derived from creatinine + cystatin C)\n\n")

cat("--- VARIABLE COUNTS BY SOURCE ---\n")
print(source_summary)
cat("  TOTAL:", sum(source_summary$n_vars), "exposure variables\n\n")

cat("--- NO MISSINGNESS-BASED EXCLUSION ---\n")
cat("  All variables retained at cleaning stage. Case-count thresholds\n")
cat("  (n_cases >= 100 per exposure) are applied downstream in the\n")
cat("  regression script, as in the follow-up and disease scans.\n\n")

cat("--- PER-GROUP SAMPLE AND PER-SOURCE BREAKDOWN ---\n")
for (gname in names(group_reports)) {
  gr <- group_reports[[gname]]
  cat("\n  ", toupper(gname), ":\n")
  cat("    Phenotype file:", gr$n_pheno, "participants\n")
  cat("    Analysis sample:", gr$n_analysis,
      "(cases:", gr$n_cases, "| controls:", gr$n_controls, ")\n")
  cat("    Per-source breakdown:\n")
  cat(sprintf("      %-14s  %-6s  %-10s  %-12s  %s\n",
              "Source", "N_vars", "from_raw", "from_pheno", "median_miss%"))
  for (src in names(gr$source_stats)) {
    ss <- gr$source_stats[[src]]
    cat(sprintf("      %-14s  %-6d  %-10d  %-12d  %.1f\n",
                src, ss$n_vars, ss$n_from_raw, ss$n_from_pheno, ss$median_miss))
  }
}
cat("\n")

cat("--- FEMALE-ONLY VARIABLES ---\n")
female_vars <- var_dict[female_only == TRUE]
cat("  Count:", nrow(female_vars), "\n")
if (nrow(female_vars) > 0) print(female_vars[, .(var_name, domain, source)])
cat("  Analysis: sex-stratified regression (handled at regression stage)\n\n")

cat("--- OUTPUT FILES ---\n")
cat("  Per-source per-group data:       baseline_exwas_<source>_<group>.csv\n")
cat("  Per-source per-group missingness: baseline_exwas_<source>_<group>_missingness.csv\n")
cat("  Per-source variable dictionary:   baseline_exwas_<source>_variable_dict.csv\n")
cat("  Unified variable dictionary:      baseline_exwas_variable_dict_all.csv\n")
cat("  QC report:                        baseline_exwas_cleaning_report.txt\n\n")

cat("--- NEXT: assembly (sections 14-19 below) ---\n")
cat("    1. Audit and aggregate the per-family per-group files\n")
cat("    2. Produce one analysis-ready matrix per cohort\n")
cat("    3. Compute per-variable case counts for regression threshold planning\n")
cat("    4. Write the summary report\n")

cat("\n================================================================\n")
cat("  END OF REPORT\n")
cat("================================================================\n")

sink()
cat("\n  Saved:", report_file, "\n")

cat("\n================================================================\n")
cat("  baseline ExWAS cleaning pipeline completed successfully\n")
cat("  Total exposure variables:", nrow(var_dict), "\n")
cat("  Sources output:", length(unique_sources), "\n")
cat("  Groups output:", paste(names(group_reports), collapse = ", "), "\n")
cat("================================================================\n")

# ============================================================================
# UTILITY: count non-NA per row (robust to any column count)
# ============================================================================

n_nonNA_per_row <- function(dt) {
  if (is.null(dt) || ncol(dt) == 0) return(integer(0))
  Reduce("+", lapply(dt, function(x) as.integer(!is.na(x))))
}

# ============================================================================
# ASSEMBLY CONFIGURATION
# ============================================================================

GROUPS <- c("group1", "group2", "group3")
GROUP_LABELS <- c(
  group1 = "Pure White British",
  group2 = "Other White",
  group3 = "Non-White"
)

PHENO_VAR <- "taste_2w_strict"
META_COLS <- c("eid")

# ============================================================================
# STEP 0: LOAD PHENOTYPE
# ============================================================================

cat("================================================================\n")
cat("  baseline ExWAS Post-Processing\n")
cat("================================================================\n")
cat("Data input  (per-source CSVs):", PER_SOURCE_DIR, "\n")
cat("Pheno input (phenotype CSVs): ", PHENO_DIR, "\n")
cat("Output      (merged datasets):", DERIVE_DIR, "\n")

cat("\nLoading phenotype files from:", PHENO_DIR, "\n")
pheno_list <- list()
for (g in GROUPS) {
  pf <- PHENO_FILES[[g]]
  if (file.exists(pf)) {
    pheno_list[[g]] <- fread(pf, select = c("eid", PHENO_VAR))
    cat("  Loaded ", g, ": ",
        format(nrow(pheno_list[[g]]), big.mark = ","), " rows\n", sep = "")
  } else {
    cat("  ERROR: ", pf, " not found\n", sep = "")
  }
}

# Cohort sizes & case counts
cohort_stats <- list()
for (g in GROUPS) {
  if (is.null(pheno_list[[g]])) next
  p <- pheno_list[[g]]
  if (!PHENO_VAR %in% names(p)) next
  p_with_outcome <- p[!is.na(get(PHENO_VAR))]
  cohort_stats[[g]] <- list(
    n_total    = nrow(p_with_outcome),
    n_cases    = sum(p_with_outcome[[PHENO_VAR]] == 1),
    n_controls = sum(p_with_outcome[[PHENO_VAR]] == 0)
  )
  cohort_stats[[g]]$case_rate <- cohort_stats[[g]]$n_cases / cohort_stats[[g]]$n_total
}

cat("\nCohort sizes (with phenotype):\n")
for (g in GROUPS) {
  cs <- cohort_stats[[g]]
  if (is.null(cs)) next
  cat(sprintf("  %s: total=%s, cases=%s (%.2f%%), controls=%s\n",
              g,
              format(cs$n_total, big.mark = ","),
              format(cs$n_cases, big.mark = ","),
              100 * cs$case_rate,
              format(cs$n_controls, big.mark = ",")))
}

# ============================================================================
# STEP 1: FILE AUDIT
# ============================================================================

cat("\n--- Step 1: File audit ---\n")

audit <- data.table()
for (src in FAMILIES) {
  for (g in GROUPS) {
    df <- file.path(PER_SOURCE_DIR, paste0("baseline_exwas_", family_slug(src), "_", g, ".csv"))
    mf <- file.path(PER_SOURCE_DIR, paste0("baseline_exwas_", family_slug(src), "_", g, "_missingness.csv"))
    audit <- rbind(audit, data.table(
      source = src, group = g,
      data_file = df, miss_file = mf,
      data_exists = file.exists(df),
      miss_exists = file.exists(mf)
    ))
  }
}

n_data_missing <- sum(!audit$data_exists)
if (n_data_missing > 0) {
  cat("  WARNING:", n_data_missing, "data files missing:\n")
  for (i in which(!audit$data_exists)) {
    cat("    -", basename(audit$data_file[i]), "\n")
  }
}
cat("  Files found:", sum(audit$data_exists), "/", nrow(audit), "\n")

# ============================================================================
# STEP 2: AGGREGATE PER-SOURCE MISSINGNESS
# ============================================================================

cat("\n--- Step 2: Aggregate per-source missingness ---\n")

miss_list <- list()
for (i in seq_len(nrow(audit))) {
  if (!audit$miss_exists[i]) next
  m <- fread(audit$miss_file[i])
  # Ensure source/group columns exist (the cleaning script writes them)
  if (!"source" %in% names(m)) m[, source := audit$source[i]]
  if (!"group"  %in% names(m)) m[, group  := audit$group[i]]
  miss_list[[paste(audit$source[i], audit$group[i], sep = "_")]] <- m
}

if (length(miss_list) > 0) {
  combined_miss <- rbindlist(miss_list, fill = TRUE)
  fwrite(combined_miss,
         file.path(DERIVE_DIR, "baseline_exwas_combined_missingness.csv"))
  cat("  Aggregated missingness rows:", nrow(combined_miss), "\n")
  cat("  Saved: baseline_exwas_combined_missingness.csv\n")
} else {
  combined_miss <- data.table()
  cat("  WARNING: no missingness files found\n")
}

# ============================================================================
# STEP 3: LOAD PER-SECTION DATA + COMPUTE CASE COUNTS + SOURCE OVERLAP
# ============================================================================

cat("\n--- Step 3: Load per-source data, compute case counts ---\n")

per_section_data <- list()      # key = "source_group" → data.table
case_counts      <- data.table()
source_overlap_raw <- list()    # per group: per-eid source count

for (g in GROUPS) {
  source_overlap_raw[[g]] <- data.table(eid = integer(0))
}

for (i in seq_len(nrow(audit))) {
  if (!audit$data_exists[i]) next
  src <- audit$source[i]
  g   <- audit$group[i]
  key <- paste(src, g, sep = "_")
  
  dat <- fread(audit$data_file[i])
  per_section_data[[key]] <- dat
  
  var_cols <- setdiff(names(dat), META_COLS)
  if (length(var_cols) == 0) next
  
  # Per-variable case counts (against phenotype)
  if (!is.null(pheno_list[[g]])) {
    pheno <- pheno_list[[g]][!is.na(get(PHENO_VAR))]
    m <- merge(dat, pheno, by = "eid")
    
    for (v in var_cols) {
      complete_idx <- !is.na(m[[v]])
      n_complete   <- sum(complete_idx)
      if (n_complete == 0) {
        n_cases_c  <- 0L
        case_rate  <- NA_real_
      } else {
        n_cases_c <- sum(m[complete_idx, get(PHENO_VAR)] == 1)
        case_rate <- 100 * n_cases_c / n_complete
      }
      case_counts <- rbind(case_counts, data.table(
        variable           = v,
        source             = src,
        group              = g,
        n_complete         = n_complete,
        n_cases_complete   = n_cases_c,
        case_rate_complete = round(case_rate, 3)
      ))
    }
  }
  
  # Source overlap tracking: eids with any non-missing value in this source
  any_data <- n_nonNA_per_row(dat[, ..var_cols]) > 0
  eids_with_data <- dat$eid[any_data]
  if (length(eids_with_data) > 0) {
    contrib <- data.table(eid = eids_with_data, src = src)
    source_overlap_raw[[g]] <- rbind(source_overlap_raw[[g]], contrib, fill = TRUE)
  }
}

if (nrow(case_counts) > 0) {
  fwrite(case_counts,
         file.path(DERIVE_DIR, "baseline_exwas_variable_case_counts.csv"))
  cat("  Per-variable case counts computed:", nrow(case_counts), "rows\n")
  cat("  Saved: baseline_exwas_variable_case_counts.csv\n")
}

# Source overlap summary per group
source_overlap <- data.table()
for (g in GROUPS) {
  raw <- source_overlap_raw[[g]]
  if (nrow(raw) == 0) next
  counts <- raw[, .(n_sources = uniqueN(src)), by = eid]
  summary_counts <- counts[, .(N = .N), by = n_sources][order(n_sources)]
  summary_counts[, group := g]
  source_overlap <- rbind(source_overlap, summary_counts, fill = TRUE)
}

if (nrow(source_overlap) > 0) {
  fwrite(source_overlap,
         file.path(DERIVE_DIR, "baseline_exwas_source_overlap.csv"))
  cat("  Source overlap computed\n")
  cat("  Saved: baseline_exwas_source_overlap.csv\n")
}

# ============================================================================
# STEP 4: MERGE PER-GROUP ANALYSIS-READY DATASETS
# ============================================================================

cat("\n--- Step 4: Merge per-group analysis-ready datasets ---\n")

merge_stats <- list()

for (g in GROUPS) {
  # Collect all per-source data for this group
  group_keys <- grep(paste0("_", g, "$"), names(per_section_data), value = TRUE)
  if (length(group_keys) == 0) {
    cat("  ", g, ": no data found, skipping\n")
    next
  }
  
  # Start with the first source's eid scope
  merged <- per_section_data[[group_keys[1]]]
  
  for (k in group_keys[-1]) {
    nxt <- per_section_data[[k]]
    # Full outer join — preserve all eids from all sources (should be same
    # scope since all sources share the same group-filtered eids, but use
    # all = TRUE defensively)
    merged <- merge(merged, nxt, by = "eid", all = TRUE)
  }
  
  var_cols     <- setdiff(names(merged), META_COLS)
  n_total_eids <- nrow(merged)
  n_any_data   <- sum(n_nonNA_per_row(merged[, ..var_cols]) > 0)
  
  out_file <- file.path(DERIVE_DIR,
                        paste0("exwas_baseline_", g, ".csv"))
  fwrite(merged, out_file)
  
  merge_stats[[g]] <- list(
    n_total_eids = n_total_eids,
    n_any_b_data = n_any_data,
    n_total_vars = length(var_cols)
  )
  
  cat("  ", g, ": N=", format(n_total_eids, big.mark = ","),
      " vars=", length(var_cols),
      " any_data=", format(n_any_data, big.mark = ","), "\n", sep = "")
  cat("    Saved:", basename(out_file), "\n")
}

# ============================================================================
# STEP 5: COPY UNIFIED VARIABLE DICTIONARY
# ============================================================================

cat("\n--- Step 5: Copy unified variable dictionary ---\n")

unified_dict_src <- file.path(PER_SOURCE_DIR, "baseline_exwas_variable_dict_all.csv")
unified_dict_dst <- file.path(DERIVE_DIR, "baseline_exwas_variable_dictionary.csv")

if (file.exists(unified_dict_src)) {
  var_dict_all <- fread(unified_dict_src)
  fwrite(var_dict_all, unified_dict_dst)
  cat("  Copied:", basename(unified_dict_dst), "\n")
  cat("  Total variables:", nrow(var_dict_all), "\n")
} else {
  cat("  WARNING: unified dictionary not found at", unified_dict_src, "\n")
  var_dict_all <- NULL
}

# ============================================================================
# STEP 6: HUMAN-READABLE SUMMARY REPORT
# ============================================================================

cat("\n--- Step 6: Writing summary report ---\n")

report_file <- file.path(DERIVE_DIR, "baseline_exwas_summary_report.txt")
con <- file(report_file, "w")
wln <- function(...) writeLines(paste0(...), con = con)

wln("================================================================")
wln("  ExWAS baseline ExWAS — Post-Processing Summary Report")
wln("  Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
wln("================================================================")
wln("")

# ----- Section 1: Cohort sizes -----
wln("[1] COHORT SIZES (with valid taste_2w_strict)")
wln("--------------------------------------------------------")
wln(sprintf("%-10s %-22s %12s %12s %10s",
            "Group", "Label", "Total", "Cases", "Case%"))
for (g in GROUPS) {
  cs <- cohort_stats[[g]]
  if (is.null(cs)) next
  wln(sprintf("%-10s %-22s %12s %12s %9.2f%%",
              g, GROUP_LABELS[g],
              format(cs$n_total, big.mark = ","),
              format(cs$n_cases, big.mark = ","),
              100 * cs$case_rate))
}
wln("")

# ----- Section 2: Per-source × per-group -----
wln("[2] PER-SOURCE x PER-GROUP SUMMARY")
wln("--------------------------------------------------------")

for (src in FAMILIES) {
  src_keys <- grep(paste0("^", src, "_"), names(per_section_data), value = TRUE)
  if (length(src_keys) == 0) next
  
  wln("")
  wln("  ", src)
  wln("  ", strrep("-", 56))
  wln(sprintf("  %-10s %12s %12s %12s %8s",
              "Group", "Respondents", "Resp_w_pheno", "Cases_resp", "N_vars"))
  
  for (g in GROUPS) {
    key <- paste(src, g, sep = "_")
    dat <- per_section_data[[key]]
    if (is.null(dat)) {
      wln(sprintf("  %-10s %12s %12s %12s %8s",
                  g, "-", "-", "-", "-"))
      next
    }
    
    var_cols <- setdiff(names(dat), META_COLS)
    n_resp   <- sum(n_nonNA_per_row(dat[, ..var_cols]) > 0)
    
    if (!is.null(pheno_list[[g]])) {
      pheno <- pheno_list[[g]][!is.na(get(PHENO_VAR)),
                               .(eid, outcome = get(PHENO_VAR))]
      has_idx   <- n_nonNA_per_row(dat[, ..var_cols]) > 0
      resp_eids <- dat[has_idx, .(eid)]
      m         <- merge(resp_eids, pheno, by = "eid")
      n_resp_pheno <- nrow(m)
      n_cases_resp <- sum(m$outcome == 1)
    } else {
      n_resp_pheno <- NA_integer_
      n_cases_resp <- NA_integer_
    }
    
    wln(sprintf("  %-10s %12s %12s %12s %8d",
                g,
                format(n_resp, big.mark = ","),
                format(n_resp_pheno, big.mark = ","),
                format(n_cases_resp, big.mark = ","),
                length(var_cols)))
  }
}
wln("")

# ----- Section 3: Source overlap -----
wln("")
wln("[3] SOURCE OVERLAP (number of A sources with any non-missing data)")
wln("--------------------------------------------------------")

if (nrow(source_overlap) > 0) {
  for (g in GROUPS) {
    sub <- source_overlap[group == g]
    if (nrow(sub) == 0) next
    wln("")
    wln("  ", g, " (", GROUP_LABELS[g], "):")
    
    denom <- if (!is.null(cohort_stats[[g]])) cohort_stats[[g]]$n_total else sum(sub$N)
    
    for (i in seq_len(nrow(sub))) {
      ns  <- sub$n_sources[i]
      n   <- sub$N[i]
      pct <- 100 * n / denom
      wln(sprintf("    %2d source(s): %12s  (%5.2f%% of cohort)",
                  ns, format(n, big.mark = ","), pct))
    }
  }
}
wln("")

# ----- Section 4: Final merged datasets -----
wln("")
wln("[4] FINAL MERGED ANALYSIS-READY DATASETS")
wln("--------------------------------------------------------")
wln(sprintf("  %-10s %15s %15s %12s",
            "Group", "Total_eids", "With_any_A", "N_variables"))
for (g in GROUPS) {
  ms <- merge_stats[[g]]
  if (is.null(ms)) next
  wln(sprintf("  %-10s %15s %15s %12d",
              g,
              format(ms$n_total_eids, big.mark = ","),
              format(ms$n_any_b_data, big.mark = ","),
              ms$n_total_vars))
}
wln("")

# ----- Section 5: Case-count distribution (for regression threshold) -----
wln("")
wln("[5] VARIABLE CASE-COUNT DISTRIBUTION (group1)")
wln("--------------------------------------------------------")
wln("    How many variables have >=N cases with non-missing exposure?")
wln("    Useful for choosing the regression filter threshold.")
wln("")

if (nrow(case_counts) > 0) {
  cc_g1 <- case_counts[group == "group1"]
  thresholds <- c(20, 50, 100, 200, 500, 1000)
  
  wln(sprintf("  %-14s %12s %12s",
              "Threshold", "N_vars", "% of total"))
  total_vars_g1 <- nrow(cc_g1)
  for (th in thresholds) {
    n_pass <- sum(cc_g1$n_cases_complete >= th)
    wln(sprintf("  cases>=%-7d %12d %11.1f%%",
                th, n_pass, 100 * n_pass / total_vars_g1))
  }
  
  wln("")
  wln("  Per-source breakdown at cases>=100 threshold (group1):")
  for (src in FAMILIES) {
    sub <- cc_g1[source == src]
    if (nrow(sub) == 0) next
    n_pass <- sum(sub$n_cases_complete >= 100)
    wln(sprintf("    %-14s: %3d / %3d variables pass",
                src, n_pass, nrow(sub)))
  }
}
wln("")

# ----- Section 6: Top & bottom variables by case count -----
wln("")
wln("[6] TOP 10 VARIABLES BY CASE COUNT (group1)")
wln("--------------------------------------------------------")
if (nrow(case_counts) > 0) {
  top10 <- case_counts[group == "group1"][order(-n_cases_complete)][1:10]
  wln(sprintf("  %-32s %-14s %12s %12s %10s",
              "Variable", "Source", "N_complete", "N_cases", "Case%"))
  for (i in seq_len(nrow(top10))) {
    wln(sprintf("  %-32s %-14s %12s %12s %9.3f%%",
                substr(top10$variable[i], 1, 32),
                top10$source[i],
                format(top10$n_complete[i], big.mark = ","),
                format(top10$n_cases_complete[i], big.mark = ","),
                top10$case_rate_complete[i]))
  }
}
wln("")

wln("")
wln("[7] BOTTOM 10 VARIABLES BY CASE COUNT (group1, n_complete>0)")
wln("--------------------------------------------------------")
wln("    Candidates for downstream filtering at cases>=100 threshold")
if (nrow(case_counts) > 0) {
  bot10 <- case_counts[group == "group1" & n_complete > 0][
    order(n_cases_complete)][1:10]
  wln(sprintf("  %-32s %-14s %12s %12s %10s",
              "Variable", "Source", "N_complete", "N_cases", "Case%"))
  for (i in seq_len(nrow(bot10))) {
    wln(sprintf("  %-32s %-14s %12s %12s %9.3f%%",
                substr(bot10$variable[i], 1, 32),
                bot10$source[i],
                format(bot10$n_complete[i], big.mark = ","),
                format(bot10$n_cases_complete[i], big.mark = ","),
                bot10$case_rate_complete[i]))
  }
}
wln("")

# ----- Section 8: Files generated -----
wln("")
wln("[8] FILES GENERATED")
wln("--------------------------------------------------------")
wln("  Merged datasets:")
for (g in GROUPS) {
  wln("    exwas_baseline_", g, ".csv")
}
wln("  Variable dictionary:")
wln("    baseline_exwas_variable_dictionary.csv")
wln("  Diagnostics:")
wln("    baseline_exwas_combined_missingness.csv")
wln("    baseline_exwas_variable_case_counts.csv")
wln("    baseline_exwas_source_overlap.csv")
wln("    baseline_exwas_summary_report.txt  (this file)")
wln("")
wln("================================================================")
wln("  END OF REPORT")
wln("================================================================")

close(con)

cat("  Saved:", basename(report_file), "\n")

# ============================================================================
# DONE
# ============================================================================

cat("\n================================================================\n")
cat("  baseline ExWAS post-processing complete.\n")
cat("================================================================\n")
cat("\nAll outputs written to:", DERIVE_DIR, "\n")
cat("\nNext step: Review baseline_exwas_summary_report.txt for:\n")
cat("  - Per-source case counts (for regression n_cases >= 100 threshold)\n")
cat("  - Cohort sizes and source overlap patterns\n")
cat("  - Methods-ready numbers for manuscript drafting\n")
