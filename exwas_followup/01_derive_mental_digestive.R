#!/usr/bin/env Rscript
#> in   input/raw/followup_exwas_merged.csv
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

# 0. CONFIGURATION — paths resolved under $PROJECT_DIR (override as needed)

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
CODE_DIR    <- Sys.getenv("CODE_DIR", unset = ".")
MERGED_CSV  <- Sys.getenv("MERGED_CSV",
                          unset = file.path(PROJECT_DIR, "input", "raw",
                                            "followup_exwas_merged.csv"))
PHENO_DIR   <- Sys.getenv("PHENO_DIR", unset = file.path(PROJECT_DIR, "input", "analysis_ready"))
OUTPUT_DIR  <- Sys.getenv("DERIVE_PER_SECTION_DIR",
                          unset = file.path(PROJECT_DIR, "output",
                                            "followup_exwas", "derive", "per_section"))

MISS_THRESHOLD_PCT <- 50   # Flag variables with >50% missing in Group 1

PHENO_FILES <- list(
  group1 = file.path(PHENO_DIR, "phenotype_group1.csv"),
  group2 = file.path(PHENO_DIR, "phenotype_group2.csv"),
  group3 = file.path(PHENO_DIR, "phenotype_group3.csv")
)

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("================================================================\n")
cat("  Follow-up ExWAS — data load, digestive health, mental health\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

FIELD_LIST <- file.path(CODE_DIR, "field_list.txt")
REQUESTED_COLS   <- setdiff(trimws(readLines(FIELD_LIST)), c("", "eid"))
REQUESTED_FIELDS <- suppressWarnings(as.integer(sub("^p(\\d+).*$", "\\1", REQUESTED_COLS)))
REQUESTED_FIELDS <- sort(unique(REQUESTED_FIELDS[!is.na(REQUESTED_FIELDS)]))
cat("--- Step 0: Export specification ---\n")
cat("  field_list.txt:", length(REQUESTED_COLS), "columns,",
    length(REQUESTED_FIELDS), "field ids\n")

FIELD_REGISTER <- new.env(parent = emptyenv())

register_fields <- function(family, ids, lookup_only = FALSE) {
  ids <- sort(unique(as.integer(ids)))
  absent <- setdiff(ids, REQUESTED_FIELDS)
  if (length(absent) > 0) {
    if (!lookup_only)
      stop(family, ": declared fields are absent from field_list.txt: ",
           paste(absent, collapse = ", "))
    cat("  ", family, ": lookup carries ", length(absent),
        " field(s) outside the export, unused (",
        paste(absent, collapse = ", "), ")\n", sep = "")
  }
  assign(family, intersect(ids, REQUESTED_FIELDS), envir = FIELD_REGISTER)
  cat("  ", family, ": ", length(intersect(ids, REQUESTED_FIELDS)),
      " field ids registered against field_list.txt\n", sep = "")
}

#' Every field in field_list.txt must be claimed by some derivation step.
check_field_list_complete <- function() {
  claimed <- sort(unique(unlist(mget(ls(FIELD_REGISTER), envir = FIELD_REGISTER),
                                use.names = FALSE)))
  unclaimed <- setdiff(REQUESTED_FIELDS, claimed)
  if (length(unclaimed) > 0)
    stop("field_list.txt lists ", length(unclaimed),
         " field id(s) no derivation step reads: ",
         paste(unclaimed, collapse = ", "))
  cat("\n  field_list.txt agrees with the derivation steps: ",
      length(claimed), " field ids across ", length(ls(FIELD_REGISTER)),
      " families\n", sep = "")
}

# 1. DATA LOADING

cat("--- Step 1: Loading data ---\n")

raw <- fread(MERGED_CSV, na.strings = c("NA"))
cat("  follow-up ExWAS raw:", nrow(raw), "rows x", ncol(raw), "columns\n")

# Load phenotype files for sample filtering
pheno_list <- list()
for (gname in names(PHENO_FILES)) {
  fp <- PHENO_FILES[[gname]]
  if (file.exists(fp)) {
    pheno_list[[gname]] <- fread(fp, na.strings = c("", "NA"))
    cat("  Phenotype", gname, ":", nrow(pheno_list[[gname]]), "rows\n")
  } else {
    cat("  WARNING:", gname, "file not found — skipping\n")
  }
}

# 2. UNIVERSAL EMPTY STRING → NA CONVERSION

cat("\n--- Step 2: Universal empty string → NA conversion ---\n")

char_cols <- names(raw)[sapply(raw, is.character)]
n_converted_total <- 0L

for (cn in char_cols) {
  n_empty <- sum(raw[[cn]] == "", na.rm = TRUE)
  if (n_empty > 0) {
    raw[get(cn) == "", (cn) := NA_character_]
    n_converted_total <- n_converted_total + n_empty
  }
}

cat("  Character columns processed:", length(char_cols), "\n")
cat("  Total empty strings → NA:", format(n_converted_total, big.mark = ","), "\n")

# 3. UNIVERSAL MISSING LABEL CLEANING

cat("\n--- Step 3: Cleaning special missing labels ---\n")

MISSING_LABELS <- c("Do not know", "Prefer not to answer",
                     "do not wish to answer")

n_special_total <- 0L
for (cn in char_cols) {
  n_special <- sum(raw[[cn]] %in% MISSING_LABELS, na.rm = TRUE)
  if (n_special > 0) {
    raw[get(cn) %in% MISSING_LABELS, (cn) := NA_character_]
    n_special_total <- n_special_total + n_special
  }
}

cat("  Special missing labels → NA:", format(n_special_total, big.mark = ","), "\n")

# 4. COLUMN NAME PARSER

parse_colname <- function(cn) {
  m <- str_match(cn, "^p(\\d+)(?:_i(\\d+))?(?:_a(\\d+))?$")
  if (is.na(m[1, 1])) return(NULL)
  data.table(
    colname  = cn,
    field_id = as.integer(m[1, 2]),
    instance = ifelse(is.na(m[1, 3]), NA_integer_, as.integer(m[1, 3])),
    array_ix = ifelse(is.na(m[1, 4]), NA_integer_, as.integer(m[1, 4]))
  )
}

col_map <- rbindlist(lapply(setdiff(names(raw), "eid"), parse_colname))
cat("\n  Parsed columns:", nrow(col_map), "from", ncol(raw) - 1, "total\n")

# 5. DOMAIN ASSIGNMENT

cat("\n--- Step 5: Domain assignment ---\n")

col_map[, domain := fcase(
  field_id >= 20600 & field_id <= 20748,                            "Food preferences",
  field_id >= 26000 & field_id <= 26061,                            "Diet (24-hour recall)",
  field_id >= 26062 & field_id <= 26155,                            "Diet (24-hour recall)",
  field_id %in% c(20084, 104670),                                   "Diet (24-hour recall)",
  field_id >= 21024 & field_id <= 21068,                            "Digestive health",
  (field_id >= 20403 & field_id <= 20416) |                         # AUDIT/GAD-7
    (field_id >= 20421 & field_id <= 20446) |                       # Mania
    (field_id >= 20492 & field_id <= 20502) |                       # Unusual exp
    (field_id >= 20505 & field_id <= 20520) |                       # PHQ-9 + dep
    field_id == 20548,                                              "Mental health",
  field_id >= 120008 & field_id <= 120103,                          "Experience of pain",
  field_id %in% c(20132, 20156, 20157, 20159, 20191, 20760),       "Cognitive function",
  field_id >= 22606 & field_id <= 22651,                            "Work environment",
  default = "UNCLASSIFIED"
)]

domain_summary <- col_map[, .(n_cols = .N, n_fields = uniqueN(field_id)), by = domain]
setorder(domain_summary, domain)
cat("\n"); print(domain_summary); cat("\n")

if (any(col_map$domain == "UNCLASSIFIED")) {
  cat("  WARNING: Unclassified columns:\n")
  print(col_map[domain == "UNCLASSIFIED"])
}

#' Recode Yes/No character to binary integer
recode_yesno <- function(x) {
  fcase(x == "Yes", 1L, x == "No", 0L, default = NA_integer_)
}

#' Recode ordinal text to integer based on named vector map
recode_ordinal <- function(x, map) {
  out <- map[x]
  out[is.na(names(out))] <- NA_integer_
  as.integer(out)
}

safe_as_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

#' Compute missingness stats for a set of columns in a data.table
missingness_report <- function(dt, cols, group_eid = NULL) {
  if (!is.null(group_eid)) {
    dt <- dt[eid %in% group_eid]
  }
  rbindlist(lapply(cols, function(cn) {
    if (!cn %in% names(dt)) return(NULL)
    n_total <- nrow(dt)
    n_na <- sum(is.na(dt[[cn]]))
    data.table(
      variable = cn,
      n_total  = n_total,
      n_valid  = n_total - n_na,
      n_missing = n_na,
      pct_missing = round(100 * n_na / n_total, 1)
    )
  }))
}

#' Build variable dictionary entry
dict_entry <- function(var_name, field_id, domain, source_cat,
                       var_type, coding, taste_relevance) {
  data.table(
    var_name    = var_name,
    field_id    = field_id,
    domain      = domain,
    source_cat  = source_cat,
    var_type    = var_type,
    coding      = coding,
    taste_relevance = taste_relevance
  )
}

cat("\n================================================================\n")
cat("  Digestive health (Cat 153)\n")
cat("================================================================\n\n")

digestive_registry <- list(
  
  # --- Digestive health: IBS & Core GI Symptoms ---
  
  ibs_ever = list(
    field = 21024, col = "p21024",
    type  = "binary", 
    recode = "Yes=1, No=0",
    relevance = "Gut-brain axis sensory pathway core marker"
  ),
  
  abdom_discomfort_freq = list(
    field = 21025, col = "p21025",
    type  = "ordinal",
    recode = "Never=0, Less than one day a month=1, One day a month=2, Two to three days a month=3, One day a week=4, More than one day a week=5, Every day=6",
    relevance = "Current GI symptom burden"
  ),
  
  abdom_pain_6mo = list(
    field = 21027, col = "p21027",
    type  = "binary",
    recode = "Yes=1, No=0",
    relevance = "Chronic vs acute GI distinction"
  ),
  
  abdom_pain_severity = list(
    field = 21036, col = "p21036",
    type  = "continuous",
    recode = "0-10 NRS (numeric, kept as-is after as.numeric)",
    relevance = "Dose-response for GI pain severity"
  ),
  
  abdom_distension = list(
    field = 21038, col = "p21038",
    type  = "binary",
    recode = "Yes=1, No=0",
    relevance = "Functional GI indicator"
  ),
  
  bowel_satisfaction = list(
    field = 21040, col = "p21040",
    type  = "continuous",
    recode = "0-10 NRS (0=very satisfied, 10=very dissatisfied)",
    relevance = "Overall GI function self-assessment"
  ),
  
  gi_life_interference = list(
    field = 21041, col = "p21041",
    type  = "continuous",
    recode = "0-10 NRS (0=not at all, 10=completely)",
    relevance = "Functional impact of GI symptoms"
  ),
  
  bowel_movements_day = list(
    field = 21044, col = "p21044",
    type  = "continuous",
    recode = "Integer count per day (kept as-is after as.numeric)",
    relevance = "Diarrhea-type vs constipation-type indicator"
  ),
  
  # --- Digestive health: Visceral Sensitivity & History ---
  
  sensitive_stomach = list(
    field = 21064, col = "p21064",
    type  = "binary",
    recode = "Yes=1, No=0",
    relevance = "Visceral hypersensitivity self-report; may share sensory threshold mechanism with taste sensitivity"
  ),
  
  ibs_family_history = list(
    field = 21065, col = "p21065",
    type  = "binary",
    recode = "Yes=1, No=0",
    relevance = "Genetic predisposition proxy for GI-taste pathway"
  ),
  
  coeliac_gluten = list(
    field = 21068, col = "p21068",
    type  = "binary",
    recode = "Yes=1, No=0",
    relevance = "DIRECT: malabsorption → Zn/B12/Fe deficiency → taste dysfunction"
  ),
  
  # --- Digestive health, block c: PHQ-12 SS Somatic Symptom Items ---
  
  phq12_headache = list(
    field = 21051, col = "p21051",
    type  = "ordinal",
    recode = "Not bothered at all=0, Bothered a little=1, Bothered a lot=2",
    relevance = "PHQ-12 SS somatic: headaches"
  ),

  phq12_dizziness = list(
    field = 21053, col = "p21053",
    type  = "ordinal",
    recode = "Not bothered at all=0, Bothered a little=1, Bothered a lot=2",
    relevance = "PHQ-12 SS somatic: dizziness"
  ),

  phq12_nausea = list(
    field = 21059, col = "p21059",
    type  = "ordinal",
    recode = "Not bothered at all=0, Bothered a little=1, Bothered a lot=2",
    relevance = "PHQ-12 SS somatic: nausea"
  ),

  phq12_fatigue = list(
    field = 21060, col = "p21060",
    type  = "ordinal",
    recode = "Not bothered at all=0, Bothered a little=1, Bothered a lot=2",
    relevance = "PHQ-12 SS somatic: feeling tired all the time"
  )
  
  # DERIVED: phq12_somatic_total (sum of 4 items, range 0-8) — see below
)

cat("  Registered digestive-health fields:", length(digestive_registry), "(+ 1 derived)\n")
register_fields("Digestive health",
                vapply(digestive_registry, function(x) x$field, numeric(1)))

# Verify all fields present
digestive_cols <- sapply(digestive_registry, function(x) x$col)
digestive_missing <- digestive_cols[!digestive_cols %in% names(raw)]
if (length(digestive_missing) > 0) {
  cat("  WARNING: Missing columns:", paste(digestive_missing, collapse = ", "), "\n")
} else {
  cat("  All declared fields found in data [OK]\n")
}

cat("\n--- Digestive health, step 2: Recoding ---\n")

# Subset to working copy (eid + digestive-health columns only)
digestive <- raw[, c("eid", digestive_cols), with = FALSE]

# --- Binary Yes/No fields ---
binary_fields <- c("p21024", "p21027", "p21038", "p21064", "p21065", "p21068")

for (cn in binary_fields) {
  var_name <- names(digestive_registry)[digestive_cols == cn]
  digestive[, (var_name) := recode_yesno(get(cn))]
  n_yes <- sum(digestive[[var_name]] == 1, na.rm = TRUE)
  n_no  <- sum(digestive[[var_name]] == 0, na.rm = TRUE)
  cat("  ", var_name, ": Yes=", n_yes, " No=", n_no, 
      " NA=", sum(is.na(digestive[[var_name]])), "\n")
}

# --- Ordinal: Abdominal discomfort frequency (p21025) ---
abdom_freq_map <- c(
  "Never"                      = 0L,
  "Less than one day a month"  = 1L,
  "One day a month"            = 2L,
  "Two to three days a month"  = 3L,
  "One day a week"             = 4L,
  "More than one day a week"   = 5L,
  "Every day"                  = 6L
)
digestive[, abdom_discomfort_freq := abdom_freq_map[p21025]]
cat("  abdom_discomfort_freq: mapped", sum(!is.na(digestive$abdom_discomfort_freq)),
    "values (0-6 ordinal)\n")

# --- Numeric Rating Scales (0-10): p21036, p21040, p21041 ---

nrs_fields <- list(
  abdom_pain_severity  = "p21036",
  bowel_satisfaction    = "p21040",
  gi_life_interference = "p21041"
)

for (var_name in names(nrs_fields)) {
  cn <- nrs_fields[[var_name]]
  # Remove text annotations like "0 - Not at all" → "0"
  digestive[, (var_name) := {
    x <- get(cn)
    x <- str_extract(x, "^\\d+\\.?\\d*")
    suppressWarnings(as.numeric(x))
  }]
  vals <- digestive[[var_name]]
  cat("  ", var_name, ": range [", min(vals, na.rm = TRUE), ",",
      max(vals, na.rm = TRUE), "] n_valid=", sum(!is.na(vals)), "\n")
}

# --- Bowel movements per day (p21044, continuous count) ---
digestive[, bowel_movements_day := {
  x <- str_extract(p21044, "^\\d+\\.?\\d*")
  suppressWarnings(as.numeric(x))
}]
cat("  bowel_movements_day: range [",
    min(digestive$bowel_movements_day, na.rm = TRUE), ",",
    max(digestive$bowel_movements_day, na.rm = TRUE), "] n_valid=",
    sum(!is.na(digestive$bowel_movements_day)), "\n")

# --- PHQ-12 SS Somatic Items (p21051, p21053, p21059, p21060) ---
phq12_map <- c(
  "Not bothered at all"  = 0L,
  "Bothered a little"    = 1L,
  "Bothered a lot"       = 2L
)

phq12_items <- list(
  phq12_headache  = "p21051",
  phq12_dizziness = "p21053",
  phq12_nausea    = "p21059",
  phq12_fatigue   = "p21060"
)

for (var_name in names(phq12_items)) {
  cn <- phq12_items[[var_name]]
  digestive[, (var_name) := phq12_map[get(cn)]]
  cat("  ", var_name, ": 0=", sum(digestive[[var_name]] == 0, na.rm = TRUE),
      " 1=", sum(digestive[[var_name]] == 1, na.rm = TRUE),
      " 2=", sum(digestive[[var_name]] == 2, na.rm = TRUE),
      " NA=", sum(is.na(digestive[[var_name]])), "\n")
}

# --- DERIVED: somatic-symptom total (sum of 4 items, range 0-8) ---
phq12_var_names <- names(phq12_items)
digestive[, phq12_somatic_total := rowSums(.SD, na.rm = FALSE), 
   .SDcols = phq12_var_names]
# If any item is NA, total should be NA (strict complete-case for scale)
cat("  phq12_somatic_total: range [",
    min(digestive$phq12_somatic_total, na.rm = TRUE), ",",
    max(digestive$phq12_somatic_total, na.rm = TRUE), "] n_valid=",
    sum(!is.na(digestive$phq12_somatic_total)), "\n")

cat("\n--- Digestive health, step 3: Assembly ---\n")

# Define the final analysis variables
digestive_analysis_vars <- c(
  # IBS and core GI symptoms
  "ibs_ever", "abdom_discomfort_freq", "abdom_pain_6mo",
  "abdom_pain_severity", "abdom_distension",
  "bowel_satisfaction", "gi_life_interference", "bowel_movements_day",
  # Visceral sensitivity and history
  "sensitive_stomach", "ibs_family_history", "coeliac_gluten",
  # PHQ-12 SS somatic items
  "phq12_headache", "phq12_dizziness", "phq12_nausea", "phq12_fatigue",
  # Derived
  "phq12_somatic_total"
)

# Select only analysis variables (drop raw p-columns)
digestive_out <- digestive[, c("eid", digestive_analysis_vars), with = FALSE]

cat("  Output matrix:", nrow(digestive_out), "rows x",
    length(digestive_analysis_vars), "exposure variables\n")

# Participant-level: anyone with at least one non-NA digestive-health variable
digestive_respondents <- digestive_out[, .(has_data = any(!is.na(.SD))),
                          .SDcols = digestive_analysis_vars, by = eid]
n_respondents <- sum(digestive_respondents$has_data)
cat("  Digestive-health respondents (≥1 non-NA):", format(n_respondents, big.mark = ","),
    "(", round(100 * n_respondents / nrow(digestive_out), 1), "% of full cohort)\n")

cat("\n--- Digestive health, step 4: Per-group output ---\n")

digestive_var_dict <- rbindlist(c(
  lapply(names(digestive_registry), function(vn) {
    r <- digestive_registry[[vn]]
    dict_entry(vn, r$field, "Digestive health", "Cat 153",
               r$type, r$recode, r$relevance)
  }),
  list(dict_entry("phq12_somatic_total", NA, "Digestive health", "Cat 153",
                  "derived_continuous",
                  "Sum of phq12_headache + phq12_dizziness + phq12_nausea + phq12_fatigue (range 0-8); NA if any item missing",
                  "Composite somatic symptom burden; somatization shares central sensitization with taste threshold alterations"))
))

# Save variable dictionary (shared across groups)
fwrite(digestive_var_dict, file.path(OUTPUT_DIR, "followup_exwas_digestive_health_variable_dict.csv"))
cat("  Saved variable dictionary:", nrow(digestive_var_dict), "variables\n")

# Per-group: inner join with phenotype file, compute missingness, save
digestive_group_reports <- list()

for (gname in names(pheno_list)) {
  cat("\n  --- Group:", toupper(gname), "---\n")
  
  gpheno <- pheno_list[[gname]]
  
  # Inner join: keep only people in this ethnic group with taste phenotype
  gdata <- merge(digestive_out, gpheno[, .(eid, taste_2w_strict)], by = "eid")
  
  # Analysis sample: must have taste outcome
  gdata <- gdata[!is.na(taste_2w_strict)]
  n_cases    <- sum(gdata$taste_2w_strict == 1, na.rm = TRUE)
  n_controls <- sum(gdata$taste_2w_strict == 0, na.rm = TRUE)
  
  cat("    Analysis sample:", nrow(gdata),
      "(cases:", n_cases, "| controls:", n_controls, ")\n")
  
  # Missingness in this group's analysis sample
  miss_g <- missingness_report(gdata, digestive_analysis_vars)
  
  # Variables below missingness threshold
  n_pass <- sum(miss_g$pct_missing <= MISS_THRESHOLD_PCT)
  n_drop <- sum(miss_g$pct_missing > MISS_THRESHOLD_PCT)
  cat("    Variables ≤", MISS_THRESHOLD_PCT, "% missing:", n_pass, "\n")
  if (n_drop > 0) {
    cat("    Variables >", MISS_THRESHOLD_PCT, "% missing (flagged):", n_drop, "\n")
    print(miss_g[pct_missing > MISS_THRESHOLD_PCT])
  }
  
  # Save — include all variables; missingness filtering at regression time
  out_file <- file.path(OUTPUT_DIR,
                         paste0("followup_exwas_digestive_health_", gname, ".csv"))
  # Drop taste_2w_strict from output (lives in phenotype file)
  fwrite(gdata[, !"taste_2w_strict"], out_file)
  cat("    Saved:", out_file, "\n")
  
  miss_file <- file.path(OUTPUT_DIR,
                          paste0("followup_exwas_digestive_health_", gname, "_missingness.csv"))
  fwrite(miss_g, miss_file)
  
  digestive_group_reports[[gname]] <- list(
    n_merged = nrow(gdata), n_cases = n_cases, n_controls = n_controls,
    miss = miss_g
  )
}

cat("\n  Digestive-health cleaning complete.\n")

cat("\n================================================================\n")
cat("  Mental health — mental-health questionnaire (Category 136)\n")
cat("================================================================\n\n")

# --- Mental health, block a: PHQ-9 Depression Scale (9 items + 2 derived) ---
phq9_fields <- c(
  phq9_q1_interest    = "p20514",   # Lack of interest/pleasure
  phq9_q2_depressed   = "p20510",   # Feeling depressed
  phq9_q3_sleep       = "p20517",   # Sleep problems
  phq9_q4_tired       = "p20519",   # Tiredness/low energy
  phq9_q5_appetite    = "p20511",   # Poor appetite/overeating
  phq9_q6_inadequacy  = "p20507",   # Feelings of inadequacy
  phq9_q7_concentrate = "p20508",   # Trouble concentrating
  phq9_q8_movement    = "p20518",   # Changes in speed of moving/speaking
  phq9_q9_selfharm    = "p20513"    # Thoughts of suicide/self-harm
)

# PHQ-9 item encoding (UKB answer set 504; text labels scored 0-3)
phq9_map <- c(
  "Not at all"              = 0L,
  "Several days"            = 1L,
  "More than half the days" = 2L,
  "Nearly every day"        = 3L
)

# --- Mental health, block b: GAD-7 Anxiety Scale (7 items + 2 derived) ---
gad7_fields <- c(
  gad7_q1_nervous     = "p20506",   # Feeling nervous/anxious
  gad7_q2_control     = "p20509",   # Not being able to stop worrying
  gad7_q3_worry       = "p20520",   # Worrying too much
  gad7_q4_relax       = "p20515",   # Trouble relaxing
  gad7_q5_restless    = "p20516",   # Being restless
  gad7_q6_irritable   = "p20505",   # Becoming easily irritable
  gad7_q7_afraid      = "p20512"    # Feeling afraid
)

# --- Mental health, block c: AUDIT alcohol items (p20403-p20416; answer set differs per item) ---

# AUDIT Q1: Frequency (p20414)
audit_q1_map <- c(
  "Never"                        = 0L,
  "Monthly or less"              = 1L,
  "2 to 4 times a month"        = 2L,
  "2 to 3 times a week"         = 3L,
  "4 or more times a week"      = 4L
)

# AUDIT Q2: Typical quantity (p20403)
audit_q2_map <- c(
  "1 or 2"     = 0L,
  "3 or 4"     = 1L,
  "5 or 6"     = 2L,
  "7, 8 or 9"  = 3L,   # Note: comma in original label
  "7, 8, or 9" = 3L,   # Variant without Oxford comma
  "10 or more"  = 4L
)

# AUDIT Q3-Q8: Frequency-based (p20407-p20413, p20416)
audit_freq_map <- c(
  "Never"              = 0L,
  "Less than monthly"  = 1L,
  "Monthly"            = 2L,
  "Weekly"             = 3L,
  "Daily or almost daily" = 4L
)

# AUDIT Q9-Q10: Yes/No with time qualifier (p20405, p20411)
audit_yesno_map <- c(
  "No"                        = 0L,
  "Yes, but not in the last year" = 2L,
  "Yes, during the last year"     = 4L
)

# Full AUDIT field list
audit_fields <- c(
  audit_q1_freq     = "p20414",   # How often drink
  audit_q2_amount   = "p20403",   # Typical number of drinks
  audit_q3_6plus    = "p20416",   # 6+ drinks on one occasion
  audit_q4_control  = "p20413",   # Inability to cease drinking in last year
  audit_q5_fail     = "p20407",   # Failure to fulfil normal expectations
  audit_q6_morning  = "p20412",   # Needed morning drink after heavy session
  audit_q7_guilt    = "p20409",   # Guilt or remorse after drinking
  audit_q8_memory   = "p20408",   # Memory loss due to drinking
  audit_q9_injury   = "p20411",   # Injured self or someone else through drinking
  audit_q10_concern = "p20405"    # Known person concerned / recommended reduction
)

# --- Mental health, block d: Mania (Cat 139) ---
mania_fields <- c(
  mania_excitable     = "p20501",   # Ever had period of mania/excitability
  mania_irritable     = "p20502",   # Ever had period of extreme irritability
  mania_duration      = "p20492",   # Longest period
  mania_severity      = "p20493",   # Severity of problems
  mania_manifestation = "p20548"    # Manifestations (multi-select)
)

# --- Mental health, block e: Unusual Experiences (Cat 144) ---
unusual_fields <- c(
  unusual_exp_mania_dur  = "p20492",  # Already in mania — shared field
  unusual_exp_mania_sev  = "p20493"   # Already in mania — shared field
)

cat("  PHQ-9 items:", length(phq9_fields), "\n")
cat("  GAD-7 items:", length(gad7_fields), "\n")
cat("  AUDIT items:", length(audit_fields), "\n")
cat("  Mania items:", length(mania_fields), "\n")

# Verify all fields present
mental_all_cols <- unique(c(phq9_fields, gad7_fields, audit_fields, mania_fields,
                            unusual_fields))
register_fields("Mental health",
                sub("^p(\\d+).*$", "\\1", mental_all_cols))
mental_missing <- mental_all_cols[!mental_all_cols %in% names(raw)]
if (length(mental_missing) > 0) {
  cat("  WARNING: Missing columns:", paste(mental_missing, collapse = ", "), "\n")
} else {
  cat("  All declared fields found in data [OK]\n")
}

cat("\n--- Mental health, step 2: Recoding ---\n")

# Working copy
mental <- raw[, c("eid", unique(mental_all_cols)), with = FALSE]

# --- PHQ-9: recode 9 items → derive total + binary ---
for (var_name in names(phq9_fields)) {
  cn <- phq9_fields[[var_name]]
  mental[, (var_name) := phq9_map[get(cn)]]
}

phq9_item_vars <- names(phq9_fields)
mental[, phq9_total := rowSums(.SD, na.rm = FALSE), .SDcols = phq9_item_vars]
mental[, phq9_binary := fifelse(phq9_total >= 10, 1L, 0L)]
# If total is NA, binary is NA
mental[is.na(phq9_total), phq9_binary := NA_integer_]

cat("  PHQ-9 total: range [", min(mental$phq9_total, na.rm = TRUE), ",",
    max(mental$phq9_total, na.rm = TRUE), "] n_valid=",
    sum(!is.na(mental$phq9_total)), "\n")
cat("  PHQ-9 binary (≥10):", sum(mental$phq9_binary == 1, na.rm = TRUE),
    "positive /", sum(!is.na(mental$phq9_binary)), "valid\n")

# --- GAD-7: recode 7 items → derive total + binary ---
for (var_name in names(gad7_fields)) {
  cn <- gad7_fields[[var_name]]
  mental[, (var_name) := phq9_map[get(cn)]]  # Same answer set 504
}

gad7_item_vars <- names(gad7_fields)
mental[, gad7_total := rowSums(.SD, na.rm = FALSE), .SDcols = gad7_item_vars]
mental[, gad7_binary := fifelse(gad7_total >= 10, 1L, 0L)]
mental[is.na(gad7_total), gad7_binary := NA_integer_]

cat("  GAD-7 total: range [", min(mental$gad7_total, na.rm = TRUE), ",",
    max(mental$gad7_total, na.rm = TRUE), "] n_valid=",
    sum(!is.na(mental$gad7_total)), "\n")
cat("  GAD-7 binary (≥10):", sum(mental$gad7_binary == 1, na.rm = TRUE),
    "positive /", sum(!is.na(mental$gad7_binary)), "valid\n")

# --- AUDIT: recode 10 items with per-item maps → derive total + binary ---
mental[, audit_q1_freq := audit_q1_map[get(audit_fields["audit_q1_freq"])]]

# Q2 (amount)
mental[, audit_q2_amount := audit_q2_map[get(audit_fields["audit_q2_amount"])]]

# Q3-Q8 (frequency-based)
for (var_name in c("audit_q3_6plus", "audit_q4_control", "audit_q5_fail",
                    "audit_q6_morning", "audit_q7_guilt", "audit_q8_memory")) {
  cn <- audit_fields[[var_name]]
  mental[, (var_name) := audit_freq_map[get(cn)]]
}

# Q9-Q10 (yes/no with time qualifier)
for (var_name in c("audit_q9_injury", "audit_q10_concern")) {
  cn <- audit_fields[[var_name]]
  mental[, (var_name) := audit_yesno_map[get(cn)]]
}

audit_item_vars <- names(audit_fields)
mental[, audit_total := rowSums(.SD, na.rm = FALSE), .SDcols = audit_item_vars]
mental[, audit_binary := fifelse(audit_total >= 8, 1L, 0L)]
mental[is.na(audit_total), audit_binary := NA_integer_]

cat("  AUDIT total: range [", min(mental$audit_total, na.rm = TRUE), ",",
    max(mental$audit_total, na.rm = TRUE), "] n_valid=",
    sum(!is.na(mental$audit_total)), "\n")
cat("  AUDIT binary (≥8):", sum(mental$audit_binary == 1, na.rm = TRUE),
    "positive /", sum(!is.na(mental$audit_binary)), "valid\n")

# --- Mania: binary + ordinal + manifestation count ---
mental[, mania_excitable := recode_yesno(p20501)]
mental[, mania_irritable := recode_yesno(p20502)]

# Duration (ordinal)
mania_dur_map <- c(
  "Less than 24 hours"                   = 1L,
  "At least a day, but less than a week" = 2L,
  "A week or more"                       = 3L
)
mental[, mania_duration := mania_dur_map[p20492]]

# Severity
mania_sev_map <- c(
  "No problems" = 0L,
  "Needed treatment or caused problems with work, relationships, finances, the law or other aspects of life." = 1L
)
mental[, mania_severity := mania_sev_map[p20493]]

# Manifestations: pipe-delimited multi-select → count
mental[, mania_manifestation_count := {
  x <- p20548
  fifelse(
    is.na(x), NA_integer_,
    as.integer(str_count(x, "\\|") + 1L)
  )
}]

cat("  mania_excitable:", sum(mental$mania_excitable == 1, na.rm = TRUE), "Yes\n")
cat("  mania_irritable:", sum(mental$mania_irritable == 1, na.rm = TRUE), "Yes\n")
cat("  mania_manifestation_count: range [",
    min(mental$mania_manifestation_count, na.rm = TRUE), ",",
    max(mental$mania_manifestation_count, na.rm = TRUE), "] n_valid=",
    sum(!is.na(mental$mania_manifestation_count)), "\n")

cat("\n--- Mental health, step 3: Assembly ---\n")

mental_analysis_vars <- c(
  # PHQ-9
  "phq9_total", "phq9_binary",
  # GAD-7
  "gad7_total", "gad7_binary",
  # AUDIT
  "audit_total", "audit_binary",
  # Mania
  "mania_excitable", "mania_irritable",
  "mania_duration", "mania_severity", "mania_manifestation_count"
)

mental_out <- mental[, c("eid", mental_analysis_vars), with = FALSE]

cat("  Output matrix:", nrow(mental_out), "rows x",
    length(mental_analysis_vars), "exposure variables\n")

mental_respondents <- mental_out[, .(has_data = any(!is.na(.SD))),
                          .SDcols = mental_analysis_vars, by = eid]
cat("  Mental-health respondents:", format(sum(mental_respondents$has_data), big.mark = ","), "\n")

cat("\n--- Mental health, step 4: Per-group output ---\n")

mental_var_dict <- rbindlist(list(
  dict_entry("phq9_total", "20514+", "Mental health", "Cat 136/138",
             "derived_continuous", "Sum of 9 PHQ-9 items (each 0-3, range 0-27)",
             "Depression severity; serotonergic pathway → taste receptor modulation"),
  dict_entry("phq9_binary", "20514+", "Mental health", "Cat 136/138",
             "derived_binary", "PHQ-9 total ≥10",
             "Clinical depression threshold"),
  dict_entry("gad7_total", "20506+", "Mental health", "Cat 136/140",
             "derived_continuous", "Sum of 7 GAD-7 items (each 0-3, range 0-21)",
             "Anxiety severity; autonomic nervous system → salivary flow → taste"),
  dict_entry("gad7_binary", "20506+", "Mental health", "Cat 136/140",
             "derived_binary", "GAD-7 total ≥10",
             "Clinical anxiety threshold"),
  dict_entry("audit_total", "20414+", "Mental health", "Cat 136/142",
             "derived_continuous", "Sum of 10 AUDIT items (range 0-40)",
             "Problem drinking pattern; AUDIT ≠ the baseline alcohol-frequency exposure (different construct)"),
  dict_entry("audit_binary", "20414+", "Mental health", "Cat 136/142",
             "derived_binary", "AUDIT total ≥8",
             "Hazardous drinking threshold"),
  dict_entry("mania_excitable", 20501, "Mental health", "Cat 136/139",
             "binary", "Yes=1, No=0",
             "Bipolar spectrum; lithium/valproate → taste side effects"),
  dict_entry("mania_irritable", 20502, "Mental health", "Cat 136/139",
             "binary", "Yes=1, No=0",
             "Dysphoric mania; may have different taste association from euphoric"),
  dict_entry("mania_duration", 20492, "Mental health", "Cat 136/139",
             "ordinal", "<24h=1, <1wk=2, ≥1wk=3",
             "Dose-response: longer mania → more likely mood stabilizer treatment"),
  dict_entry("mania_severity", 20493, "Mental health", "Cat 136/139",
             "ordinal", "No problems=0, Treatment/problems=1",
             "Severity of mania consequences"),
  dict_entry("mania_manifestation_count", 20548, "Mental health", "Cat 136/139",
             "derived_continuous", "Count of pipe-separated manifestations (range 1-7+)",
             "Breadth of manic symptoms as severity proxy")
))

fwrite(mental_var_dict, file.path(OUTPUT_DIR, "followup_exwas_mental_health_variable_dict.csv"))

for (gname in names(pheno_list)) {
  cat("  --- Group:", toupper(gname), "---\n")
  gpheno <- pheno_list[[gname]]
  gdata <- merge(mental_out, gpheno[, .(eid, taste_2w_strict)], by = "eid")
  gdata <- gdata[!is.na(taste_2w_strict)]
  cat("    N=", nrow(gdata), " cases=",
      sum(gdata$taste_2w_strict == 1, na.rm = TRUE), "\n")
  
  out_file <- file.path(OUTPUT_DIR,
                         paste0("followup_exwas_mental_health_", gname, ".csv"))
  fwrite(gdata[, !"taste_2w_strict"], out_file)
  cat("    Saved:", out_file, "\n")
  
  miss_g <- missingness_report(gdata, mental_analysis_vars)
  fwrite(miss_g, file.path(OUTPUT_DIR,
                            paste0("followup_exwas_mental_health_", gname, "_missingness.csv")))
}

cat("\n  Mental-health cleaning complete.\n")

cat("\n================================================================\n")
cat("  follow-up ExWAS — Cleaning QC Report\n")
cat("================================================================\n\n")

report_file <- file.path(OUTPUT_DIR, "followup_exwas_cleaning_report.txt")
sink(report_file)

cat("================================================================\n")
cat("  Follow-up ExWAS — Data Cleaning QC Report\n")
cat("  Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

cat("--- INPUT ---\n")
cat("  Raw CSV:", nrow(raw), "participants x", ncol(raw), "columns\n")
cat("  Empty string → NA conversions:", format(n_converted_total, big.mark = ","), "\n")
cat("  Special missing labels → NA:", format(n_special_total, big.mark = ","), "\n\n")

cat("--- DOMAIN COLUMN COUNTS ---\n")
print(domain_summary)

cat("\n--- DIGESTIVE HEALTH ---\n")
cat("  Raw fields: 15 | Analysis variables: 16 (incl. 1 derived)\n")
for (gname in names(digestive_group_reports)) {
  gr <- digestive_group_reports[[gname]]
  cat("  ", toupper(gname), ": N=", gr$n_merged,
      " cases=", gr$n_cases, " controls=", gr$n_controls, "\n")
  high_miss <- gr$miss[pct_missing > MISS_THRESHOLD_PCT]
  if (nrow(high_miss) > 0) {
    cat("    High missingness variables:\n")
    print(high_miss[, .(variable, pct_missing)])
  }
}

cat("\n--- MENTAL HEALTH (MHQ1) ---\n")
cat("  Raw fields: 31 | Analysis variables: 11\n")
cat("  Instruments: PHQ-9 (total + binary), GAD-7 (total + binary),\n")
cat("               AUDIT (total + binary), Mania (5 vars)\n")

sink()
cat("  QC report saved:", report_file, "\n")

cat("\n================================================================\n")
cat("  Digestive health and mental health complete.\n")
cat("  `raw` and `pheno_list` stay in the session for 02 to 06.\n")
cat("================================================================\n")
