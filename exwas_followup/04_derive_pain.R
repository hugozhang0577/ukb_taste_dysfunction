#!/usr/bin/env Rscript
# ============================================================================
# Follow-up ExWAS — experience of pain (Category 154)
# ============================================================================
#
# Operates on the in-memory `raw` (+ `pheno_list`) loaded by
# 01_derive_mental_digestive.R.
#
# Input:   `raw` data.table (from shared framework Steps 1-3)
# Output:  followup_exwas_experience_of_pain_{group}.csv
#          followup_exwas_experience_of_pain_variable_dict.csv
#
# Strategy: Literature-validated instrument-based subscales + condition flags.
#
# Variables (16 total):
#   Chronic pain screen:    pain_chronic, pain_widespread
#   BPI:                    pain_avg (single item), pain_interference_mean (7-item subscale)
#   DN4 (UKB 7-item):       pain_dn4_score, pain_dn4_positive
#   Site-specific:          pain_orofacial, pain_n_sites
#   Pain conditions:        cond_neuropathy, cond_fibromyalgia, cond_cfs, cond_migraine
#   ACR-related severity:   pain_fatigue, pain_unrefreshed, pain_cog_symptoms
#   Global health:          pain_health_today
#
# Key references:
#   - UKB Pain Questionnaire v1.0 (Wynick, Smith, Bennett, Macfarlane)
#   - Bouhassira et al. 2005 Pain — DN4 development
#   - Cleeland & Ryan 1994 — BPI development
#   - Zelman et al. 2005 — BPI 2-factor structure (severity + interference)
#   - Liu et al. 2023 PMC — UKB neuropathic pain epidemiology
#   - Macfarlane et al. 2023 PAIN Reports — PAINSTORM/UKB phenotyping
#   - Tanguay-Sabourin et al. 2023 Nat Med — pain spreading risk score
#   - Wolfe et al. 2010 Arthritis Care Res — ACR fibromyalgia criteria
# ============================================================================

cat("\n================================================================\n")
cat("  Experience of pain (Category 154)\n")
cat("================================================================\n\n")

# ============================================================================
# 1. FIELD REGISTRY
# ============================================================================

# All 42 Cat 154 fields are from a single 2019 administration (no instances)
# Column naming: p{field_id} (no _i suffix)

mkcol <- function(fid) paste0("p", fid)

pain_fields <- list(
  # Screen
  chronic       = list(fid = 120019, col = mkcol(120019)),
  pain_length   = list(fid = 120020, col = mkcol(120020)),
  widespread    = list(fid = 120021, col = mkcol(120021)),
  
  # Pain conditions (binary)
  neuropathy    = list(fid = 120008, col = mkcol(120008)),
  fibromyalgia  = list(fid = 120009, col = mkcol(120009)),
  cfs           = list(fid = 120010, col = mkcol(120010)),
  migraine      = list(fid = 120016, col = mkcol(120016)),
  
  # BPI: average pain (single severity item)
  bpi_avg       = list(fid = 120088, col = mkcol(120088)),
  
  # BPI Interference subscale (7 items)
  bpi_int_activity = list(fid = 120091, col = mkcol(120091)),
  bpi_int_mood     = list(fid = 120092, col = mkcol(120092)),
  bpi_int_walking  = list(fid = 120093, col = mkcol(120093)),
  bpi_int_work     = list(fid = 120094, col = mkcol(120094)),
  bpi_int_relations= list(fid = 120095, col = mkcol(120095)),
  bpi_int_sleep    = list(fid = 120096, col = mkcol(120096)),
  bpi_int_enjoy    = list(fid = 120097, col = mkcol(120097)),
  
  # DN4 self-report items (UKB has 7 of 10 DN4 items)
  dn4_burning    = list(fid = 120046, col = mkcol(120046)),
  dn4_cold       = list(fid = 120047, col = mkcol(120047)),
  dn4_shocks     = list(fid = 120048, col = mkcol(120048)),
  dn4_tingling   = list(fid = 120049, col = mkcol(120049)),
  dn4_pins       = list(fid = 120050, col = mkcol(120050)),
  dn4_numbness   = list(fid = 120051, col = mkcol(120051)),
  dn4_itching    = list(fid = 120052, col = mkcol(120052)),
  
  # Body sites (14 sites — combined "presence + rating" fields)
  site_allover   = list(fid = 120022, col = mkcol(120022)),
  site_head      = list(fid = 120023, col = mkcol(120023)),
  site_face      = list(fid = 120024, col = mkcol(120024)),
  site_neck      = list(fid = 120025, col = mkcol(120025)),
  site_back      = list(fid = 120026, col = mkcol(120026)),
  site_stomach   = list(fid = 120027, col = mkcol(120027)),
  site_hip       = list(fid = 120028, col = mkcol(120028)),
  site_knee      = list(fid = 120029, col = mkcol(120029)),
  site_arm       = list(fid = 120030, col = mkcol(120030)),
  site_hand      = list(fid = 120031, col = mkcol(120031)),
  site_bothfeet  = list(fid = 120032, col = mkcol(120032)),
  site_feet      = list(fid = 120033, col = mkcol(120033)),
  site_leg       = list(fid = 120034, col = mkcol(120034)),
  site_chest     = list(fid = 120035, col = mkcol(120035)),
  site_other     = list(fid = 120036, col = mkcol(120036)),
  
  # ACR-related severity (Wolfe 2010 fibromyalgia criteria-aligned)
  acr_fatigue    = list(fid = 120040, col = mkcol(120040)),
  acr_unrefresh  = list(fid = 120041, col = mkcol(120041)),
  acr_cognitive  = list(fid = 120042, col = mkcol(120042)),
  
  # Treatment relief + global health
  treatment_relief = list(fid = 120090, col = mkcol(120090)),
  health_today     = list(fid = 120103, col = mkcol(120103))
)

cat("  Registered fields:", length(pain_fields), "\n")
register_fields("Experience of pain",
                vapply(pain_fields, function(x) x$fid, numeric(1)))

# Check column presence
all_pain_cols <- sapply(pain_fields, `[[`, "col")
pain_missing <- all_pain_cols[!all_pain_cols %in% names(raw)]
if (length(pain_missing) > 0) {
  cat("  WARNING: Missing columns (", length(pain_missing), "):\n")
  for (m in pain_missing) cat("    -", m, "\n")
} else {
  cat("  All declared columns found in data [OK]\n")
}

# Working copy — only present columns
present_cols <- intersect(all_pain_cols, names(raw))
pain <- raw[, c("eid", present_cols), with = FALSE]

# ============================================================================
# 2. HELPER: Coerce categorical Yes/No fields to 0/1
# ============================================================================

cat("\n--- Step 2: Yes/No coercion ---\n")

# Standard UKB Yes/No mapping (REPLACE mode text labels)
yn_map <- c(
  "Yes" = 1L, "yes" = 1L, "1" = 1L,
  "No" = 0L, "no" = 0L, "0" = 0L,
  "Do not know" = NA_integer_, "Prefer not to answer" = NA_integer_,
  "do not know" = NA_integer_, "prefer not to answer" = NA_integer_
)

recode_yn <- function(x) {
  if (is.numeric(x)) {
    # Already numeric — keep 0/1, NA otherwise
    return(fifelse(x %in% c(0L, 1L), as.integer(x), NA_integer_))
  }
  yn_map[as.character(x)]
}

# ============================================================================
# 3. SCREEN VARIABLES
# ============================================================================

cat("\n--- Step 3: Screen variables ---\n")

pain[, pain_chronic := recode_yn(get(pain_fields$chronic$col))]
cat("  pain_chronic: yes=", sum(pain$pain_chronic == 1, na.rm = TRUE),
    " no=", sum(pain$pain_chronic == 0, na.rm = TRUE), "\n")

pain[, pain_widespread := recode_yn(get(pain_fields$widespread$col))]
cat("  pain_widespread: yes=", sum(pain$pain_widespread == 1, na.rm = TRUE), "\n")

# ============================================================================
# 4. PAIN CONDITIONS (binary)
# ============================================================================

cat("\n--- Step 4: Pain conditions ---\n")

pain[, cond_neuropathy   := recode_yn(get(pain_fields$neuropathy$col))]
pain[, cond_fibromyalgia := recode_yn(get(pain_fields$fibromyalgia$col))]
pain[, cond_cfs          := recode_yn(get(pain_fields$cfs$col))]
pain[, cond_migraine     := recode_yn(get(pain_fields$migraine$col))]

for (cn in c("cond_neuropathy", "cond_fibromyalgia", "cond_cfs", "cond_migraine")) {
  n_yes <- sum(pain[[cn]] == 1, na.rm = TRUE)
  cat("  ", cn, ": yes=", n_yes, "\n")
}

# ============================================================================
# 5. BPI: AVERAGE PAIN (SINGLE SEVERITY ITEM)
# ============================================================================
#
# Note: UKB only collected 1 of 4 BPI severity items (average pain).
# We do NOT label this as "BPI severity score" — see Methods.
# ============================================================================

cat("\n--- Step 5: Average pain (BPI single item) ---\n")

pain[, pain_avg := suppressWarnings(as.numeric(get(pain_fields$bpi_avg$col)))]
pain[pain_avg < 0 | pain_avg > 10, pain_avg := NA_real_]
cat("  pain_avg: n_valid=", sum(!is.na(pain$pain_avg)),
    " mean=", round(mean(pain$pain_avg, na.rm = TRUE), 2), "\n")

# ============================================================================
# 6. BPI INTERFERENCE SUBSCALE (7-ITEM MEAN)
# ============================================================================
#
# Validated 7-item BPI interference subscale (Cleeland & Ryan 1994;
# Zelman et al. 2005 confirmed 2-factor structure).
# Subscale score = mean of valid items (requires ≥4/7 items per BPI manual).
# ============================================================================

cat("\n--- Step 6: BPI Interference subscale ---\n")

bpi_int_cols <- c(
  pain_fields$bpi_int_activity$col,
  pain_fields$bpi_int_mood$col,
  pain_fields$bpi_int_walking$col,
  pain_fields$bpi_int_work$col,
  pain_fields$bpi_int_relations$col,
  pain_fields$bpi_int_sleep$col,
  pain_fields$bpi_int_enjoy$col
)
bpi_int_cols <- intersect(bpi_int_cols, names(pain))

# Coerce to numeric (0-10 scale)
for (cn in bpi_int_cols) {
  pain[, (cn) := suppressWarnings(as.numeric(get(cn)))]
  pain[get(cn) < 0 | get(cn) > 10, (cn) := NA_real_]
}

# Compute mean (require ≥4 of 7 valid items, per BPI scoring guide)
pain[, bpi_n_valid := rowSums(!is.na(.SD)), .SDcols = bpi_int_cols]
pain[, pain_interference_mean := rowMeans(.SD, na.rm = TRUE),
   .SDcols = bpi_int_cols]
pain[bpi_n_valid < 4, pain_interference_mean := NA_real_]
pain[is.nan(pain_interference_mean), pain_interference_mean := NA_real_]
pain[, bpi_n_valid := NULL]

n_int <- sum(!is.na(pain$pain_interference_mean))
cat("  pain_interference_mean: n_valid=", n_int,
    " mean=", round(mean(pain$pain_interference_mean, na.rm = TRUE), 2), "\n")

# ============================================================================
# 7. DN4 SCORE (7 SELF-REPORT ITEMS)
# ============================================================================
#
# UKB collected 7 of 10 DN4 items (the 3 clinician-administered items
# for sensory examination are not in UKB self-report).
#
# Score = sum of "Yes" responses (range 0-7).
# Cutoff: ≥3/7 → "possible neuropathic pain" (Liu 2023, Macfarlane 2023).
# Note: Original DN4 cutoff is ≥4/10; UKB-adapted cutoff is ≥3/7.
# ============================================================================

cat("\n--- Step 7: DN4 score (UKB 7-item self-report) ---\n")

dn4_cols <- c(
  pain_fields$dn4_burning$col, pain_fields$dn4_cold$col,
  pain_fields$dn4_shocks$col, pain_fields$dn4_tingling$col,
  pain_fields$dn4_pins$col, pain_fields$dn4_numbness$col,
  pain_fields$dn4_itching$col
)
dn4_cols <- intersect(dn4_cols, names(pain))

# Recode each to 0/1
for (cn in dn4_cols) {
  pain[, (cn) := recode_yn(get(cn))]
}

pain[, dn4_n_valid := rowSums(!is.na(.SD)), .SDcols = dn4_cols]
pain[, pain_dn4_score := rowSums(.SD, na.rm = TRUE), .SDcols = dn4_cols]
# Require all 7 items present for valid score (DN4 is short, no imputation)
pain[dn4_n_valid < 7, pain_dn4_score := NA_integer_]

pain[, pain_dn4_positive := fifelse(
  is.na(pain_dn4_score), NA_integer_,
  fifelse(pain_dn4_score >= 3L, 1L, 0L)
)]
pain[, dn4_n_valid := NULL]

n_dn4 <- sum(!is.na(pain$pain_dn4_score))
n_dn4_pos <- sum(pain$pain_dn4_positive == 1, na.rm = TRUE)
cat("  pain_dn4_score: n_valid=", n_dn4, "\n")
cat("  pain_dn4_positive (≥3/7): n=", n_dn4_pos, "\n")

# ============================================================================
# 8. BODY SITE PROCESSING
# ============================================================================
#
# Site fields combine presence + rating in a single field. Encoding varies:
#   - May be "Yes" / "No" (categorical)
#   - May be numeric rating where 0/NA = no pain
#   - May be character "Yes, rating X"
#
# Strategy: treat any non-NA non-"No" non-zero value as positive.
# ============================================================================

cat("\n--- Step 8: Body site processing ---\n")

site_keys <- c("site_head", "site_face", "site_neck", "site_back",
               "site_stomach", "site_hip", "site_knee", "site_arm",
               "site_hand", "site_bothfeet", "site_feet", "site_leg",
               "site_chest", "site_other")

site_present_flag <- function(x) {
  if (is.numeric(x)) {
    return(fifelse(!is.na(x) & x > 0, 1L, fifelse(is.na(x), NA_integer_, 0L)))
  }
  # Character
  xc <- as.character(x)
  fifelse(
    is.na(xc) | xc %in% c("", "Prefer not to answer", "Do not know"),
    NA_integer_,
    fifelse(xc %in% c("No", "no", "0"), 0L, 1L)
  )
}

# Compute per-site presence flags (temporary)
site_flag_cols <- character(0)
for (sk in site_keys) {
  raw_col <- pain_fields[[sk]]$col
  if (!raw_col %in% names(pain)) next
  flag_col <- paste0("__flag_", sk)
  pain[, (flag_col) := site_present_flag(get(raw_col))]
  site_flag_cols <- c(site_flag_cols, flag_col)
}

# Orofacial pain: head OR face
head_flag_col <- "__flag_site_head"
face_flag_col <- "__flag_site_face"
head_avail <- head_flag_col %in% names(pain)
face_avail <- face_flag_col %in% names(pain)

if (head_avail && face_avail) {
  pain[, pain_orofacial := fifelse(
    is.na(get(head_flag_col)) & is.na(get(face_flag_col)),
    NA_integer_,
    fifelse(
      (get(head_flag_col) == 1L & !is.na(get(head_flag_col))) |
        (get(face_flag_col) == 1L & !is.na(get(face_flag_col))),
      1L, 0L
    )
  )]
} else if (head_avail) {
  pain[, pain_orofacial := get(head_flag_col)]
} else if (face_avail) {
  pain[, pain_orofacial := get(face_flag_col)]
} else {
  pain[, pain_orofacial := NA_integer_]
}
n_orofac <- sum(pain$pain_orofacial == 1, na.rm = TRUE)
cat("  pain_orofacial (head OR face): n=", n_orofac, "\n")

# Number of pain sites: count of positive site flags
pain[, pain_n_sites := rowSums(.SD == 1L, na.rm = TRUE), .SDcols = site_flag_cols]
# If ALL site flags are NA, set count to NA (not 0)
pain[, n_sites_valid := rowSums(!is.na(.SD)), .SDcols = site_flag_cols]
pain[n_sites_valid == 0, pain_n_sites := NA_integer_]
pain[, n_sites_valid := NULL]

cat("  pain_n_sites: n_valid=", sum(!is.na(pain$pain_n_sites)),
    " mean=", round(mean(pain$pain_n_sites, na.rm = TRUE), 2), "\n")

# Drop temporary site flags
pain[, (site_flag_cols) := NULL]

# ============================================================================
# 9. ACR-RELATED SEVERITY (FATIGUE / UNREFRESHED / COGNITIVE)
# ============================================================================
#
# These three items align with the Wolfe et al. 2010 ACR fibromyalgia
# preliminary diagnostic criteria Symptom Severity Scale (SS).
# Each rated 0-3 in the original ACR criteria; UKB scale may differ.
# Treated as continuous ordinal exposures.
# ============================================================================

cat("\n--- Step 9: ACR-related severity items ---\n")

acr_recode <- function(x) {
  if (is.numeric(x)) return(x)
  # Common UKB severity scale labels
  m <- c(
    "No problem"     = 0,
    "Slight problem" = 1,
    "Moderate problem" = 2,
    "Severe problem" = 3,
    "0" = 0, "1" = 1, "2" = 2, "3" = 3,
    "Prefer not to answer" = NA_real_,
    "Do not know" = NA_real_
  )
  m[as.character(x)]
}

pain[, pain_fatigue       := acr_recode(get(pain_fields$acr_fatigue$col))]
pain[, pain_unrefreshed   := acr_recode(get(pain_fields$acr_unrefresh$col))]
pain[, pain_cog_symptoms  := acr_recode(get(pain_fields$acr_cognitive$col))]

cat("  pain_fatigue: n_valid=", sum(!is.na(pain$pain_fatigue)), "\n")
cat("  pain_unrefreshed: n_valid=", sum(!is.na(pain$pain_unrefreshed)), "\n")
cat("  pain_cog_symptoms: n_valid=", sum(!is.na(pain$pain_cog_symptoms)), "\n")

# ============================================================================
# 10. GLOBAL HEALTH TODAY
# ============================================================================

cat("\n--- Step 10: Global health today ---\n")

pain[, pain_health_today := suppressWarnings(
  as.numeric(get(pain_fields$health_today$col))
)]
pain[pain_health_today < 0 | pain_health_today > 100, pain_health_today := NA_real_]
cat("  pain_health_today (0-100): n_valid=", sum(!is.na(pain$pain_health_today)), "\n")

# ============================================================================
# 11. ASSEMBLE + PER-GROUP OUTPUT
# ============================================================================

cat("\n--- Step 11: Assembly ---\n")

pain_analysis_vars <- c(
  "pain_chronic", "pain_widespread",
  "pain_avg", "pain_interference_mean",
  "pain_dn4_score", "pain_dn4_positive",
  "pain_orofacial", "pain_n_sites",
  "cond_neuropathy", "cond_fibromyalgia", "cond_cfs", "cond_migraine",
  "pain_fatigue", "pain_unrefreshed", "pain_cog_symptoms",
  "pain_health_today"
)

pain_out <- pain[, c("eid", pain_analysis_vars), with = FALSE]

cat("  Output:", nrow(pain_out), "rows x",
    length(pain_analysis_vars), "variables\n")

pain_respondents <- pain_out[, .(has_data = any(!is.na(.SD))),
                         .SDcols = pain_analysis_vars, by = eid]
cat("  Respondents:", format(sum(pain_respondents$has_data), big.mark = ","), "\n")

pain_var_dict <- rbindlist(list(
  dict_entry("pain_chronic", 120019, "Experience of pain", "Cat 154", "binary",
             "Pain or discomfort >3 months", "Chronic pain screen"),
  dict_entry("pain_widespread", 120021, "Experience of pain", "Cat 154", "binary",
             "Pain all over body in last 3 months",
             "Widespread pain (sensitization phenotype)"),
  dict_entry("pain_avg", 120088, "Experience of pain", "Cat 154", "continuous",
             "Average pain rating (0-10)",
             "BPI single severity item (UKB only collected 1/4 BPI severity items)"),
  dict_entry("pain_interference_mean", "120091-120097", "Experience of pain", "Cat 154",
             "continuous",
             "Mean of 7 BPI Interference items (0-10); requires ≥4 valid items",
             "Validated BPI Interference subscale (Cleeland & Ryan 1994; Zelman 2005)"),
  dict_entry("pain_dn4_score", "120046-120052", "Experience of pain", "Cat 154",
             "continuous",
             "Sum of 7 DN4 self-report neuropathic items (0-7); all 7 required",
             "DN4 self-report subset (Bouhassira 2005; UKB excludes 3 clinician items)"),
  dict_entry("pain_dn4_positive", "120046-120052", "Experience of pain", "Cat 154",
             "binary",
             "DN4 self-report ≥3/7 cutoff",
             "Possible neuropathic pain per Liu 2023; Macfarlane 2023"),
  dict_entry("pain_orofacial", "120023+120024", "Experience of pain", "Cat 154",
             "binary",
             "Head OR face pain in last 3 months",
             "Trigeminal-mediated pain — anatomically relevant to gustatory pathways"),
  dict_entry("pain_n_sites", "120023-120036", "Experience of pain", "Cat 154",
             "continuous",
             "Count of body sites with pain in last 3 months (0-14)",
             "Pain spreading index (Tanguay-Sabourin 2023 Nat Med)"),
  dict_entry("cond_neuropathy", 120008, "Experience of pain", "Cat 154", "binary",
             "Ever had nerve damage/neuropathy (non-diabetic)",
             "Direct sensory nerve pathology"),
  dict_entry("cond_fibromyalgia", 120009, "Experience of pain", "Cat 154", "binary",
             "Ever had fibromyalgia",
             "Central sensitization syndrome — known association with taste dysfunction"),
  dict_entry("cond_cfs", 120010, "Experience of pain", "Cat 154", "binary",
             "Ever had CFS/ME", "Post-viral/central syndrome"),
  dict_entry("cond_migraine", 120016, "Experience of pain", "Cat 154", "binary",
             "Ever had migraine",
             "Trigeminal sensory disorder — shares cranial nerve pathways"),
  dict_entry("pain_fatigue", 120040, "Experience of pain", "Cat 154", "ordinal",
             "Fatigue severity past week (0-3)",
             "Aligned with Wolfe 2010 ACR fibromyalgia symptom severity scale"),
  dict_entry("pain_unrefreshed", 120041, "Experience of pain", "Cat 154", "ordinal",
             "Waking unrefreshed past week (0-3)",
             "ACR fibromyalgia SS item"),
  dict_entry("pain_cog_symptoms", 120042, "Experience of pain", "Cat 154", "ordinal",
             "Cognitive symptoms severity (0-3)",
             "ACR fibromyalgia SS item"),
  dict_entry("pain_health_today", 120103, "Experience of pain", "Cat 154", "continuous",
             "Self-rated health today (0-100)",
             "EQ-5D-VAS-like global health rating")
))

fwrite(pain_var_dict, file.path(OUTPUT_DIR, "followup_exwas_experience_of_pain_variable_dict.csv"))
cat("  Saved variable dictionary:", nrow(pain_var_dict), "variables\n")

for (gname in names(pheno_list)) {
  cat("  --- Group:", toupper(gname), "---\n")
  gpheno <- pheno_list[[gname]]
  gdata <- merge(pain_out, gpheno[, .(eid, taste_2w_strict)], by = "eid")
  gdata <- gdata[!is.na(taste_2w_strict)]
  n_cases <- sum(gdata$taste_2w_strict == 1, na.rm = TRUE)
  cat("    N=", nrow(gdata), " cases=", n_cases, "\n")
  
  out_file <- file.path(OUTPUT_DIR,
                        paste0("followup_exwas_experience_of_pain_", gname, ".csv"))
  fwrite(gdata[, !"taste_2w_strict"], out_file)
  cat("    Saved:", out_file, "\n")
  
  miss_g <- missingness_report(gdata, pain_analysis_vars)
  fwrite(miss_g, file.path(OUTPUT_DIR,
                           paste0("followup_exwas_experience_of_pain_", gname, "_missingness.csv")))
}

cat("\n  Pain cleaning complete.\n")
