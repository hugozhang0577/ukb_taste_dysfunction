#!/usr/bin/env Rscript
# ============================================================================
# Follow-up ExWAS — diet by 24-hour recall (Categories 100117 and 100112)
# ============================================================================
#
# Operates on the in-memory `raw` (+ `pheno_list`) loaded by
# 01_derive_mental_digestive.R.
#
# Input:   `raw` data.table (from shared framework Steps 1-3)
#          `pheno_list` (for sex variable needed in energy exclusion)
# Output:  followup_exwas_diet_24_hour_recall_{group}.csv
#          followup_exwas_diet_24_hour_recall_variable_dict.csv
#
# Strategy (Option A — literature-driven, minimal work):
#   Tier 1: Cat 100117 derived nutrients (61 fields × 5 instances → ~50 vars)
#   Tier 3: Cat 100112 supplements (p20084 + p104670 → ~2 vars)
#   Tier 2 (Cat 100118 food groups) SKIPPED — see Methods justification.
#
# Processing pipeline (per Greenwood 2019, Perez-Cornago 2021):
#   1. Per-instance extreme energy exclusion (sex-specific thresholds)
#   2. Cross-instance mean (habitual intake estimation)
#   3. Flag single-instance records
#   4. High-correlation filter (|r| > 0.8)
#
# Key references:
#   - Greenwood et al. 2019 Am J Epidemiol — biomarker validation of WebQ
#   - Perez-Cornago et al. 2021 Eur J Nutr — updated nutrient calculation
#   - Liu et al. 2011 Public Health Nutr — Oxford WebQ development
#   - Galante et al. 2016 Br J Nutr — repeat WebQ acceptability
# ============================================================================

cat("\n================================================================\n")
cat("  Diet (24-hour recall) — nutrients (Category 100117), supplements (Category 100112)\n")
cat("================================================================\n\n")

# ============================================================================
# CONFIGURATION
# ============================================================================

# Energy field for extreme intake exclusion
# Field 26002 = "Energy" in kJ/day from Cat 100117
ENERGY_FIELD_ID <- 26002

# Sex-specific extreme energy thresholds in kJ (Perez-Cornago et al. 2021)
# Men:   <3347 or >17573 kJ/day  (equivalent to <800 or >4200 kcal/day)
# Women: <2092 or >14644 kJ/day  (equivalent to <500 or >3500 kcal/day)
ENERGY_EXCL <- list(
  male   = c(low = 3347, high = 17573),
  female = c(low = 2092, high = 14644)
)

# High-correlation threshold for nutrient de-redundancy
CORR_THRESHOLD <- 0.80

# Number of instances
N_INSTANCES <- 5  # i0 through i4

# Nutrient field → name lookup (from UKB Showcase Cat 100117).
# A lookup, not the export specification: the nutrients actually processed are
# whichever Category 100117 fields the export carries (see below). The lookup
# is deliberately a superset, so a nutrient can never reach the dictionary
# unnamed; field_list.txt reports any entry it does not export.
NUTRIENT_NAMES <- c(
  "26000" = "total_weight_food_bev_g",
  "26001" = "total_weight_bev_g",
  "26002" = "energy_kj",
  "26003" = "energy_from_bev_kj",
  "26004" = "energy_density",
  "26005" = "protein_g",
  "26006" = "vegetable_protein_g",
  "26007" = "animal_protein_g",
  "26008" = "fat_g",
  "26009" = "vegetable_fat_g",
  "26010" = "animal_fat_g",
  "26011" = "total_sugars_g",
  "26012" = "free_sugar_g",
  "26013" = "carbohydrate_g",
  "26014" = "saturated_fa_g",
  "26015" = "n3_fa_g",
  "26016" = "n6_fa_g",
  "26017" = "englyst_fibre_g",
  "26018" = "calcium_mg",
  "26019" = "iron_mg",
  "26020" = "vitamin_b6_mg",
  "26021" = "vitamin_b12_ug",
  "26022" = "folate_ug",
  "26023" = "vitamin_c_mg",
  "26024" = "potassium_mg",
  "26025" = "magnesium_mg",
  "26026" = "retinol_ug",
  "26027" = "total_carotene_ug",
  "26028" = "vitamin_e_mg",
  "26029" = "vitamin_d_ug",
  "26030" = "alcohol_g",
  "26031" = "starch_g",
  "26032" = "mufa_g",
  "26033" = "zinc_mg",
  "26034" = "thiamin_mg",
  "26035" = "riboflavin_mg",
  "26036" = "phosphorus_mg",
  "26037" = "cholesterol_mg",
  "26038" = "alpha_carotene_ug",
  "26039" = "beta_carotene_ug",
  "26040" = "beta_cryptoxanthin_ug",
  "26041" = "biotin_ug",
  "26042" = "chloride_mg",
  "26043" = "copper_mg",
  "26044" = "fructose_g",
  "26045" = "glucose_g",
  "26046" = "haem_iron_mg",
  "26047" = "iodine_ug",
  "26048" = "lactose_g",
  "26049" = "maltose_g",
  "26050" = "intrinsic_milk_sugars_g",
  "26051" = "manganese_mg",
  "26052" = "sodium_mg",
  "26053" = "non_haem_iron_mg",
  "26054" = "niacin_eq_mg",
  "26055" = "nmes_g",
  "26056" = "other_sugars_g",
  "26057" = "pantothenic_acid_mg",
  "26058" = "selenium_ug",
  "26059" = "sucrose_g",
  "26060" = "total_nitrogen_g",
  "26061" = "vitamin_a_re_ug",
  "26155" = "trans_fa_g"
)

# ============================================================================
# 1. IDENTIFY TIER 1 (NUTRIENT) AND TIER 3 (SUPPLEMENT) COLUMNS
# ============================================================================

cat("--- Step 1: Column identification ---\n")

# Tier 1: Cat 100117 derived nutrients, field ids 26000-26061.
# The id range, not the family name, is what selects them: the family also
# carries the supplement fields (character, multi-select) and — where an export
# includes them — the Category 100118 food-group fields, and neither belongs in
# the cross-instance numeric averaging below.
NUTRIENT_ID_RANGE <- c(26000, 26061)
t1_field_ids <- sort(unique(
  col_map[domain == "Diet (24-hour recall)" &
            field_id >= NUTRIENT_ID_RANGE[1] &
            field_id <= NUTRIENT_ID_RANGE[2], field_id]
))
cat("  Tier 1 nutrient fields:", length(t1_field_ids), "\n")
register_fields("Diet (24-hour recall) nutrients",
                as.integer(names(NUTRIENT_NAMES)), lookup_only = TRUE)

# Build instance column matrix: rows = fields, cols = instances
t1_inst_cols <- matrix(NA_character_, nrow = length(t1_field_ids), ncol = N_INSTANCES)
rownames(t1_inst_cols) <- as.character(t1_field_ids)
colnames(t1_inst_cols) <- paste0("i", 0:(N_INSTANCES - 1))

for (i in seq_along(t1_field_ids)) {
  fid <- t1_field_ids[i]
  for (j in 0:(N_INSTANCES - 1)) {
    cn <- paste0("p", fid, "_i", j)
    if (cn %in% names(raw)) {
      t1_inst_cols[i, j + 1] <- cn
    }
  }
}

n_found <- sum(!is.na(t1_inst_cols))
cat("  Instance columns found:", n_found, "/", length(t1_inst_cols), "\n")

# Tier 3: Supplements
# p20084 = multi-select supplement list (pipe-delimited, 5 instances)
# p104670 = binary "Did you take any vitamin/mineral supplements?" (5 instances)
SUPPLEMENT_FIELD_IDS <- c(20084, 104670)
register_fields("Diet (24-hour recall) supplements", SUPPLEMENT_FIELD_IDS)

t3_supp_cols <- paste0("p", SUPPLEMENT_FIELD_IDS[1], "_i", 0:(N_INSTANCES - 1))
t3_any_cols  <- paste0("p", SUPPLEMENT_FIELD_IDS[2], "_i", 0:(N_INSTANCES - 1))

t3_supp_present <- t3_supp_cols[t3_supp_cols %in% names(raw)]
t3_any_present  <- t3_any_cols[t3_any_cols %in% names(raw)]

cat("  Tier 3 supplement list columns:", length(t3_supp_present), "\n")
cat("  Tier 3 any-supplement columns:", length(t3_any_present), "\n")

# ============================================================================
# 2. EXTRACT WORKING COPY + MERGE SEX
# ============================================================================

cat("\n--- Step 2: Working copy + sex merge ---\n")

# Collect all 24-hour-recall columns
all_diet_cols <- c(
  as.vector(t1_inst_cols[!is.na(t1_inst_cols)]),
  t3_supp_present,
  t3_any_present
)

diet <- raw[, c("eid", all_diet_cols), with = FALSE]

# Merge sex from primary phenotype file (needed for energy exclusion)
# sex must be coded: 0 = female, 1 = male (or "Female"/"Male")
pheno_primary <- pheno_list[["group1"]]
if (!"sex" %in% names(pheno_primary)) {
  # Try common alternatives
  sex_col <- intersect(c("sex", "Sex", "genetic_sex", "p31"), names(pheno_primary))
  if (length(sex_col) > 0) {
    cat("  Sex column found as:", sex_col[1], "\n")
  } else {
    cat("  WARNING: No sex column found — energy exclusion will use combined thresholds\n")
  }
}

# Build sex lookup from ALL group phenotype files
sex_lookup <- rbindlist(lapply(pheno_list, function(p) {
  sx <- intersect(c("sex", "Sex", "genetic_sex", "p31"), names(p))
  if (length(sx) > 0) {
    p[, .(eid, sex_val = get(sx[1]))]
  } else {
    data.table(eid = integer(0), sex_val = character(0))
  }
}), fill = TRUE)
sex_lookup <- unique(sex_lookup, by = "eid")

# Standardise sex coding
sex_lookup[, is_male := fcase(
  sex_val %in% c(1, "1", "Male", "male"), TRUE,
  sex_val %in% c(0, "0", "Female", "female"), TRUE,
  default = NA
)]
sex_lookup[, is_male := fcase(
  sex_val %in% c(1, "1", "Male", "male"), TRUE,
  sex_val %in% c(0, "0", "Female", "female"), FALSE,
  default = NA
)]

diet <- merge(diet, sex_lookup[, .(eid, is_male)], by = "eid", all.x = TRUE)
cat("  Sex merged:", sum(!is.na(diet$is_male)), "/", nrow(diet), "participants\n")

# ============================================================================
# 3. PER-INSTANCE EXTREME ENERGY EXCLUSION
# ============================================================================
#
# For each instance: if total energy is outside sex-specific bounds,
# set ALL nutrient values for that instance to NA.
# This must happen BEFORE cross-instance averaging (Perez-Cornago 2021).
# ============================================================================

cat("\n--- Step 3: Per-instance extreme energy exclusion ---\n")

n_excluded_total <- 0L

for (inst in 0:(N_INSTANCES - 1)) {
  energy_col <- paste0("p", ENERGY_FIELD_ID, "_i", inst)
  
  if (!energy_col %in% names(diet)) next
  
  # Determine per-person threshold based on sex
  # Male: <800 or >4200; Female: <500 or >3500
  energy_vals <- diet[[energy_col]]
  
  extreme_flag <- fifelse(
    diet$is_male == TRUE,
    !is.na(energy_vals) & (energy_vals < ENERGY_EXCL$male["low"] |
                             energy_vals > ENERGY_EXCL$male["high"]),
    fifelse(
      diet$is_male == FALSE,
      !is.na(energy_vals) & (energy_vals < ENERGY_EXCL$female["low"] |
                               energy_vals > ENERGY_EXCL$female["high"]),
      FALSE  # Unknown sex: don't exclude
    )
  )
  
  n_extreme <- sum(extreme_flag, na.rm = TRUE)
  n_excluded_total <- n_excluded_total + n_extreme
  
  if (n_extreme > 0) {
    # Set ALL nutrient columns for this instance to NA
    inst_cols <- t1_inst_cols[, inst + 1]
    inst_cols <- inst_cols[!is.na(inst_cols)]
    
    for (cn in inst_cols) {
      diet[extreme_flag == TRUE, (cn) := NA_real_]
    }
    
    cat("  Instance", inst, ": excluded", n_extreme, "extreme energy records\n")
  }
}

cat("  Total instance-records excluded:", n_excluded_total, "\n")

# ============================================================================
# 4. CROSS-INSTANCE AVERAGING
# ============================================================================
#
# For each nutrient field: compute mean across all valid (non-NA) instances.
# This reduces within-person day-to-day variation and better approximates
# habitual intake (Greenwood et al. 2019).
#
# Also record: number of valid instances per person (for QC flagging).
# ============================================================================

cat("\n--- Step 4: Cross-instance averaging ---\n")

# Compute mean and instance count for each nutrient field
nutrient_mean_names <- sapply(as.character(t1_field_ids), function(fid) {
  nm <- NUTRIENT_NAMES[fid]
  if (is.na(nm)) paste0("nutr_", fid) else paste0("nutr_", nm)
})
names(nutrient_mean_names) <- NULL
n_inst_col <- "nutr_n_valid_instances"

# First: count valid instances using energy field as proxy
energy_inst_cols <- paste0("p", ENERGY_FIELD_ID, "_i", 0:(N_INSTANCES - 1))
energy_inst_cols <- intersect(energy_inst_cols, names(diet))

diet[, (n_inst_col) := rowSums(!is.na(.SD)), .SDcols = energy_inst_cols]

cat("  Instance count distribution:\n")
print(diet[get(n_inst_col) > 0, .N, by = get(n_inst_col)][order(get)])

# Compute per-nutrient mean across instances
for (i in seq_along(t1_field_ids)) {
  fid <- t1_field_ids[i]
  inst_cols <- t1_inst_cols[i, ]
  inst_cols <- inst_cols[!is.na(inst_cols)]
  
  if (length(inst_cols) == 0) next
  
  mean_name <- nutrient_mean_names[i]
  diet[, (mean_name) := rowMeans(.SD, na.rm = TRUE), .SDcols = inst_cols]
  # rowMeans with all NA → NaN; convert to NA
  diet[is.nan(get(mean_name)), (mean_name) := NA_real_]
}

n_with_diet <- sum(diet[[n_inst_col]] > 0)
cat("  Participants with ≥1 valid dietary recall:", 
    format(n_with_diet, big.mark = ","), "\n")
cat("  Nutrient variables computed:", length(nutrient_mean_names), "\n")

# ============================================================================
# 5. HIGH-CORRELATION FILTER
# ============================================================================
#
# Remove one variable from each pair with |r| > 0.8 to reduce redundancy.
# Retain the variable with lower missingness (or arbitrary if equal).
# This prevents near-collinear exposures inflating ExWAS Stage 2.
# ============================================================================

cat("\n--- Step 5: High-correlation filter (|r| >", CORR_THRESHOLD, ") ---\n")

# Compute correlation on participants with data
nutr_mat <- as.matrix(diet[get(n_inst_col) > 0, ..nutrient_mean_names])
nutr_cor <- cor(nutr_mat, use = "pairwise.complete.obs")

# Find correlated pairs
high_cor_pairs <- which(abs(nutr_cor) > CORR_THRESHOLD & upper.tri(nutr_cor),
                        arr.ind = TRUE)

# Track which variables to drop
drop_vars <- character(0)

if (nrow(high_cor_pairs) > 0) {
  cat("  High-correlation pairs found:", nrow(high_cor_pairs), "\n")
  
  for (k in seq_len(nrow(high_cor_pairs))) {
    v1 <- nutrient_mean_names[high_cor_pairs[k, 1]]
    v2 <- nutrient_mean_names[high_cor_pairs[k, 2]]
    r_val <- nutr_cor[high_cor_pairs[k, 1], high_cor_pairs[k, 2]]
    
    # Skip if one already marked for drop
    if (v1 %in% drop_vars || v2 %in% drop_vars) next
    
    # Drop the one with more missingness
    miss1 <- sum(is.na(diet[[v1]]))
    miss2 <- sum(is.na(diet[[v2]]))
    to_drop <- ifelse(miss1 >= miss2, v1, v2)
    drop_vars <- c(drop_vars, to_drop)
    
    cat("    |r| =", round(abs(r_val), 3), ":", v1, "vs", v2,
        "→ drop", to_drop, "\n")
  }
}

nutrient_retained <- setdiff(nutrient_mean_names, drop_vars)
cat("  Nutrients retained:", length(nutrient_retained),
    "/ dropped:", length(drop_vars), "\n")

# ============================================================================
# 6. SUPPLEMENT PROCESSING (TIER 3)
# ============================================================================
#
# p104670: "Did you take any vitamin/mineral supplements yesterday?"
#   → any_supplement_ever: 1 if "Yes" in ANY instance
#
# p20084: pipe-delimited supplement list
#   → Extract key taste-relevant supplements as binary flags
#   → "ever reported" across all instances
# ============================================================================

cat("\n--- Step 6: Supplement processing ---\n")

# --- p104670: any supplement binary ---
diet[, supp_any_ever := {
  vals <- .SD
  any_yes <- apply(vals, 1, function(row) {
    any(row == "Yes", na.rm = TRUE)
  })
  fifelse(
    any_yes, 1L,
    fifelse(
      apply(vals, 1, function(row) all(is.na(row))), NA_integer_,
      0L
    )
  )
}, .SDcols = t3_any_present]

cat("  supp_any_ever: Yes=", sum(diet$supp_any_ever == 1, na.rm = TRUE),
    " No=", sum(diet$supp_any_ever == 0, na.rm = TRUE),
    " NA=", sum(is.na(diet$supp_any_ever)), "\n")

# --- p20084: specific supplements ---
# Target supplements with taste-pathway relevance:
#   Zinc → gustin/CA-VI → taste bud maintenance
#   Vitamin B12 → nerve integrity
#   Folic acid → mucosal cell turnover
#   Iron → taste receptor function
#   Vitamin D → immune-mediated taste effects
#   Calcium → taste signal transduction

target_supplements <- c(
  "Zinc", "Vitamin B12", "Folic acid", "Iron",
  "Vitamin D", "Calcium"
)

# Search across all instances for each target
for (supp_name in target_supplements) {
  var_name <- paste0("supp_", tolower(gsub("[^a-zA-Z0-9]", "_", supp_name)))
  
  diet[, (var_name) := {
    vals <- .SD
    # Check if supplement name appears in any instance's pipe-delimited list
    found <- apply(vals, 1, function(row) {
      any(grepl(supp_name, row, fixed = TRUE), na.rm = TRUE)
    })
    # Determine if person has any supplement data at all
    has_data <- apply(vals, 1, function(row) any(!is.na(row)))
    
    fifelse(found, 1L,
            fifelse(has_data, 0L, NA_integer_))
  }, .SDcols = t3_supp_present]
  
  n_yes <- sum(diet[[var_name]] == 1, na.rm = TRUE)
  cat("  ", var_name, ":", n_yes, "users\n")
}

# ============================================================================
# 7. ASSEMBLE FINAL EXPOSURE MATRIX
# ============================================================================

cat("\n--- Step 7: Assembly ---\n")

# Supplement variable names
supp_vars <- c(
  "supp_any_ever",
  paste0("supp_", tolower(gsub("[^a-zA-Z0-9]", "_", target_supplements)))
)

# Single-instance flag (for sensitivity analysis)
diet[, nutr_single_instance := fifelse(get(n_inst_col) == 1, 1L, 0L)]

diet_analysis_vars <- c(nutrient_retained, supp_vars)

diet_out <- diet[, c("eid", diet_analysis_vars, n_inst_col, "nutr_single_instance"),
             with = FALSE]

cat("  Output matrix:", nrow(diet_out), "rows x",
    length(diet_analysis_vars), "exposure variables\n")
cat("    Nutrients:", length(nutrient_retained),
    "| Supplements:", length(supp_vars), "\n")

diet_respondents <- diet_out[get(n_inst_col) > 0]
cat("  Respondents (≥1 valid recall):",
    format(nrow(diet_respondents), big.mark = ","), "\n")

# ============================================================================
# 8. PER-GROUP OUTPUT
# ============================================================================

cat("\n--- Step 8: Per-group output ---\n")

# Variable dictionary
diet_var_dict <- rbindlist(c(
  # Nutrient variables
  lapply(seq_along(nutrient_retained), function(i) {
    vn <- nutrient_retained[i]
    # Reverse-lookup field ID from nutrient name
    fid <- t1_field_ids[which(nutrient_mean_names %in% vn)]
    fid_str <- if (length(fid) > 0) fid[1] else NA
    desc <- NUTRIENT_NAMES[as.character(fid_str)]
    if (is.na(desc)) desc <- vn
    dict_entry(
      vn, fid_str, "Diet (24-hour recall)", "Cat 100117",
      "continuous",
      paste0("Mean across valid 24h-recall instances (range 1-5); ",
             "extreme energy days excluded per Perez-Cornago 2021"),
      paste0("Dietary nutrient: ", desc)
    )
  }),
  # Supplement variables
  list(
    dict_entry("supp_any_ever", 104670, "Diet (24-hour recall)", "Cat 100112",
               "binary", "1 if 'Yes' in any instance, 0 if all 'No', NA if no data",
               "Any vitamin/mineral supplement use"),
    dict_entry("supp_zinc", 20084, "Diet (24-hour recall)", "Cat 100112",
               "binary", "1 if zinc reported in any instance",
               "DIRECT: Zn → gustin/carbonic anhydrase VI → taste bud turnover"),
    dict_entry("supp_vitamin_b12", 20084, "Diet (24-hour recall)", "Cat 100112",
               "binary", "1 if Vitamin B12 reported in any instance",
               "DIRECT: B12 → nerve integrity → gustatory nerve function"),
    dict_entry("supp_folic_acid", 20084, "Diet (24-hour recall)", "Cat 100112",
               "binary", "1 if Folic acid reported in any instance",
               "Folate → mucosal cell turnover → taste papillae renewal"),
    dict_entry("supp_iron", 20084, "Diet (24-hour recall)", "Cat 100112",
               "binary", "1 if Iron reported in any instance",
               "Iron → taste receptor protein synthesis"),
    dict_entry("supp_vitamin_d", 20084, "Diet (24-hour recall)", "Cat 100112",
               "binary", "1 if Vitamin D reported in any instance",
               "Vitamin D → immune modulation → inflammatory taste effects"),
    dict_entry("supp_calcium", 20084, "Diet (24-hour recall)", "Cat 100112",
               "binary", "1 if Calcium reported in any instance",
               "Ca2+ → taste signal transduction cascade")
  )
))

fwrite(diet_var_dict, file.path(OUTPUT_DIR, "followup_exwas_diet_24_hour_recall_variable_dict.csv"))
cat("  Saved variable dictionary:", nrow(diet_var_dict), "variables\n")

# Save correlation filter log
if (length(drop_vars) > 0) {
  corr_log <- data.table(dropped_variable = drop_vars,
                         reason = "High correlation (|r| > 0.8) with retained variable")
  fwrite(corr_log, file.path(OUTPUT_DIR, "followup_exwas_diet_24_hour_recall_correlation_drops.csv"))
  cat("  Saved correlation filter log\n")
}

# Per-group output
for (gname in names(pheno_list)) {
  cat("  --- Group:", toupper(gname), "---\n")
  gpheno <- pheno_list[[gname]]
  gdata <- merge(diet_out, gpheno[, .(eid, taste_2w_strict)], by = "eid")
  gdata <- gdata[!is.na(taste_2w_strict)]
  n_cases <- sum(gdata$taste_2w_strict == 1, na.rm = TRUE)
  n_with <- sum(gdata[[n_inst_col]] > 0)
  cat("    N=", nrow(gdata), " cases=", n_cases,
      " with_diet_data=", n_with, "\n")
  
  out_file <- file.path(OUTPUT_DIR,
                        paste0("followup_exwas_diet_24_hour_recall_", gname, ".csv"))
  fwrite(gdata[, !"taste_2w_strict"], out_file)
  cat("    Saved:", out_file, "\n")
  
  miss_g <- missingness_report(gdata, diet_analysis_vars)
  fwrite(miss_g, file.path(OUTPUT_DIR,
                           paste0("followup_exwas_diet_24_hour_recall_", gname, "_missingness.csv")))
}

cat("\n  Diet-recall cleaning complete.\n")
cat("  Nutrients:", length(nutrient_retained), "after correlation filter\n")
cat("  Supplements:", length(supp_vars), "\n")
cat("  Extreme energy exclusions:", n_excluded_total, "instance-records\n")
