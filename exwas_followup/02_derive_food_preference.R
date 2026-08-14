#!/usr/bin/env Rscript
# ============================================================================
# Follow-up ExWAS — food preferences (Category 1039)
# ============================================================================
#
# Operates on the in-memory `raw` (+ `pheno_list`) loaded by
# 01_derive_mental_digestive.R.
#
# Input:   `raw` data.table (loaded by 01_derive_mental_digestive.R)
# Output:  followup_exwas_food_preferences_{group}.csv
#          followup_exwas_food_preferences_variable_dict.csv
#          followup_exwas_food_preferences_pca_loadings.csv
#
# Strategy (Option A — conservative, literature-supported):
#   1. Taste-quality items (7 items): FPQ items that EXPLICITLY name a
#      gustatory/trigeminal quality. Enter ExWAS individually as continuous 1-9.
#   2. General food preference items (133 items): All remaining food items.
#      PCA → PC1-PC5 enter ExWAS as continuous composite scores.
#   3. Non-food items excluded from export (already absent from data).
#
# Literature support:
#   - 7 taste-quality items are self-evident probes: the questionnaire text
#     directly names the sensory quality (e.g., "Liking for bitter foods").
#   - PCA on remaining items follows May-Wilson et al. (2022, Nat Commun)
#     data-driven approach to food-liking dimensionality.
#   - Fatty taste inclusion: Keast & Costanzo (2015) Prog Lipid Res.
#   - Spicy/chemesthetic: Lawless & Heymann (2010) Sensory Evaluation.
# ============================================================================

cat("\n================================================================\n")
cat("  Food preferences — food-preference questionnaire (Category 1039)\n")
cat("================================================================\n\n")

# ============================================================================
# 1. FIELD REGISTRY + TASTE-QUALITY CLASSIFICATION
# ============================================================================
#
# Classification principle (Option A):
#   taste_quality = TRUE:  ONLY items whose questionnaire text explicitly names
#                          a gustatory or chemesthetic quality (bitter, sweet,
#                          salty, fatty, spicy). These are unambiguous taste
#                          modality probes requiring no inferential classification.
#   taste_quality = FALSE: All other food/beverage items — integrated flavour
#                          experiences reflecting taste + smell + texture +
#                          cultural factors → PCA reduction.
#
# This conservative strategy avoids the methodological risk of classifying
# specific foods (e.g., dark chocolate, grapefruit) as single-modality probes
# without consumer-panel validation (cf. Bawajeeh et al. 2022; TasteLQ 2024).
# ============================================================================

# --- All 140 FPQ food item field IDs and descriptions ---
fpq_field_ids <- c(
  20600, 20601, 20602, 20603, 20604, 20605, 20606, 20607, 20608, 20609,
  20610, 20611, 20612, 20613, 20615, 20616, 20617, 20618, 20619, 20620,
  20621, 20622, 20623, 20624, 20625, 20626, 20627, 20628, 20629, 20630,
  20631, 20632, 20633, 20634, 20635, 20636, 20637, 20638, 20639, 20640,
  20642, 20643, 20644, 20645, 20646, 20647, 20648, 20649, 20650, 20651,
  20652, 20653, 20654, 20655, 20658, 20659, 20660, 20661, 20662, 20663,
  20664, 20665, 20666, 20667, 20671, 20672, 20673, 20674, 20675, 20676,
  20677, 20678, 20679, 20680, 20681, 20682, 20683, 20684, 20685, 20686,
  20687, 20688, 20689, 20690, 20691, 20692, 20693, 20694, 20695, 20696,
  20697, 20698, 20699, 20700, 20701, 20702, 20703, 20704, 20705, 20706,
  20707, 20708, 20709, 20710, 20711, 20712, 20713, 20714, 20715, 20716,
  20717, 20718, 20719, 20720, 20721, 20722, 20723, 20724, 20725, 20726,
  20727, 20728, 20729, 20730, 20731, 20732, 20734, 20735, 20736, 20737,
  20738, 20739, 20740, 20742, 20743, 20744, 20745, 20746, 20747, 20748
)

fpq_descriptions <- c(
  "adding salt to foods", "aniseed", "apple juice", "apples", "asparagus",
  "aubergine", "avocados", "bacon", "baked/steamed fish", "bananas",
  "barbequed/grilled meat", "beef steak", "beetroot", "bell pepper",
  "biscuits", "bitter foods", "bitter/ale", "black olives", "black pepper",
  "blue cheese", "bolognese sauce", "broad beans", "broccoli", "brown rice",
  "brussel sprouts", "burgers (meat)", "burn of spicy foods", "butter on bread",
  "butternut squash", "cabbage", "cake", "cake icing", "capers", "cauliflower",
  "cereal/granola bar", "cheesecake", "cherries", "chicken",
  "chilli pepper", "chips/french fries", "cod", "coffee with sugar",
  "coffee without sugar", "coriander", "corn flakes", "cream", "croissant",
  "cucumber", "curry", "dairy products", "dark chocolate",
  "diet fizzy drinks", "dried fruit", "eggs", "extra virgin olive oil",
  "fatty foods", "fresh tomatoes", "fried chicken", "fried/battered fish",
  "fruit", "garlic", "gherkins", "globe artichoke", "goat's cheese",
  "grapefruit", "green olives", "haddock", "ham", "hard cheese",
  "herring", "honey", "horseradish/wasabi", "ice cream", "jam",
  "kiwi fruit", "lager", "lamb", "lemons", "lentils/beans", "liver",
  "mackerel", "marzipan", "mayonnaise", "melon", "milk chocolate",
  "mushrooms", "onions", "orange juice", "oranges", "pasta", "pears",
  "pizza", "plain yogurt", "plums", "pollock", "pork chop", "porridge",
  "potato crisps", "potatoes", "prawns", "raw carrots", "red meat",
  "red wine", "regular fizzy drinks", "roast chicken", "salad dressing",
  "salad leaves", "salami", "salmon", "salty foods", "salty pretzels",
  "sardines", "sausages (meat)", "savoury biscuits", "shellfish",
  "skimmed milk", "smoked fish", "soft cheese", "soy sauce", "soya milk",
  "spicy foods", "spinach", "spirits", "strawberries",
  "sweet coffee house drinks", "sweet foods", "tea with sugar",
  "tea without sugar", "tomato ketchup", "tinned tuna", "turnip",
  "vegetables", "vinegar", "whisky", "white bread", "white rice",
  "white wine", "whole grain breakfast cereal", "whole milk", "wholemeal bread"
)

# --- 7 taste-quality items: explicit gustatory/chemesthetic quality in text ---
taste_quality_ids <- c(
  20616,  # "Liking for bitter foods"          → bitter
  20732,  # "Liking for sweet foods"           → sweet
  20716,  # "Liking for salty foods"            → salty
  20600,  # "Liking for adding salt to foods"   → salty
  20659,  # "Liking for fatty foods"            → fatty
  20727,  # "Liking for spicy foods"            → spicy
  20627   # "Liking for burn of spicy foods"    → spicy/chemesthetic
)

taste_quality_modality <- c(
  "bitter", "sweet", "salty", "salty", "fatty", "spicy", "spicy"
)

taste_quality_names <- c(
  "fpq_taste_bitter",       # 20616
  "fpq_taste_sweet",        # 20732
  "fpq_taste_salty",        # 20716
  "fpq_taste_salt_adding",  # 20600
  "fpq_taste_fatty",        # 20659
  "fpq_taste_spicy",        # 20727
  "fpq_taste_spicy_burn"    # 20627
)

# Build registry
fpq_fields <- data.table(
  field_id    = fpq_field_ids,
  description = fpq_descriptions,
  colname     = paste0("p", fpq_field_ids),
  taste_quality = fpq_field_ids %in% taste_quality_ids
)

n_taste <- sum(fpq_fields$taste_quality)
n_pca   <- sum(!fpq_fields$taste_quality)
register_fields("Food preferences", fpq_field_ids)

cat("  Total FPQ food items:", nrow(fpq_fields), "\n")
cat("  Taste-quality items (enter ExWAS individually):", n_taste, "\n")
cat("  General food items (PCA reduction):", n_pca, "\n")

# Verify all columns present
fpq_missing <- fpq_fields$colname[!fpq_fields$colname %in% names(raw)]
if (length(fpq_missing) > 0) {
  cat("  WARNING: Missing columns:", paste(fpq_missing, collapse = ", "), "\n")
} else {
  cat("  All declared fields found in data [OK]\n")
}

# ============================================================================
# 2. NUMERIC RECODING (Data-Coding 96)
# ============================================================================
#
# The 1-9 hedonic scale has partial REPLACE encoding:
#   "extremely dislike" = 1,  "2"-"4" = as-is,
#   "neither like nor dislike" = 5,  "6"-"8" = as-is,
#   "extremely like" = 9
# Special: "never tried" → NA  (Data-Coding -121)
#          "do not wish to answer" → NA  (Data-Coding -818, already cleaned)
# ============================================================================

cat("\n--- Step 2: Numeric recoding (Data-Coding 96) ---\n")

fpq_recode_map <- c(
  "extremely dislike"         = 1L,
  "2"                         = 2L,
  "3"                         = 3L,
  "4"                         = 4L,
  "neither like nor dislike"  = 5L,
  "6"                         = 6L,
  "7"                         = 7L,
  "8"                         = 8L,
  "extremely like"            = 9L,
  "never tried"               = NA_integer_
)

# Working copy: eid + all food-preference columns
fpq <- raw[, c("eid", fpq_fields$colname), with = FALSE]

# Recode all 140 columns to integer 1-9
n_never_tried <- 0L
for (cn in fpq_fields$colname) {
  n_nt <- sum(fpq[[cn]] == "never tried", na.rm = TRUE)
  n_never_tried <- n_never_tried + n_nt
  fpq[, (cn) := fpq_recode_map[get(cn)]]
}

cat("  Recoded", nrow(fpq_fields), "columns to integer 1-9\n")
cat("  'never tried' → NA:", format(n_never_tried, big.mark = ","),
    "total entries across all items\n")

# Validation: check all values are in 1-9 or NA
range_check <- fpq[, lapply(.SD, function(x) {
  r <- range(x, na.rm = TRUE)
  if (any(is.infinite(r))) return(c(NA, NA))
  r
}), .SDcols = fpq_fields$colname]
cat("  Global value range: [", min(unlist(range_check), na.rm = TRUE), ",",
    max(unlist(range_check), na.rm = TRUE), "]\n")

# ============================================================================
# 3. TASTE-QUALITY VARIABLE EXTRACTION
# ============================================================================

cat("\n--- Step 3: Taste-quality variable extraction ---\n")

taste_core_cols <- paste0("p", taste_quality_ids)

# Create semantically named copies
for (i in seq_along(taste_quality_ids)) {
  cn <- paste0("p", taste_quality_ids[i])
  fpq[, (taste_quality_names[i]) := get(cn)]
  n_valid <- sum(!is.na(fpq[[taste_quality_names[i]]]))
  cat("  ", taste_quality_names[i], " (p", taste_quality_ids[i], " — ",
      fpq_fields[field_id == taste_quality_ids[i], description],
      "): n_valid=", n_valid, "\n", sep = "")
}

# ============================================================================
# 4. PCA ON GENERAL FOOD PREFERENCE ITEMS
# ============================================================================
#
# PCA on the 133 non-taste-quality food items to extract composite dietary
# preference pattern scores (PC1-PC5). This data-driven approach follows
# May-Wilson et al. (2022, Nat Commun) who used genetic correlations +
# structural equation modelling to identify hierarchical food-liking
# dimensions from the same UKB FPQ data.
#
# Technical notes:
#   - Uses pairwise-complete correlation matrix (avoids listwise deletion
#     which would lose participants with any "never tried" response)
#   - Eigen decomposition on correlation matrix (equivalent to PCA on
#     standardized variables)
#   - PCA computed on full cohort; per-group projection uses same loadings
# ============================================================================

cat("\n--- Step 4: PCA on general food preference items ---\n")

pca_pool_cols <- fpq_fields[taste_quality == FALSE, colname]
pca_mat <- as.matrix(fpq[, ..pca_pool_cols])

cat("  PCA input:", ncol(pca_mat), "variables\n")
cat("  Complete cases:", format(sum(complete.cases(pca_mat)), big.mark = ","), "\n")

# Check for zero-variance columns (can occur if all responses are "never tried")
col_vars <- apply(pca_mat, 2, var, na.rm = TRUE)
zero_var_cols <- names(col_vars)[is.na(col_vars) | col_vars < 1e-10]
if (length(zero_var_cols) > 0) {
  cat("  Dropping", length(zero_var_cols), "zero-variance columns from PCA\n")
  pca_pool_cols <- setdiff(pca_pool_cols, zero_var_cols)
  pca_mat <- as.matrix(fpq[, ..pca_pool_cols])
}

# Pairwise-complete correlation matrix
pca_cor <- cor(pca_mat, use = "pairwise.complete.obs")

# Check for NA in correlation matrix (occurs when a pair has no common obs)
n_na_cor <- sum(is.na(pca_cor))
if (n_na_cor > 0) {
  cat("  WARNING:", n_na_cor, "NA values in correlation matrix;",
      "replacing with 0\n")
  pca_cor[is.na(pca_cor)] <- 0
}

# Eigen decomposition
pca_eigen <- eigen(pca_cor)
n_pc <- 5L

# Explained variance
var_explained <- pca_eigen$values / sum(pca_eigen$values) * 100
cum_var <- cumsum(var_explained)

cat("  Variance explained:\n")
for (k in 1:n_pc) {
  cat("    PC", k, ": ", round(var_explained[k], 1), "%",
      " (cumulative: ", round(cum_var[k], 1), "%)\n", sep = "")
}

# Compute PC scores via standardized projection
col_means <- colMeans(pca_mat, na.rm = TRUE)
col_sds   <- apply(pca_mat, 2, sd, na.rm = TRUE)

pca_std <- sweep(pca_mat, 2, col_means, "-")
pca_std <- sweep(pca_std, 2, col_sds, "/")

loadings <- pca_eigen$vectors[, 1:n_pc]
pc_scores <- pca_std %*% loadings

# Add PC scores to fpq
pc_names <- paste0("fpq_pca_pc", 1:n_pc)
for (k in 1:n_pc) {
  fpq[, (pc_names[k]) := pc_scores[, k]]
}

cat("  PC scores computed:", paste(pc_names, collapse = ", "), "\n")

# Save loadings for reproducibility
pca_loadings_dt <- data.table(
  variable = pca_pool_cols,
  description = fpq_fields[colname %in% pca_pool_cols, description]
)
for (k in 1:n_pc) {
  pca_loadings_dt[, (paste0("PC", k)) := round(loadings[, k], 5)]
}
fwrite(pca_loadings_dt,
       file.path(OUTPUT_DIR, "followup_exwas_food_preferences_pca_loadings.csv"))
cat("  Saved PCA loadings\n")

# Save PCA parameters for Methods reporting
pca_params <- data.table(
  PC = paste0("PC", 1:n_pc),
  eigenvalue = round(pca_eigen$values[1:n_pc], 3),
  var_explained_pct = round(var_explained[1:n_pc], 1),
  cum_var_pct = round(cum_var[1:n_pc], 1)
)
fwrite(pca_params,
       file.path(OUTPUT_DIR, "followup_exwas_food_preferences_pca_variance.csv"))

# ============================================================================
# 5. ASSEMBLE FINAL EXPOSURE MATRIX
# ============================================================================

cat("\n--- Step 5: Assembly ---\n")

fpq_analysis_vars <- c(taste_quality_names, pc_names)

fpq_out <- fpq[, c("eid", fpq_analysis_vars), with = FALSE]

cat("  Output matrix:", nrow(fpq_out), "rows x",
    length(fpq_analysis_vars), "exposure variables\n")
cat("    Taste-quality:", length(taste_quality_names),
    "| PCA:", length(pc_names), "\n")

fpq_respondents <- fpq_out[, .(has_data = any(!is.na(.SD))),
                         .SDcols = fpq_analysis_vars, by = eid]
cat("  Respondents:", format(sum(fpq_respondents$has_data), big.mark = ","), "\n")

# ============================================================================
# 6. PER-GROUP OUTPUT
# ============================================================================

cat("\n--- Step 6: Per-group output ---\n")

# Variable dictionary
fpq_var_dict <- rbindlist(c(
  # 7 taste-quality items
  lapply(seq_along(taste_quality_ids), function(i) {
    dict_entry(
      taste_quality_names[i],
      taste_quality_ids[i],
      "Food preferences",
      "Cat 1039",
      "continuous",
      "1-9 hedonic scale (1=extremely dislike, 9=extremely like); 'never tried'=NA",
      paste0("Direct taste-quality probe: ",
             taste_quality_modality[i], " — ",
             fpq_fields[field_id == taste_quality_ids[i], description])
    )
  }),
  # PCA scores
  lapply(1:n_pc, function(k) {
    dict_entry(
      pc_names[k],
      NA,
      "Food preferences",
      "Cat 1039",
      "derived_continuous",
      paste0("PC", k, " from PCA on ", length(pca_pool_cols),
             " general food preference items; variance explained: ",
             round(var_explained[k], 1), "%"),
      "Composite dietary preference pattern (data-driven, cf. May-Wilson 2022)"
    )
  })
))

fwrite(fpq_var_dict, file.path(OUTPUT_DIR, "followup_exwas_food_preferences_variable_dict.csv"))
cat("  Saved variable dictionary:", nrow(fpq_var_dict), "variables\n")

# Per-group output
for (gname in names(pheno_list)) {
  cat("  --- Group:", toupper(gname), "---\n")
  gpheno <- pheno_list[[gname]]
  gdata <- merge(fpq_out, gpheno[, .(eid, taste_2w_strict)], by = "eid")
  gdata <- gdata[!is.na(taste_2w_strict)]
  n_cases <- sum(gdata$taste_2w_strict == 1, na.rm = TRUE)
  cat("    N=", nrow(gdata), " cases=", n_cases, "\n")
  
  out_file <- file.path(OUTPUT_DIR,
                        paste0("followup_exwas_food_preferences_", gname, ".csv"))
  fwrite(gdata[, !"taste_2w_strict"], out_file)
  cat("    Saved:", out_file, "\n")
  
  miss_g <- missingness_report(gdata, fpq_analysis_vars)
  fwrite(miss_g, file.path(OUTPUT_DIR,
                           paste0("followup_exwas_food_preferences_", gname, "_missingness.csv")))
}

cat("\n  Food-preference cleaning complete.\n")
cat("  Taste-quality items:", length(taste_quality_names), "\n")
cat("  PCA components:", n_pc,
    "(explaining", round(cum_var[n_pc], 1), "% variance)\n")
