#!/usr/bin/env Rscript
# Follow-up questionnaire exposures, block 3: 24-hour dietary recall. Nutrient intakes (Category
# 100117) are averaged over valid recalls after sex-specific energy exclusion and de-duplicated at
# |r| > 0.8; supplement use (Category 100112) is any report across recalls. Runs in the R session of 01.
library(data.table)
library(stringr)

NUTRIENTS <- c(
  "26000" = "total_weight_food_bev_g", "26001" = "total_weight_bev_g",
  "26002" = "energy_kj", "26003" = "energy_from_bev_kj", "26004" = "energy_density",
  "26005" = "protein_g", "26006" = "vegetable_protein_g", "26007" = "animal_protein_g",
  "26008" = "fat_g", "26009" = "vegetable_fat_g", "26010" = "animal_fat_g",
  "26011" = "total_sugars_g", "26012" = "free_sugar_g", "26013" = "carbohydrate_g",
  "26014" = "saturated_fa_g", "26015" = "n3_fa_g", "26016" = "n6_fa_g",
  "26017" = "englyst_fibre_g", "26018" = "calcium_mg", "26019" = "iron_mg",
  "26020" = "vitamin_b6_mg", "26021" = "vitamin_b12_ug", "26022" = "folate_ug",
  "26023" = "vitamin_c_mg", "26024" = "potassium_mg", "26025" = "magnesium_mg",
  "26026" = "retinol_ug", "26027" = "total_carotene_ug", "26028" = "vitamin_e_mg",
  "26029" = "vitamin_d_ug", "26030" = "alcohol_g", "26031" = "starch_g",
  "26032" = "mufa_g", "26033" = "zinc_mg", "26034" = "thiamin_mg",
  "26035" = "riboflavin_mg", "26036" = "phosphorus_mg", "26037" = "cholesterol_mg",
  "26038" = "alpha_carotene_ug", "26039" = "beta_carotene_ug",
  "26040" = "beta_cryptoxanthin_ug", "26041" = "biotin_ug", "26042" = "chloride_mg",
  "26043" = "copper_mg", "26044" = "fructose_g", "26045" = "glucose_g",
  "26046" = "haem_iron_mg", "26047" = "iodine_ug", "26048" = "lactose_g",
  "26049" = "maltose_g", "26050" = "intrinsic_milk_sugars_g", "26051" = "manganese_mg",
  "26052" = "sodium_mg", "26053" = "non_haem_iron_mg", "26054" = "niacin_eq_mg",
  "26055" = "nmes_g", "26056" = "other_sugars_g", "26057" = "pantothenic_acid_mg",
  "26058" = "selenium_ug", "26059" = "sucrose_g", "26060" = "total_nitrogen_g",
  "26061" = "vitamin_a_re_ug")
# energy intake (kJ) outside these limits marks an implausible recall
ENERGY_LIMITS <- list(male = c(3347, 17573), female = c(2092, 14644))
SUPPLEMENTS <- c(supp_zinc = "Zinc", supp_vitamin_b12 = "Vitamin B12",
                 supp_folic_acid = "Folic acid", supp_iron = "Iron",
                 supp_vitamin_d = "Vitamin D", supp_calcium = "Calcium")

# nutrient fields present in the extract, with their recall-instance columns (i0-i4)
ids <- as.integer(str_match(names(raw), "^p(\\d+)(?:_i\\d+)?(?:_a\\d+)?$")[, 2])
nutr_ids <- sort(unique(ids[!is.na(ids) & ids >= 26000 & ids <= 26061]))
nutr_cols <- setNames(lapply(nutr_ids, function(id) intersect(paste0("p", id, "_i", 0:4), names(raw))),
                      paste0("nutr_", NUTRIENTS[as.character(nutr_ids)]))
supp_list_cols <- paste0("p20084_i", 0:4)
supp_any_cols  <- paste0("p104670_i", 0:4)

diet <- raw[, c("eid", unlist(nutr_cols), supp_list_cols, supp_any_cols), with = FALSE]
sex <- unique(rbindlist(lapply(pheno_list, function(p) p[, .(eid, sex)])), by = "eid")
diet <- merge(diet, sex[, .(eid, is_male = fcase(sex == 1, TRUE, sex == 0, FALSE))],
              by = "eid", all.x = TRUE)

# Drop every nutrient of a recall whose energy intake is implausible (unknown sex: kept)
lo <- fifelse(diet$is_male, ENERGY_LIMITS$male[1], ENERGY_LIMITS$female[1])
hi <- fifelse(diet$is_male, ENERGY_LIMITS$male[2], ENERGY_LIMITS$female[2])
for (j in 0:4) {
  energy <- diet[[paste0("p26002_i", j)]]
  diet[which(energy < lo | energy > hi),
       (intersect(paste0("p", nutr_ids, "_i", j), names(diet))) := NA_real_]
}

diet[, n_valid_recalls := rowSums(!is.na(.SD)), .SDcols = patterns("^p26002_i")]
for (v in names(nutr_cols))
  diet[, (v) := { m <- rowMeans(.SD, na.rm = TRUE); fifelse(is.nan(m), NA_real_, m) },
       .SDcols = nutr_cols[[v]]]

# De-duplicate correlated nutrients: of each pair with |r| > 0.8, drop the one with more missing values
nutr <- names(nutr_cols)
R <- cor(as.matrix(diet[n_valid_recalls > 0, ..nutr]), use = "pairwise.complete.obs")
pairs <- which(abs(R) > 0.8 & upper.tri(R), arr.ind = TRUE)
drop_vars <- character(0)
for (k in seq_len(nrow(pairs))) {
  v <- nutr[pairs[k, ]]
  if (any(v %in% drop_vars)) next
  drop_vars <- c(drop_vars, if (sum(is.na(diet[[v[1]]])) >= sum(is.na(diet[[v[2]]]))) v[1] else v[2])
}
if (length(drop_vars) > 0)
  fwrite(data.table(dropped_variable = drop_vars,
                    reason = "High correlation (|r| > 0.8) with retained variable"),
         file.path(OUT_DIR, "followup_exwas_diet_24_hour_recall_correlation_drops.csv"))

# Supplements: 1 if reported in any recall, 0 if answered without it, NA if never answered
any_m  <- as.matrix(diet[, ..supp_any_cols])
list_m <- as.matrix(diet[, ..supp_list_cols])
diet[, supp_any_ever := fifelse(rowSums(any_m == "Yes", na.rm = TRUE) > 0, 1L,
                                fifelse(rowSums(!is.na(any_m)) == 0, NA_integer_, 0L))]
for (v in names(SUPPLEMENTS)) {
  found <- rowSums(matrix(grepl(SUPPLEMENTS[[v]], list_m, fixed = TRUE), nrow(list_m))) > 0
  diet[, (v) := fifelse(found, 1L, fifelse(rowSums(!is.na(list_m)) > 0, 0L, NA_integer_))]
}

keep <- setdiff(nutr, drop_vars)
save_section(diet, "Diet (24-hour recall)",
             c(setNames(rep("continuous", length(keep)), keep),
               setNames(rep("binary", 1 + length(SUPPLEMENTS)), c("supp_any_ever", names(SUPPLEMENTS)))))
