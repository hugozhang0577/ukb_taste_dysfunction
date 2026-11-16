# Build the master feature manifest from the discovery cohort (group1): FDR-significant
# results of the exposure, diagnosis, protein and metabolite scans, the APOE e4 dose,
# and the declared covariates, outcome and smell items.
library(data.table)
setwd(Sys.getenv("PROJECT_DIR", "."))

IN_RESULTS <- "input/assoc_results"
IN_READY   <- "input/analysis_ready"
OUT_DIR    <- "output/feature_manifest"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Per-cohort value sources; {group} is expanded to group1/2/3 at extraction
RAW_PROTEOMICS   <- file.path(IN_READY, "proteomics_{group}.csv")
RAW_METABOLOMICS <- file.path(IN_READY, "metabolomics_{group}.csv")
RAW_EXWAS_A      <- file.path(IN_READY, "exwas_baseline_{group}.csv")
RAW_EXWAS_B      <- file.path(IN_READY, "exwas_followup_{group}.csv")
RAW_DIAGNOSIS    <- file.path(IN_READY, "icd10_matrix_{group}.rds")
RAW_PHENOTYPE    <- file.path(IN_READY, "phenotype_{group}.csv")
RAW_APOE         <- file.path(IN_READY, "apoe_{group}.csv")

FDR_THRESHOLD  <- 0.05
FEMALE_ANCHORS <- c("age_menarche", "age_menopause", "menopause")

# data_requirement: data a model needs to use the feature
#   A = questionnaire and routine biochemistry (A_female: female-only items),
#   B = recorded diagnoses, C = genotype, D = plasma proteins / metabolic biomarkers

# Declared covariates, outcome and smell items
L0 <- data.table(
  feature_id = c("age_baseline", "age", "sex", "townsend", "BMI",
                 "smoking", "drink", "assess_centre_id",
                 paste0("PC", 1:10), "taste_2w_strict",
                 "smell_any", "smell_time", "smell_extent"),
  data_requirement = c(rep("A", 8), rep("C", 10), "outcome", rep("A", 3)),
  var_role = c(rep("covariate_demographic", 5), rep("covariate_lifestyle", 3),
               rep("covariate_genetic", 10), "outcome_primary", rep("feature_sensory", 3)))
L0[, `:=`(layer = "L0", source_analysis = "covariate", original_name = feature_id,
          raw_file = RAW_PHENOTYPE, raw_column = feature_id, selection_reason = "forced",
          pval = NA_real_, direction = NA_character_, female_only = FALSE)]

# APOE e4 dose (genome-wide significant locus)
apoe <- data.table(
  feature_id = "APOE_e4_dose", layer = "L1", source_analysis = "GWAS",
  original_name = "APOE_e4_dose", raw_file = RAW_APOE, raw_column = "APOE_e4_dose",
  selection_reason = "GWS_anchor", pval = 2.71e-08, direction = "+",
  data_requirement = "C", female_only = FALSE)

# Baseline questionnaire and biochemistry exposures (FDR within domain)
exa <- rbind(fread(file.path(IN_RESULTS, "exwas_baseline_primary.csv"))[, female_only := FALSE],
             fread(file.path(IN_RESULTS, "exwas_baseline_female_primary.csv"))[, female_only := TRUE],
             fill = TRUE)
setorder(exa, variable)  # variable-name order sets the column order of the modeling matrices
exa_sel <- exa[significant == TRUE & converged == TRUE, .(
  feature_id = variable, layer = "L1", source_analysis = paste0("ExWAS_A_", fdr_group),
  original_name = variable, raw_file = RAW_EXWAS_A, raw_column = variable,
  selection_reason = "FDR_within_domain", pval = p_value, direction = fifelse(beta > 0, "+", "-"),
  data_requirement = fifelse(female_only, "A_female", "A"), domain = fdr_group,
  var_type = var_type, female_only = female_only)]
# Reproductive anchors are kept for women even when not FDR-selected
to_force <- setdiff(FEMALE_ANCHORS, exa_sel[female_only == TRUE, feature_id])
exa_forced <- exa[variable %in% to_force, .(
  feature_id = variable, layer = "L1", source_analysis = "ExWAS_A_Reproductive (female-only)",
  original_name = variable, raw_file = RAW_EXWAS_A, raw_column = variable,
  selection_reason = "forced_female_anchor", pval = p_value, direction = fifelse(beta > 0, "+", "-"),
  data_requirement = "A_female", domain = "Reproductive (female-only)", var_type = var_type,
  female_only = TRUE)]
exa_sel <- unique(rbind(exa_sel, exa_forced, fill = TRUE), by = "feature_id")

# Follow-up questionnaire exposures (FDR within source)
exb_sel <- fread(file.path(IN_RESULTS, "exwas_followup_primary.csv"))[
  significant == TRUE & converged == TRUE, .(
  feature_id = variable, layer = "L1", source_analysis = paste0("ExWAS_B_", fdr_group),
  original_name = variable, raw_file = RAW_EXWAS_B, raw_column = variable,
  selection_reason = "FDR_within_domain", pval = p_value, direction = fifelse(beta > 0, "+", "-"),
  data_requirement = "A", domain = fdr_group, var_type = var_type, female_only = FALSE)]

# Prior diagnoses: three-character ICD-10 codes from the disease-wide scan
dwas_sel <- fread(file.path(IN_RESULTS, "dwas_icd10_primary.csv"))[
  !is.na(pval_fdr) & pval_fdr < FDR_THRESHOLD & converged == TRUE, .(
  feature_id = paste0("icd_", icd10), layer = "L1", source_analysis = "ExWAS_C",
  original_name = icd10, raw_file = RAW_DIAGNOSIS, raw_column = icd10,
  selection_reason = "FDR", pval = pval, direction = fifelse(or > 1, "+", "-"),
  data_requirement = "B", icd10 = icd10, description = description, or_value = or,
  n_exposed = n_exposed, female_only = FALSE)]

# Metabolic biomarkers: FDR-significant cluster representatives (certified measures first)
mwas <- fread(file.path(IN_RESULTS, "mwas_cluster_representatives.csv"))
ann  <- unique(fread("input/reference/nmr_metabolite_annotation.csv")[
  , .(metabolite, Display_Category, Nightingale_Group)], by = "metabolite")
mwas <- merge(mwas, ann, by.x = "protein", by.y = "metabolite", all.x = TRUE)
mwas_sel <- mwas[, .(
  feature_id = protein, layer = "L4", source_analysis = "MWAS", original_name = protein,
  raw_file = RAW_METABOLOMICS, raw_column = protein,
  selection_reason = fifelse(is_ce, "CE_sig", "nonCE_sig"), pval = pval,
  pval_fdr = pval_fdr_global, direction = fifelse(log2or > 0, "+", "-"),
  data_requirement = "D", cluster_id = cluster_id, is_ce = is_ce, sig_fdr = sig_fdr,
  display_category = Display_Category, nightingale_group = Nightingale_Group,
  or = or, or_lower = or_lower, or_upper = or_upper)]
mwas_sel <- mwas_sel[order(selection_reason != "CE_sig", pval)]

# Plasma proteins
pwas_sel <- fread(file.path(IN_RESULTS, "pwas_primary.csv"))[
  !is.na(pval_fdr) & pval_fdr < FDR_THRESHOLD & converged == TRUE, .(
  feature_id = protein, layer = "L3", source_analysis = "PWAS", original_name = protein,
  raw_file = RAW_PROTEOMICS, raw_column = protein, selection_reason = "FDR", pval = pval,
  direction = fifelse(or > 1, "+", "-"), data_requirement = "D")]

manifest <- rbindlist(list(L0, apoe, exa_sel, exb_sel, dwas_sel, mwas_sel, pwas_sel), fill = TRUE)
# A feature selected by more than one scan keeps the row with the smallest P
if (anyDuplicated(manifest$feature_id))
  manifest <- manifest[order(feature_id, pval, na.last = TRUE)][, .SD[1], by = feature_id]

fwrite(manifest, file.path(OUT_DIR, "master_feature_manifest.csv"))
