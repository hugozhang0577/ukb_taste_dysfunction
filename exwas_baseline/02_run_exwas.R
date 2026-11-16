# Baseline ExWAS: one logistic model per baseline exposure in the discovery cohort (group1) and
# the two held-out cohorts (group2, group3), for the primary model and six sensitivity variants.
source("exwas_baseline/exwas_common.R")
setwd(PROJECT_DIR)

OMICS_FILES <- setNames(
  file.path(READY_DIR, sprintf("exwas_baseline_%s.csv", c("group1", "group2", "group3"))),
  c("1", "2", "3"))
DICT        <- file.path(READY_DIR, "baseline_exwas_variable_dictionary.csv")
RESULTS_DIR <- file.path(OUTPUT_DIR, "baseline_exwas")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

METADATA_MAIN   <- file.path(RESULTS_DIR, "baseline_exwas_metadata_main.csv")
METADATA_BMI    <- file.path(RESULTS_DIR, "baseline_exwas_metadata_bmi_sens.csv")
METADATA_FEMALE <- file.path(RESULTS_DIR, "baseline_exwas_metadata_female.csv")

# FDR family = exposure family; female-only exposures form one family, tested in women only
d <- fread(DICT)
d[, female_only := as.logical(female_only)]
d[is.na(female_only), female_only := FALSE]
meta_main <- d[female_only == FALSE, .(variable = var_name, var_type, fdr_group = source)]
fwrite(meta_main, METADATA_MAIN)
# BMI enters as a covariate, so the anthropometric family is not tested in that variant
fwrite(meta_main[fdr_group != "Anthropometric & physiological"], METADATA_BMI)
fwrite(d[female_only == TRUE, .(variable = var_name, var_type, fdr_group = "Reproductive (female-only)")],
       METADATA_FEMALE)

PHENO_FEMALE <- setNames(
  file.path(RESULTS_DIR, sprintf("phenotype_group%s_female.csv", 1:3)), c("1", "2", "3"))
for (g in c("1", "2", "3")) prep_sex_subset(PHENO_FILES[[g]], PHENO_FEMALE[[g]])

VARIANTS <- data.table(
  variant    = c("primary", "bmi_sens", "smell_sens", "minimal", "cell5", "cell15", "female_primary"),
  covariates = c(COVARIATES_BASE, COVARIATES_BMI, COVARIATES_SMELL, COVARIATES_MINIMAL,
                 COVARIATES_BASE, COVARIATES_BASE, COVARIATES_FEMALE),
  metadata   = c(METADATA_MAIN, METADATA_BMI, METADATA_MAIN, METADATA_MAIN,
                 METADATA_MAIN, METADATA_MAIN, METADATA_FEMALE),
  cell_drop  = c(CELL_DROP_DEFAULT, CELL_DROP_DEFAULT, CELL_DROP_DEFAULT, CELL_DROP_DEFAULT,
                 5, 15, CELL_DROP_DEFAULT),
  cell_firth = c(CELL_FIRTH_DEFAULT, CELL_FIRTH_DEFAULT, CELL_FIRTH_DEFAULT, CELL_FIRTH_DEFAULT,
                 10, 30, CELL_FIRTH_DEFAULT),
  female     = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE))

# All three cohorts; group3 models also adjust for ethnicity
for (g in c("1", "2", "3")) {
  for (i in seq_len(nrow(VARIANTS))) {
    v <- VARIANTS[i]
    call_exwas(OMICS_FILES[[g]],
               if (v$female) PHENO_FEMALE[[g]] else PHENO_FILES[[g]],
               add_ethnicity(v$covariates, g), v$metadata,
               file.path(RESULTS_DIR, sprintf("baseline_exwas_results_group%s_%s.csv", g, v$variant)),
               cell_drop = v$cell_drop, cell_firth = v$cell_firth)
  }
}
