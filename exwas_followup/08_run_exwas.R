#!/usr/bin/env Rscript
# Follow-up questionnaire ExWAS: logistic regression of the outcome on each follow-up exposure in
# group1 (discovery) and the held-out group2 and group3 (group3 also adjusted for ethnicity), for the
# primary model and five sensitivity variants. BH FDR within each questionnaire family.
library(data.table)
source("exwas_baseline/exwas_common.R")
setwd(Sys.getenv("PROJECT_DIR", "."))

RESULTS_DIR <- "output/followup_exwas"
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
METADATA <- file.path(RESULTS_DIR, "followup_exwas_metadata.csv")

dict <- fread("input/analysis_ready/followup_exwas_variable_dictionary.csv")
fwrite(dict[exclude_primary == FALSE, .(variable = var_name, var_type, fdr_group = domain)], METADATA)

VARIANTS <- data.table(
  variant    = c("primary", "bmi_sens", "smell_sens", "minimal", "cell5", "cell15"),
  covariates = c(COVARIATES_BASE, COVARIATES_BMI, COVARIATES_SMELL, COVARIATES_MINIMAL,
                 COVARIATES_BASE, COVARIATES_BASE),
  cell_drop  = c(rep(CELL_DROP_DEFAULT, 4), 5, 15),
  cell_firth = c(rep(CELL_FIRTH_DEFAULT, 4), 10, 30))

# The primary results are copied to input/assoc_results/ as exwas_followup_primary.csv (group1)
# and exwas_followup_group2_primary.csv / exwas_followup_group3_primary.csv (held-out).
for (g in c("1", "2", "3")) {
  for (i in seq_len(nrow(VARIANTS))) {
    v   <- VARIANTS[i]
    out <- file.path(RESULTS_DIR, sprintf("followup_exwas_results_group%s_%s.csv", g, v$variant))
    call_exwas(sprintf("input/analysis_ready/exwas_followup_group%s.csv", g), PHENO_FILES[[g]],
               add_ethnicity(v$covariates, g), METADATA, out, sub("\\.csv$", ".log", out),
               cell_drop = v$cell_drop, cell_firth = v$cell_firth)
  }
}
