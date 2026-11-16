# Disease-wide association scan (DWAS) of taste dysfunction on prior ICD-10 diagnoses:
# all covariate models in the discovery cohort, the primary model in the held-out cohorts.
CODE_DIR <- getwd()
setwd(Sys.getenv("PROJECT_DIR", "."))

N_JOBS      <- 10
RESULTS_DIR <- "output/dwas/results"
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

MODELS <- c(
  primary         = "age_baseline,sex,assess_centre_id,townsend,smoking,drink,BMI,dx_count_pre",
  minimal         = "age_baseline,sex,dx_count_pre",
  no_bmi          = "age_baseline,sex,assess_centre_id,townsend,smoking,drink,dx_count_pre",
  no_centre       = "age_baseline,sex,townsend,smoking,drink,BMI",   # also without the diagnosis count
  `primary+smell` = "age_baseline,sex,assess_centre_id,townsend,smoking,drink,BMI,dx_count_pre,smell_any"
)

# Inputs from 02-03, placed in input/analysis_ready/
run_dwas <- function(group, covariates, min_exposed, output) {
  system2("Rscript", c(
    shQuote(file.path(CODE_DIR, "dwas", "dwas_regression.R")),
    "--matrix",      sprintf("input/analysis_ready/icd10_matrix_%s.rds", group),
    "--phenotype",   sprintf("input/analysis_ready/phenotype_%s.csv", group),
    "--dx-count",    sprintf("input/analysis_ready/dx_count_pre_%s.csv", group),
    "--dictionary",  "input/analysis_ready/icd10_dictionary.csv",
    "--outcome",     "taste_2w_strict",
    "--covariates",  covariates,
    "--factor-vars", "sex,assess_centre_id,smoking,drink,smell_any",
    "--min-exposed", min_exposed,
    "--output",      shQuote(output),
    "--n-jobs",      N_JOBS
  ))
}

# Discovery cohort: codes with >= 100 exposed participants, every model
for (m in names(MODELS))
  run_dwas("group1", MODELS[[m]], 100, file.path(RESULTS_DIR, sprintf("group1_%s.csv", m)))

# Held-out cohorts: primary covariates, codes with >= 50 exposed participants
for (g in c("group2", "group3"))
  run_dwas(g, MODELS[["primary"]], 50, file.path(RESULTS_DIR, sprintf("%s_primary.csv", g)))

system2("Rscript", c(shQuote(file.path(CODE_DIR, "dwas", "summarise_results.R")), RESULTS_DIR))
# group{1,2,3}_primary.csv are copied to input/assoc_results/ as
# dwas_icd10_primary.csv, dwas_icd10_group2_primary.csv and dwas_icd10_group3_primary.csv
