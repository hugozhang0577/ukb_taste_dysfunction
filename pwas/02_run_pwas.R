# Proteome-wide association scan (PWAS) in the discovery cohort: Firth logistic regression
# of taste dysfunction on each QC'd Olink protein, primary covariate set (17 terms).
CODE_DIR <- getwd()
setwd(Sys.getenv("PROJECT_DIR", "."))

N_JOBS <- 16
dir.create("output/pwas/results", recursive = TRUE, showWarnings = FALSE)

COVARIATES  <- paste0("age_baseline,sex,olink_batch_number,", paste0("PC", 1:10, collapse = ","),
                      ",smoking,drink,surg_taste_affecting_full,townsend")
FACTOR_VARS <- "sex,olink_batch_number,smoking,drink"

system2("Rscript", c(
  shQuote(file.path(CODE_DIR, "pwas", "pwas_firth_analysis.R")),
  "--olink",       "input/analysis_ready/proteomics_group1.csv",
  "--phenotype",   "input/analysis_ready/phenotype_group1.csv,input/analysis_ready/extra_covariates.csv",
  "--outcome",     "taste_2w_strict",
  "--covariates",  COVARIATES,
  "--factor-vars", FACTOR_VARS,
  "--output",      "output/pwas/results/primary.csv",
  "--n-jobs",      N_JOBS
))
# primary.csv is copied to input/assoc_results/pwas_primary.csv
