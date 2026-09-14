#!/usr/bin/env Rscript
#> ship output/mwas/batch3_validation/V1_group2_primary_20cov.csv -> input/assoc_results/mwas_group2_primary.csv
#> ship output/mwas/batch3_validation/V3_group3_primary_20cov.csv -> input/assoc_results/mwas_group3_primary.csv
source(file.path(Sys.getenv("CODE_DIR", unset = "."), "mwas_config.R"))

OUTPUT_DIR <- file.path(MWAS_RESULTS_DIR, "batch3_validation")
BATCH_NAME <- "MWAS batch 3: held-out cross-population validation"

BASE_COV <- paste0("age_baseline,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,",
                   "smoking,drink,assess_centre_id,fasting_hours,townsend,",
                   "blood_hour,assess_month,surg_taste_affecting_full")

MODELS <- data.table(
  name       = c("V1_group2_primary", "V2_group2_metabolic", "V3_group3_primary"),
  covariates = c(BASE_COV, paste0(BASE_COV, ",BMI,diabetes,hypertension"), BASE_COV),
  nmr        = c(NMR_FILES[["group2"]], NMR_FILES[["group2"]], NMR_FILES[["group3"]]),
  pheno      = c(paste(c(PHENO_FILES[["group2"]], PHENO_SUPP), collapse = ","),
                 paste(c(PHENO_FILES[["group2"]], PHENO_SUPP), collapse = ","),
                 paste(c(PHENO_FILES[["group3"]], PHENO_SUPP), collapse = ",")),
  description = c(
    "V1: Group 2, primary model - transportability of the discovery result",
    "V2: Group 2, plus the metabolic block",
    "V3: Group 3, primary model - fairness audit"
  )
)

t0 <- start_batch(BATCH_NAME, OUTPUT_DIR, nrow(MODELS),
                  extra = c(paste("Group 2 NMR:", NMR_FILES[["group2"]]),
                            paste("Group 3 NMR:", NMR_FILES[["group3"]])))

for (i in seq_len(nrow(MODELS))) {
  run_model(MODELS$name[i], MODELS$covariates[i], MODELS$description[i],
            OUTPUT_DIR, MODELS$nmr[i], MODELS$pheno[i], i, nrow(MODELS))
}

results <- summarise_batch(OUTPUT_DIR, BATCH_NAME)

if (!is.null(results)) {
  cat("\n", strrep("=", 72), "\n  discovery versus held-out\n", strrep("=", 72), "\n", sep = "")
  discovery_file <- file.path(MWAS_RESULTS_DIR, "primary",
                              "batch_models_comparison.csv")
  if (file.exists(discovery_file)) {
    g1 <- fread(discovery_file)
    g1_primary <- g1[grepl("^primary", Model)][1]
    v1 <- results[grepl("^V1_", Model)]
    if (nrow(g1_primary) > 0 && nrow(v1) > 0) {
      cat(sprintf("\n  Group 1 (discovery): lambda=%.3f  FDR=%d  top=%s\n",
                  g1_primary$Lambda, g1_primary$N_FDR, g1_primary$Top_Metabolite))
      cat(sprintf("  Group 2 (held out) : lambda=%.3f  FDR=%d  top=%s\n",
                  v1$Lambda, v1$N_FDR, v1$Top_Metabolite))
      if (tolower(v1$Top_Metabolite) == tolower(g1_primary$Top_Metabolite)) {
        cat("  the strongest measure is the same in both cohorts\n")
      } else {
        cat("  the strongest measure differs; compare effect directions across the\n",
            "  full panel before drawing any conclusion from that\n", sep = "")
      }
    }
  } else {
    cat("\n  [INFO] batch 1 results not found; run it first to enable this comparison\n")
  }
}

end_batch(BATCH_NAME, t0, OUTPUT_DIR)
