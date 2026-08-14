#!/usr/bin/env Rscript
# =============================================================================
# Baseline exposure-wide association scan — batch runner
# =============================================================================
#
# Three cohorts x seven variants = 21 runs.
#
#   primary         the six-variable covariate set
#   bmi_sens        primary + BMI. The anthropometric family is excluded
#                   from this variant: adjusting body-size exposures for BMI
#                   would be adjusting them for themselves
#   smell_sens      primary + current smell change; bounds how much of the
#                   signal is chemosensory rather than taste-specific
#   minimal         age and sex only; the upper bound on effect magnitude
#   cell5 / cell15  the same primary model at a permissive and a strict
#                   sparse-cell threshold, so no reported hit depends on where
#                   that threshold was drawn
#   female_primary  the three reproductive variables, refitted on female
#                   participants without a sex term and given their own FDR
#                   family, since testing them in a mixed-sex sample is
#                   meaningless
#
# The ancestrally heterogeneous cohort additionally adjusts for ethnicity in
# every variant.
#
# FDR is Benjamini-Hochberg within each of eleven pre-specified exposure
# families, never pooled across them: the families differ enormously in size,
# and pooling would let one large family set the threshold for all the others.
#   Demographics & SES                 Lifestyle & sleep
#   Dietary intake & preferences       Anthropometric & physiological
#   Oral health                        General & mental health
#   Sensory function & pain            Clinical screening & treatment
#   Blood biochemistry & haematology
#   Environmental exposures            Handedness & laterality
#
# Usage (locally, or as a RAP job):
#   CODE_DIR=. PROJECT_DIR=/path/to/project \
#     nohup Rscript 02_run_exwas.R > baseline_exwas_run.log 2>&1 &
# =============================================================================

source(file.path(Sys.getenv("CODE_DIR", unset = "."), "exwas_common.R"))

EXPOSURE_DIR <- Sys.getenv("EXPOSURE_DIR", unset = READY_DIR)
OMICS_FILES  <- setNames(
  file.path(EXPOSURE_DIR, sprintf("exwas_baseline_%s.csv", c("group1", "group2", "group3"))),
  c("1", "2", "3"))
DICT <- file.path(EXPOSURE_DIR, "baseline_exwas_variable_dictionary.csv")

RESULTS_DIR <- file.path(OUTPUT_DIR, "baseline_exwas")
LOG_DIR     <- file.path(RESULTS_DIR, "logs")
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

METADATA_MAIN   <- file.path(RESULTS_DIR, "baseline_exwas_metadata_main.csv")
METADATA_BMI    <- file.path(RESULTS_DIR, "baseline_exwas_metadata_bmi_sens.csv")
METADATA_FEMALE <- file.path(RESULTS_DIR, "baseline_exwas_metadata_female.csv")

ts     <- format(Sys.time(), "%Y%m%d_%H%M%S")
REPORT <- file.path(RESULTS_DIR, sprintf("baseline_exwas_report_%s.txt", ts))

# -----------------------------------------------------------------------------
# Step 1: metadata files and the female phenotype subsets
# -----------------------------------------------------------------------------
# The engine needs variable, var_type and fdr_group. Three metadata files are
# built from one dictionary because three run families need different variable
# sets: the main sweep, the BMI variant with the anthropometric family removed,
# and the female-only sub-batch collapsed into its own FDR family.

cat("[PREP] building metadata from", DICT, "\n")
d <- fread(DICT)
req <- c("var_name", "var_type", "source", "female_only")
miss <- setdiff(req, names(d))
if (length(miss) > 0) stop("dictionary missing columns: ", paste(miss, collapse = ","))

d[, female_only := as.logical(female_only)]
d[is.na(female_only), female_only := FALSE]

meta_main   <- d[female_only == FALSE, .(variable = var_name, var_type, fdr_group = source)]
meta_bmi    <- meta_main[fdr_group != "Anthropometric & physiological"]
meta_female <- d[female_only == TRUE, .(variable = var_name, var_type, fdr_group = "Reproductive (female-only)")]

fwrite(meta_main,   METADATA_MAIN)
fwrite(meta_bmi,    METADATA_BMI)
fwrite(meta_female, METADATA_FEMALE)
cat(sprintf("[PREP] main:   %d variables across %d families (%s)\n",
            nrow(meta_main), uniqueN(meta_main$fdr_group),
            paste(sort(unique(meta_main$fdr_group)), collapse = ",")))
cat(sprintf("[PREP] bmi_sens: %d variables (anthropometric family excluded)\n", nrow(meta_bmi)))
cat(sprintf("[PREP] female: %d variables [Reproductive (female-only)]\n", nrow(meta_female)))

PHENO_FEMALE <- setNames(
  file.path(RESULTS_DIR, sprintf("phenotype_group%s_female.csv", 1:3)), c("1", "2", "3"))
for (g in c("1", "2", "3")) prep_sex_subset(PHENO_FILES[[g]], PHENO_FEMALE[[g]])

# -----------------------------------------------------------------------------
# Step 2: report header
# -----------------------------------------------------------------------------
writeLines(c(
  strrep("=", 62),
  " Baseline ExWAS - run report",
  strrep("=", 62),
  paste("Timestamp:      ", ts),
  paste("Outcome:        ", OUTCOME),
  paste("Method:         ", METHOD_DEFAULT),
  paste("Minimum cases:  ", MIN_CASES_DEFAULT),
  paste("Results:        ", RESULTS_DIR),
  "",
  "Covariate sets (mixed-sex variants):",
  paste("  primary:      ", COVARIATES_BASE),
  paste("  bmi_sens:     ", COVARIATES_BMI, " [anthropometric family excluded]"),
  paste("  smell_sens:   ", COVARIATES_SMELL),
  paste("  minimal:      ", COVARIATES_MINIMAL),
  "",
  "Female sub-batch (reproductive, female-only):",
  paste("  female_primary:", COVARIATES_FEMALE, sprintf("[sex == %d]", SEX_FEMALE_CODE)),
  "",
  paste("Group 3 additionally adjusts for:", GROUP3_ETHNICITY_VAR),
  paste("Factor variables:", FACTOR_VARS),
  ""
), REPORT)

# -----------------------------------------------------------------------------
# Step 3: run matrix
# -----------------------------------------------------------------------------
VARIANTS <- data.table(
  variant    = c("primary", "bmi_sens", "smell_sens", "minimal", "cell5", "cell15", "female_primary"),
  covariates = c(COVARIATES_BASE, COVARIATES_BMI, COVARIATES_SMELL, COVARIATES_MINIMAL,
                 COVARIATES_BASE, COVARIATES_BASE, COVARIATES_FEMALE),
  metadata   = c(METADATA_MAIN, METADATA_BMI, METADATA_MAIN, METADATA_MAIN,
                 METADATA_MAIN, METADATA_MAIN, METADATA_FEMALE),
  cell_drop  = c(CELL_DROP_DEFAULT, CELL_DROP_DEFAULT, CELL_DROP_DEFAULT,
                 CELL_DROP_DEFAULT, 5, 15, CELL_DROP_DEFAULT),
  cell_firth = c(CELL_FIRTH_DEFAULT, CELL_FIRTH_DEFAULT, CELL_FIRTH_DEFAULT,
                 CELL_FIRTH_DEFAULT, 10, 30, CELL_FIRTH_DEFAULT),
  female     = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)
)

for (g in c("1", "2", "3")) {
  for (i in seq_len(nrow(VARIANTS))) {
    v      <- VARIANTS[i]
    covars <- add_ethnicity(v$covariates, g)
    pheno  <- if (v$female) PHENO_FEMALE[[g]] else PHENO_FILES[[g]]
    out    <- file.path(RESULTS_DIR,
                        sprintf("baseline_exwas_results_group%s_%s.csv", g, v$variant))
    log    <- file.path(LOG_DIR, sprintf("group%s_%s_%s.log", g, v$variant, ts))
    label  <- sprintf("Group %s / %s", g, v$variant)

    cat(strrep("=", 62), "\n>>> ", label, "\n", sep = "")
    cat("    metadata  : ", basename(v$metadata), "\n", sep = "")
    cat("    cells     : ", v$cell_drop, "/", v$cell_firth, "\n", sep = "")
    cat("    phenotype : ", basename(pheno), "\n", sep = "")
    cat("    covariates: ", covars, "\n", strrep("=", 62), "\n", sep = "")

    log_lines <- call_exwas(OMICS_FILES[[g]], pheno, covars, v$metadata, out, log,
                            cell_drop = v$cell_drop, cell_firth = v$cell_firth)
    append_summary(log_lines, label, REPORT)
  }
}
# Export and dx upload to RAP  (the 21 result CSVs, the run report and the logs)

cat("\n", strrep("=", 62), "\n Baseline ExWAS complete\n", strrep("=", 62), "\n", sep = "")
cat(" runs   : 3 cohorts x ", nrow(VARIANTS), " variants = ", 3 * nrow(VARIANTS), "\n", sep = "")
cat(" results: ", RESULTS_DIR, "\n", sep = "")
cat(" report : ", REPORT, "\n", sep = "")
