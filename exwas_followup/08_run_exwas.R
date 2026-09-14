#!/usr/bin/env Rscript
#> ship output/followup_exwas/followup_exwas_results_group1_primary.csv -> input/assoc_results/exwas_followup_primary.csv
CODE_DIR <- Sys.getenv("CODE_DIR", unset = ".")
source(file.path(CODE_DIR, "..", "exwas_baseline", "exwas_common.R"))

EXPOSURE_DIR <- Sys.getenv("EXPOSURE_DIR", unset = READY_DIR)
OMICS_FILES  <- setNames(
  file.path(EXPOSURE_DIR, sprintf("exwas_followup_%s.csv", c("group1", "group2", "group3"))),
  c("1", "2", "3"))
DICT <- file.path(EXPOSURE_DIR, "followup_exwas_variable_dictionary.csv")

FAMILIES <- c(
  "Food preferences", "Diet (24-hour recall)", "Digestive health",
  "Mental health", "Experience of pain", "Cognitive function",
  "Work environment")

RESULTS_DIR <- file.path(OUTPUT_DIR, "followup_exwas")
LOG_DIR     <- file.path(RESULTS_DIR, "logs")
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

METADATA <- file.path(RESULTS_DIR, "followup_exwas_metadata.csv")

ts     <- format(Sys.time(), "%Y%m%d_%H%M%S")
REPORT <- file.path(RESULTS_DIR, sprintf("followup_exwas_report_%s.txt", ts))

# Step 1: metadata

cat("[PREP] building metadata from", DICT, "\n")
d <- fread(DICT)
if (!"exclude_primary" %in% names(d)) d[, exclude_primary := FALSE]

n_total    <- nrow(d)
n_excluded <- sum(d$exclude_primary, na.rm = TRUE)
d_keep     <- d[exclude_primary == FALSE | is.na(exclude_primary)]

unknown_family <- setdiff(unique(d_keep$domain), FAMILIES)
if (length(unknown_family) > 0)
  stop("dictionary carries exposure families this runner does not know: ",
       paste(unknown_family, collapse = " | "))

meta <- data.table(variable  = d_keep$var_name,
                   var_type  = d_keep$var_type,
                   fdr_group = d_keep$domain)
fwrite(meta, METADATA)

cat(sprintf("[PREP] dictionary: %d total, %d excluded, %d kept\n",
            n_total, n_excluded, nrow(meta)))
cat(sprintf("[PREP] FDR families: %s\n",
            paste(sort(unique(meta$fdr_group)), collapse = ",")))
if (n_excluded > 0) {
  cat("[PREP] excluded variables:\n")
  print(d[exclude_primary == TRUE, .(var_name, exclude_reason)])
}

# Step 2: report header
writeLines(c(
  strrep("=", 62),
  " Follow-up questionnaire ExWAS - run report",
  strrep("=", 62),
  paste("Timestamp:      ", ts),
  paste("Outcome:        ", OUTCOME),
  paste("Method:         ", METHOD_DEFAULT),
  paste("Minimum cases:  ", MIN_CASES_DEFAULT),
  paste("Metadata:       ", METADATA),
  paste("Results:        ", RESULTS_DIR),
  "",
  "Covariate sets:",
  paste("  primary:      ", COVARIATES_BASE),
  paste("  bmi_sens:     ", COVARIATES_BMI),
  paste("  smell_sens:   ", COVARIATES_SMELL),
  paste("  minimal:      ", COVARIATES_MINIMAL),
  "",
  paste("Group 3 additionally adjusts for:", GROUP3_ETHNICITY_VAR),
  paste("Factor variables:", FACTOR_VARS),
  ""
), REPORT)

# Step 3: run matrix
VARIANTS <- data.table(
  variant    = c("primary", "bmi_sens", "smell_sens", "minimal", "cell5", "cell15"),
  covariates = c(COVARIATES_BASE, COVARIATES_BMI, COVARIATES_SMELL,
                 COVARIATES_MINIMAL, COVARIATES_BASE, COVARIATES_BASE),
  cell_drop  = c(CELL_DROP_DEFAULT, CELL_DROP_DEFAULT, CELL_DROP_DEFAULT,
                 CELL_DROP_DEFAULT, 5, 15),
  cell_firth = c(CELL_FIRTH_DEFAULT, CELL_FIRTH_DEFAULT, CELL_FIRTH_DEFAULT,
                 CELL_FIRTH_DEFAULT, 10, 30)
)

for (g in c("1", "2", "3")) {
  for (i in seq_len(nrow(VARIANTS))) {
    v      <- VARIANTS[i]
    covars <- add_ethnicity(v$covariates, g)
    out    <- file.path(RESULTS_DIR,
                        sprintf("followup_exwas_results_group%s_%s.csv", g, v$variant))
    log    <- file.path(LOG_DIR, sprintf("group%s_%s_%s.log", g, v$variant, ts))
    label  <- sprintf("Group %s / %s", g, v$variant)

    cat(strrep("=", 62), "\n>>> ", label, "\n", sep = "")
    cat("    cells     : ", v$cell_drop, "/", v$cell_firth, "\n", sep = "")
    cat("    covariates: ", covars, "\n", strrep("=", 62), "\n", sep = "")

    log_lines <- call_exwas(OMICS_FILES[[g]], PHENO_FILES[[g]], covars, METADATA, out, log,
                            cell_drop = v$cell_drop, cell_firth = v$cell_firth)
    append_summary(log_lines, label, REPORT)
  }
}
# Export and dx upload to RAP  (the 18 result CSVs, the run report and the logs)

cat("\n", strrep("=", 62), "\n Follow-up ExWAS complete\n", strrep("=", 62), "\n", sep = "")
cat(" runs   : 3 cohorts x ", nrow(VARIANTS), " variants = ", 3 * nrow(VARIANTS), "\n", sep = "")
cat(" results: ", RESULTS_DIR, "\n", sep = "")
cat(" report : ", REPORT, "\n", sep = "")
