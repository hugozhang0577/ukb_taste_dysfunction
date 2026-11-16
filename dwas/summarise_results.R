# Called by 04_run_dwas.R with the results directory: per-model summary of the DWAS
# result tables, and the table of FDR-significant codes in the primary discovery model.
library(data.table)

RESULTS_DIR <- commandArgs(trailingOnly = TRUE)[1]

files <- list.files(RESULTS_DIR, pattern = "^group[123]_.*\\.csv$", full.names = TRUE)
comparison <- rbindlist(lapply(files, function(f) {
  dt  <- fread(f)
  p   <- dt$pval[!is.na(dt$pval)]
  top <- dt[which.min(pval)]
  data.table(
    Model           = sub("\\.csv$", "", basename(f)),
    N_Codes         = nrow(dt),
    N_Converged     = sum(dt$converged, na.rm = TRUE),
    Converge_Rate   = round(mean(dt$converged, na.rm = TRUE) * 100, 1),
    N_Firth         = sum(dt$method == "firth", na.rm = TRUE),
    N_Nominal       = sum(dt$pval < 0.05, na.rm = TRUE),
    N_Suggestive    = sum(dt$pval < 0.001, na.rm = TRUE),
    N_FDR           = sum(dt$pval_fdr < 0.05, na.rm = TRUE),
    N_FDR_Positive  = sum(dt$pval_fdr < 0.05 & dt$beta > 0, na.rm = TRUE),
    N_FDR_Inverse   = sum(dt$pval_fdr < 0.05 & dt$beta < 0, na.rm = TRUE),
    N_Bonferroni    = sum(dt$pval_bonf < 0.05, na.rm = TRUE),
    Lambda          = round(median(qchisq(1 - p, df = 1)) / qchisq(0.5, df = 1), 3),
    Top_Code        = top$icd10,
    Top_Description = top$description,
    Top_OR          = round(top$or, 3),
    Top_P           = format(top$pval, scientific = TRUE, digits = 2)
  )
}), fill = TRUE)
setorder(comparison, Model)
fwrite(comparison, file.path(RESULTS_DIR, "dwas_models_comparison.csv"))

# FDR < 0.05 codes of the primary discovery model
r <- fread(file.path(RESULTS_DIR, "group1_primary.csv"))
fdr_table <- r[pval_fdr < 0.05][order(pval), .(
  ICD10       = icd10,
  Description = description,
  Chapter     = chapter,
  N_Exposed   = n_exposed,
  N_Unexposed = n_unexposed,
  OR          = sprintf("%.2f", or),
  CI_95       = sprintf("%.2f-%.2f", or_lower, or_upper),
  P_value     = formatC(pval, format = "e", digits = 2),
  FDR         = formatC(pval_fdr, format = "e", digits = 2),
  Direction   = fifelse(beta > 0, "Risk", "Protective")
)]
fwrite(fdr_table, file.path(RESULTS_DIR, "dwas_fdr_significant_table.csv"))
