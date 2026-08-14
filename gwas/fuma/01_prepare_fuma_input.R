#!/usr/bin/env Rscript
# =============================================================================
# FUMA — prepare SNP2GENE input from SAIGE summary statistics
# =============================================================================
# Converts the merged SAIGE summary statistics (../saige/03_merge_sumstats.sh)
# to the column format FUMA expects. SAIGE BETA is on the Allele2 (effect) allele, so
# A1 = Allele2, A2 = Allele1. Rows with missing P/BETA/SE are dropped (FUMA
# rejects NA); scientific-notation BP is avoided by writing integer positions.
#
# Defaults to the primary scan; override GWAS_FILE/OUT_FILE for another run.
# =============================================================================
suppressPackageStartupMessages(library(data.table))

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
GWAS_FILE <- Sys.getenv("GWAS_FILE",
  file.path(PROJECT_DIR, "gwas/cohort_primary/SAIGE/step2/gwas_sumstats.tsv.gz"))
OUT_FILE  <- Sys.getenv("OUT_FILE",
  file.path(PROJECT_DIR, "gwas/cohort_primary/SAIGE/fuma/gwas_for_fuma.txt"))

gwas <- fread(GWAS_FILE)
fuma_input <- gwas[, .(
  SNP = MarkerID, CHR = CHR, BP = as.integer(POS),
  A1 = Allele2, A2 = Allele1,             # SAIGE effect allele = Allele2
  BETA = BETA, SE = SE, P = p.value,
  N = N_case + N_ctrl)]

n0 <- nrow(fuma_input)
fuma_input <- fuma_input[!is.na(P) & !is.na(BETA) & !is.na(SE)]
cat(sprintf("variants: %d -> %d after NA drop; min P = %.3e\n",
            n0, nrow(fuma_input), min(fuma_input$P, na.rm = TRUE)))

# Export and dx upload to RAP  (upload .txt.gz, then submit to FUMA SNP2GENE)
fwrite(fuma_input, OUT_FILE, sep = "\t", quote = FALSE, scipen = 999)
R.utils::gzip(OUT_FILE, overwrite = TRUE)
cat("FUMA upload file ->", paste0(OUT_FILE, ".gz"), "\n")
