#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
GWAS_FILE <- Sys.getenv("GWAS_FILE",
  file.path(PROJECT_DIR, "input/assoc_results/gwas_primary_sumstats.tsv.gz"))
OUT_FILE  <- Sys.getenv("OUT_FILE",
  file.path(PROJECT_DIR, "output/gwas_fuma/gwas_for_fuma.txt"))
dir.create(dirname(OUT_FILE), showWarnings = FALSE, recursive = TRUE)

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
