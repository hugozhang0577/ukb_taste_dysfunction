# FUMA SNP2GENE input from the primary GWAS summary statistics (effect allele = SAIGE Allele2).
library(data.table)
setwd(Sys.getenv("PROJECT_DIR", "."))
dir.create("output/gwas_fuma", showWarnings = FALSE, recursive = TRUE)

gwas <- fread("input/assoc_results/gwas_primary_sumstats.tsv.gz")
fuma_input <- gwas[, .(SNP = MarkerID, CHR = CHR, BP = as.integer(POS), A1 = Allele2, A2 = Allele1,
                       BETA = BETA, SE = SE, P = p.value, N = N_case + N_ctrl)]
fuma_input <- fuma_input[!is.na(P) & !is.na(BETA) & !is.na(SE)]
fwrite(fuma_input, "output/gwas_fuma/gwas_for_fuma.txt.gz", sep = "\t", quote = FALSE, scipen = 999)
