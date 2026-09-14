#!/usr/bin/env Rscript
#> in   output/gwas_sample_qc/taste_2w_strict_white_gwas_pheno.txt
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
suppressPackageStartupMessages(library(data.table))

SEED       <- 20260413L
DISCOVERY  <- 0.70
PHENO_COL  <- "pheno"

IN  <- "output/gwas_sample_qc/taste_2w_strict_white_gwas_pheno.txt"
OUT <- "output/prs"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(IN))
  stop("GWAS phenotype file not found -- run ../../preprocessing/04_gwas_sample_qc.R: ", IN)

ph <- fread(IN)
for (v in c("FID", "IID", PHENO_COL)) if (!v %in% names(ph))
  stop("phenotype file has no column '", v, "'")
ph <- ph[!is.na(get(PHENO_COL))]
cat(sprintf("cohort: %d (cases %d, controls %d)\n", nrow(ph),
            sum(ph[[PHENO_COL]] == 1), sum(ph[[PHENO_COL]] == 0)))

set.seed(SEED)
ph[, .disc := FALSE]
for (lvl in sort(unique(ph[[PHENO_COL]]))) {
  idx <- which(ph[[PHENO_COL]] == lvl)
  take <- sample(idx, size = floor(DISCOVERY * length(idx)))
  ph[take, .disc := TRUE]
}

disc <- ph[.disc == TRUE][, .disc := NULL]
heldout <- ph[.disc == FALSE][, .disc := NULL]

pct <- function(d) 100 * mean(d[[PHENO_COL]] == 1)
cat(sprintf("discovery : %6d (%d cases, %.2f%%)\n", nrow(disc),
            sum(disc[[PHENO_COL]] == 1), pct(disc)))
cat(sprintf("held-out  : %6d (%d cases, %.2f%%)\n", nrow(heldout),
            sum(heldout[[PHENO_COL]] == 1), pct(heldout)))

# No participant may appear in both, and together they must be the whole cohort.
if (length(intersect(disc$IID, heldout$IID)) > 0)
  stop("the two subsets overlap -- the held-out evaluation would be invalid")
stopifnot(nrow(disc) + nrow(heldout) == nrow(ph))

fwrite(disc,    file.path(OUT, "discovery_70_gwas_pheno.txt"),  sep = "\t")
fwrite(heldout, file.path(OUT, "pheno_validation_subset.txt"),  sep = "\t")
# Export and dx upload to RAP  (both subsets; the split must not be redrawn)

cat("\nwritten to", OUT, "\n")
cat("next: run ../saige/ with --phenoFile=discovery_70_gwas_pheno.txt\n")
