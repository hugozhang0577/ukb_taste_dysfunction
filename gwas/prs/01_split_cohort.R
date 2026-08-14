#!/usr/bin/env Rscript
# =============================================================================
# PRS — split the cohort into a discovery and a held-out subset
# =============================================================================
#
# A polygenic score has to be weighted by a GWAS that did not see the people it
# is scored in. This splits the primary GWAS cohort 70/30 once, and everything
# downstream follows from that split:
#
#   discovery_70_gwas_pheno.txt      -> the GWAS that supplies the PRS weights
#   pheno_validation_subset.txt      -> the only place the score is evaluated
#
# The discovery GWAS is not a separate piece of code. It is
# ../saige/02_run_saige_primary.sh with --phenoFile pointed at
# discovery_70_gwas_pheno.txt and a different --outputPrefix, then
# ../saige/03_merge_sumstats.sh. Run those before 03_format_base.sh.
#
# Why this matters more than it looks: scoring is cheap and the temptation is to
# evaluate on everybody. Doing so re-uses the discovery individuals and returns a
# near-perfect result that is entirely an artefact — see the warning in
# 06_regression.R, which is the step that would silently produce it.
#
# The split is stratified by outcome so that both subsets carry the same case
# fraction; at this prevalence an unstratified split would leave the held-out
# case count noticeably variable.
#
# Usage: Rscript 01_split_cohort.R
# Input:  output/gwas_sample_qc/cohort3_2w_strict_gwas_pheno.txt
# Output: output/prs/discovery_70_gwas_pheno.txt
#         output/prs/pheno_validation_subset.txt
# =============================================================================

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
suppressPackageStartupMessages(library(data.table))

SEED       <- 20260413L
DISCOVERY  <- 0.70
PHENO_COL  <- "pheno"

IN  <- "output/gwas_sample_qc/cohort3_2w_strict_gwas_pheno.txt"
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
