#!/usr/bin/env Rscript
# =============================================================================
# PRS — association with the outcome, in the held-out subset only
# =============================================================================
#
# Regresses the binary taste outcome on the z-scored polygenic score at each
# P-value threshold, and reports the odds ratio with a 95% CI, Nagelkerke's
# pseudo-R-squared and the area under the ROC curve. The reported threshold is
# the one with the smallest P.
#
# THE SUBSET IS THE ANALYSIS. The phenotype file must be the held-out subset from
# 01_split_cohort.R — the 30% that did not contribute to the discovery GWAS.
# Pointed at the full cohort instead, this script re-uses the discovery
# individuals, and the score then partly encodes their own outcomes: it returns
# an R-squared near 1 and an AUC near 1. Those numbers are not a strong result,
# they are the signature of the mistake. The script therefore refuses to run on
# a file larger than the split it expects, rather than trusting the caller.
#
# The score is unadjusted, with no covariates. That is the reported specification
# and it is the conservative direction here: adding the ancestry components and
# age would if anything raise the apparent performance of a score whose honest
# performance is the finding.
#
# Usage: Rscript 06_regression.R
# Input:  prs_combined_all_scores.txt          (05_merge_scores.R)
#         output/prs/pheno_validation_subset.txt  (01_split_cohort.R)
# Output: prs_heldout_performance.csv
# =============================================================================

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
suppressPackageStartupMessages({ library(data.table); library(pROC) })

PHENO_COL <- "pheno"
PRS_FILE   <- Sys.getenv("PRS_FILE",   "prs_combined_all_scores.txt")
PHENO_FILE <- Sys.getenv("PHENO_FILE",
                         file.path(PROJECT_DIR, "output/prs/pheno_validation_subset.txt"))
FULL_FILE  <- file.path(PROJECT_DIR, "output/gwas_sample_qc/cohort3_2w_strict_gwas_pheno.txt")

for (f in c(PRS_FILE, PHENO_FILE)) if (!file.exists(f)) stop("not found: ", f)

prs   <- fread(PRS_FILE,   colClasses = list(character = c("FID", "IID")))
pheno <- fread(PHENO_FILE, colClasses = list(character = c("FID", "IID")))

# Guard against being handed the whole cohort. The held-out subset is about 30%
# of it; anything much larger means the discovery individuals are back in.
if (file.exists(FULL_FILE)) {
  n_full <- nrow(fread(FULL_FILE, select = "IID"))
  if (nrow(pheno) > 0.45 * n_full)
    stop("the phenotype file has ", nrow(pheno), " rows against a cohort of ", n_full,
         " -- that is not the held-out subset. Evaluating here would re-use the ",
         "discovery individuals and report an inflated result.")
}

data <- merge(prs, pheno, by = c("FID", "IID"))
data <- data[!is.na(get(PHENO_COL))]
cat(sprintf("held-out target: %d (cases %d, controls %d)\n", nrow(data),
            sum(data[[PHENO_COL]] == 1), sum(data[[PHENO_COL]] == 0)))
if (nrow(data) == 0) stop("no overlap between the scores and the held-out subset")

score_cols <- setdiff(names(prs), c("FID", "IID"))
null_dev <- glm(reformulate("1", PHENO_COL), data = data, family = binomial())$deviance

rows <- rbindlist(lapply(score_cols, function(cc) {
  x <- scale(data[[cc]])[, 1]
  if (!is.finite(sd(data[[cc]])) || sd(data[[cc]]) == 0) return(NULL)
  d   <- data.frame(y = data[[PHENO_COL]], prs_z = x)
  fit <- glm(y ~ prs_z, data = d, family = binomial())
  s   <- summary(fit)$coefficients["prs_z", ]
  # Nagelkerke: Cox-Snell rescaled to a maximum of 1
  n   <- nrow(d)
  cs  <- 1 - exp((fit$deviance - null_dev) / n)
  nag <- cs / (1 - exp(-null_dev / n))
  data.table(threshold = cc,
             OR = exp(s[1]), lo = exp(s[1] - 1.96 * s[2]), hi = exp(s[1] + 1.96 * s[2]),
             P = s[4], nagelkerke_r2 = nag,
             auc = as.numeric(pROC::auc(pROC::roc(d$y, fitted(fit), quiet = TRUE))))
}))

setorder(rows, P)
print(rows)
fwrite(rows, "prs_heldout_performance.csv")
# Export and dx upload to RAP  (the held-out performance table)

best <- rows[1]
cat(sprintf("\nbest threshold: %s  OR %.3f (%.3f-%.3f)  P %.3g  R2 %.2e  AUC %.3f\n",
            best$threshold, best$OR, best$lo, best$hi, best$P,
            best$nagelkerke_r2, best$auc))
