#!/usr/bin/env Rscript
#> in   output/gwas_sample_qc/taste_2w_strict_white_gwas_pheno.txt
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
suppressPackageStartupMessages({ library(data.table); library(pROC) })

PHENO_COL <- "pheno"
PRS_FILE   <- Sys.getenv("PRS_FILE",   "prs_combined_all_scores.txt")
PHENO_FILE <- Sys.getenv("PHENO_FILE",
                         file.path(PROJECT_DIR, "output/prs/pheno_validation_subset.txt"))
FULL_FILE  <- file.path(PROJECT_DIR, "output/gwas_sample_qc/taste_2w_strict_white_gwas_pheno.txt")

for (f in c(PRS_FILE, PHENO_FILE)) if (!file.exists(f)) stop("not found: ", f)

prs   <- fread(PRS_FILE,   colClasses = list(character = c("FID", "IID")))
pheno <- fread(PHENO_FILE, colClasses = list(character = c("FID", "IID")))

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
