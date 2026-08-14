#!/usr/bin/env Rscript
# =============================================================================
# 07_hla_persubtype.R  (per-subtype effect at the chromosome 6 lead variant)
# =============================================================================
#
# Each subtype's cases against the SHARED control set, at the HLA class III lead
# variant, plus a heterogeneity test across the four subtype effects. Figure 5D
# draws the four odds ratios.
#
#   is_subtype_case ~ dosage + age_baseline + sex + PC1..PC10
#
# Shared controls, not one-versus-rest among cases: the question is whether the
# variant's effect on being a case differs by which subtype the case belongs to,
# so all four models are fitted against the same controls and the four betas are
# on a common footing. Comparing subtypes to each other instead would confound
# the contrast with whatever else separates the two subtypes.
#
# The variant is given, not discovered here. rs2071293, chr6:32,062,687 (GRCh37),
# effect allele A, was the lead of the subtype-stratified genome-wide scan and
# the ASSET subset-based meta-analysis (joint P = 3.0e-7); those run outside this
# package. What is reproduced here is the panel's four estimates at that locus.
#
# Heterogeneity is Cochran's Q across the four betas, inverse-variance weighted.
# Read it as descriptive: the subtypes were defined on phenotype, so a difference
# in genetic effect between them is a statement about the partition, not an
# independent test of it. A sex-stratified pass is also written, because the
# subtype labels are markedly less stable in males.
#
# CODE_DIR (env, default = current dir) must hold _subtype_map.R, captured before
# setwd(PROJECT_DIR).
#
# Input:  $COHORT_DIR/genotype/HLA_lead/hla_lead_genotypes.raw   (06, PLINK2)
#         output/subtyping/clusters/cluster_assignments_g1_{m,f}_k4.rds
#         output/ml_ready/group1_full.rds   (outcome, age_baseline, sex, PCs)
# Output: output/subtyping/evidence_genetic/hla_persubtype_effects.csv
# =============================================================================

CODE_DIR <- Sys.getenv("CODE_DIR", unset = getwd())
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
suppressPackageStartupMessages(library(data.table))

source(file.path(CODE_DIR, "_subtype_map.R"))   # SUBTYPE_MAP

COHORT_DIR <- Sys.getenv("COHORT_DIR", unset = "gwas/cohort_primary")
RAW <- file.path(COHORT_DIR, "genotype", "HLA_lead", "hla_lead_genotypes.raw")
OUT_DIR <- "output/subtyping/evidence_genetic"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

EFFECT_ALLELE <- "A"      # the allele the reported odds ratio is per-copy of
N_PC <- 10L
PC <- sprintf("PC%d", seq_len(N_PC))

# ---- [1] dosages ------------------------------------------------------------
if (!file.exists(RAW))
  stop("dosage file not found -- run 06_hla_variant_extract.sh first: ", RAW)
raw <- fread(RAW)
dose_cols <- setdiff(names(raw), c("FID","IID","PAT","MAT","SEX","PHENOTYPE"))
if (!length(dose_cols)) stop("no dosage columns in ", RAW)

# PLINK2 --export A names each column <ID>_<counted allele>. Reporting an odds
# ratio per copy of the wrong allele would invert the panel without any error, so
# the counted allele is checked rather than assumed.
counted <- sub("^.*_([ACGT]+)$", "\\1", dose_cols)
if (!all(counted == EFFECT_ALLELE))
  stop("counted allele is ", paste(unique(counted), collapse = "/"),
       " but the reported effect allele is ", EFFECT_ALLELE,
       " -- re-export with the intended allele counted, or flip deliberately")
cat("dosage columns:", paste(dose_cols, collapse = ", "), "\n")

# ---- [2] outcome, covariates, subtype labels --------------------------------
ml_cols <- c("eid", "taste_2w_strict", "age_baseline", "sex", PC)
ml <- as.data.table(readRDS("output/ml_ready/group1_full.rds"))
miss <- setdiff(ml_cols, names(ml))
if (length(miss)) stop("group1_full.rds is missing: ", paste(miss, collapse = ", "))
ml <- ml[, ..ml_cols][!is.na(taste_2w_strict)]
cat("participants with a defined outcome:", nrow(ml),
    "(cases", sum(ml$taste_2w_strict == 1), ")\n")

labels_of <- function(sx) {
  f <- sprintf("output/subtyping/clusters/cluster_assignments_g1_%s_k4.rds", sx)
  if (!file.exists(f)) stop("cluster assignments not found: ", f)
  fl <- readRDS(f)$final_labels
  data.table(eid = as.integer(names(fl)),
             subtype = unname(SUBTYPE_MAP[[sx]][as.character(as.integer(fl))]))
}
sub <- rbind(labels_of("m"), labels_of("f"))
if (anyNA(sub$subtype)) stop("unmapped cluster id -- check _subtype_map.R")

d <- merge(ml, sub, by = "eid", all.x = TRUE)
# Controls are outcome-negative participants; cases without a subtype label (not
# carried into the factor model) are dropped rather than pooled into controls.
d[, group := fifelse(taste_2w_strict == 0, "control", subtype)]
n0 <- nrow(d); d <- d[!is.na(group)]
cat("dropped unlabelled cases:", n0 - nrow(d), "\n")

n1 <- nrow(d)
d <- merge(d, raw[, c("IID", dose_cols), with = FALSE], by.x = "eid", by.y = "IID")
cat("after genotype join:", n1, "->", nrow(d), "\n")
if (nrow(d) == 0) stop("no participants have both a label and a dosage")
print(d[, .N, by = group][order(group)])

# ---- [3] one model per subtype ---------------------------------------------
fit_one <- function(dt, snp, s, use_sex) {
  x <- dt[group %in% c("control", s)]
  x[, y01 := as.integer(group == s)]
  covs <- c("age_baseline", PC, if (use_sex) "sex")
  fit <- tryCatch(glm(reformulate(c(sprintf("`%s`", snp), covs), "y01"),
                      data = x, family = binomial()),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  co <- summary(fit)$coefficients
  rn <- grep(snp, rownames(co), fixed = TRUE, value = TRUE)[1]
  if (is.na(rn)) return(NULL)
  b <- co[rn, 1]; se <- co[rn, 2]
  data.table(snp = sub("_[ACGT]+$", "", snp), subtype = s,
             n_case = sum(x$y01), n_control = sum(x$y01 == 0),
             OR = round(exp(b), 3),
             lo = round(exp(b - 1.96 * se), 3), hi = round(exp(b + 1.96 * se), 3),
             p = signif(co[rn, 4], 3), beta = b, se = se)
}

cochran_q <- function(tab) {
  w <- 1 / tab$se^2
  bbar <- sum(w * tab$beta) / sum(w)
  Q <- sum(w * (tab$beta - bbar)^2)
  pchisq(Q, df = nrow(tab) - 1, lower.tail = FALSE)
}

# sex is a covariate in the pooled pass and a stratifier in the other two, so it
# must not also be a covariate inside a single-sex stratum.
out <- list()
for (snp in dose_cols) {
  for (stratum in c("all", "female", "male")) {
    dd <- switch(stratum,
                 all    = d,
                 female = d[sex == 0],
                 male   = d[sex == 1])
    if (nrow(dd) == 0) { cat("[skip]", stratum, "- no participants\n"); next }
    rows <- rbindlist(lapply(c("A","B","C","D"), function(s)
      fit_one(dd, snp, s, use_sex = (stratum == "all"))))
    if (nrow(rows) == 4L) rows[, heter_Q_p := signif(cochran_q(rows), 3)]
    rows[, stratum := stratum]
    out[[paste(snp, stratum)]] <- rows
  }
}
res <- rbindlist(out, fill = TRUE)
if (!nrow(res)) stop("no model converged")
res[, c("beta", "se") := NULL]
setcolorder(res, c("snp", "stratum", "subtype", "n_case", "n_control",
                   "OR", "lo", "hi", "p", "heter_Q_p"))

cat("\n=== pooled (sex as covariate) ===\n"); print(res[stratum == "all"])
cat("\n=== sex-stratified ===\n"); print(res[stratum != "all"])

fwrite(res, file.path(OUT_DIR, "hla_persubtype_effects.csv"))
# Export and dx upload to RAP  (Figure 5D reads the pooled rows)
cat("\nwrote", file.path(OUT_DIR, "hla_persubtype_effects.csv"), "\n")
