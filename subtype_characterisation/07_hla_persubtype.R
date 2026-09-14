#!/usr/bin/env Rscript
#> in   input/gwas_genotype/hla_lead_genotypes.raw
#> in   output/ml_ready/group1_full.rds
#> in   output/subtyping/clusters/cluster_assignments_g1_{}_k4.rds
CODE_DIR <- Sys.getenv("CODE_DIR", unset = getwd())
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
suppressPackageStartupMessages(library(data.table))

source(file.path(CODE_DIR, "_subtype_map.R"))   # SUBTYPE_MAP

GENO_DIR <- Sys.getenv("GENO_DIR", unset = "input/gwas_genotype")
RAW <- file.path(GENO_DIR, "hla_lead_genotypes.raw")
OUT_DIR <- "output/subtyping/evidence_genetic"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

EFFECT_ALLELE <- "A"
N_PC <- 10L
PC <- sprintf("PC%d", seq_len(N_PC))

# ---- [1] dosages ------------------------------------------------------------
if (!file.exists(RAW))
  stop("dosage file not found -- run 06_hla_variant_extract.sh first: ", RAW)
raw <- fread(RAW)
dose_cols <- setdiff(names(raw), c("FID","IID","PAT","MAT","SEX","PHENOTYPE"))
if (!length(dose_cols)) stop("no dosage columns in ", RAW)

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
# Export and dx upload to RAP  (the pooled per-subtype rows)
cat("\nwrote", file.path(OUT_DIR, "hla_persubtype_effects.csv"), "\n")
