#!/usr/bin/env Rscript
# =============================================================================
# GWAS — individual-level (sample) quality control
# =============================================================================
#
# Reproduces the individual-level exclusions reported under
# "Genotype data and quality control" for the taste-dysfunction GWAS.
#
# All exclusion flags are built from UK Biobank central QC fields (no genotypes
# are re-derived here); central genotyping/imputation/QC followed the UK Biobank
# standard release (Bycroft 2018). This script applies the study-specific
# individual-level exclusions and splits the QC-passing sample by genetic
# ancestry. Variant-level QC (MAF / HWE / INFO / missingness) is performed in
# PLINK 2 separately (see README) and is not part of this script.
#
# Input  (under $PROJECT_DIR/input/gwas_qc/):
#   qc_fields.csv   — eid + UK Biobank QC fields:
#                       p22010 (recommended genomic-analysis exclusions),
#                       p22027 (heterozygosity/missingness-rate outliers),
#                       p22019 (sex-chromosome aneuploidy),
#                       p22006 (genetic ethnic grouping; "Caucasian"),
#                       p31    (reported sex), p22001 (genetic sex)
#   qc_fields2.csv  — eid + p22018 (relatedness/ancestry exclusions)
#   taste_pheno.csv — eid + taste phenotype (pure-smell cases already removed)
#
# Output (under $PROJECT_DIR/output/gwas_qc/):
#   qc_sample_flags.csv  — per-participant exclusion flags
#   qc_step_summary.csv  — per-step exclusion counts (Methods source)
#   qc_clean_eid.txt     — FID/IID of all QC-passing participants
#   qc_white_eid.txt     — FID/IID of QC-passing White British (primary cohort)
#   qc_nonwhite_eid.txt  — FID/IID of QC-passing non-White (held-out validation)
#
# Convention: data.table; relative paths under $PROJECT_DIR; fail-visibly counts.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
setwd(PROJECT_DIR)

IN_DIR  <- "input/gwas_qc"
OUT_DIR <- "output/gwas_qc"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

QC_FIELDS   <- file.path(IN_DIR, "qc_fields.csv")    # p22010/p22027/p22019/p22006/p31/p22001
QC_FIELDS2  <- file.path(IN_DIR, "qc_fields2.csv")   # p22018
TASTE_PHENO <- file.path(IN_DIR, "taste_pheno.csv")  # eid + taste phenotype

# Per-step counts as reported in the Supplementary QC table. Re-running this
# script against the same UK Biobank release should reproduce them exactly; any
# difference means the input extract differs and the numbers must not be
# reported without investigating why.
REFERENCE <- list(start = 189604L, rec = 180L, sca = 149L, het = 399L,
                  relhi = 303L, relmix = 216L, nogeno = 4247L, discsex = 99L,
                  excl_union = 5509L, clean = 184095L, white = 158867L,
                  nonwhite = 25228L)

# ---- load ------------------------------------------------------------------
cat("[1] loading QC fields + taste phenotype ...\n")
stopifnot(file.exists(QC_FIELDS), file.exists(QC_FIELDS2), file.exists(TASTE_PHENO))

qc1   <- fread(QC_FIELDS)
qc2   <- fread(QC_FIELDS2)
pheno <- fread(TASTE_PHENO)

# restrict to genotyped participants who carry a (non-pure-smell) taste phenotype
n_pheno <- nrow(pheno)
dt <- merge(pheno[, .(eid)], qc1, by = "eid", all.x = TRUE)
dt <- merge(dt, qc2[, .(eid, p22018)], by = "eid", all.x = TRUE)
cat(sprintf("    taste phenotype N            : %d\n", n_pheno))
cat(sprintf("    after join to QC fields N    : %d\n", nrow(dt)))

# ---- build exclusion flags ---------------------------------------------------
cat("\n[2] building exclusion flags ...\n")
dt[, Rec_Exclusions          := as.integer(p22010 == "poor heterozygosity/missingness")]
dt[, Hetero_missing_outliers := as.integer(p22027 == "Yes")]
dt[, High_hetero_missing     := as.integer(p22018 == 2)]   # relatedness, high het/missing
dt[, Mixed_Ancestry          := as.integer(p22018 == 1)]   # relatedness/ancestry exclusion
dt[, No_gene_participant      := as.integer(!(p22001 %in% c("Female", "Male")))]
dt[, Discordant_Sex          := as.integer(p31 != p22001 & No_gene_participant == 0)]
dt[, Sex_Chr_aneuploidy      := as.integer(p22019 == "Yes")]
dt[, is_Caucasian            := as.integer(p22006 == "Caucasian")]
# UK Biobank QC fields are blank for non-applicable participants -> treat NA as 0
flag_cols <- c("Rec_Exclusions", "Hetero_missing_outliers", "High_hetero_missing",
               "Mixed_Ancestry", "No_gene_participant", "Discordant_Sex",
               "Sex_Chr_aneuploidy", "is_Caucasian")
for (cc in flag_cols) dt[is.na(get(cc)), (cc) := 0L]

# union of exclusions (a participant may meet more than one criterion)
excl_flags <- setdiff(flag_cols, "is_Caucasian")
dt[, any_exclusion := as.integer(rowSums(.SD) > 0), .SDcols = excl_flags]

# ---- per-step summary (Methods source) -------------------------------------
cat("\n[3] per-step exclusion counts ...\n")
clean    <- dt[any_exclusion == 0]
white    <- clean[is_Caucasian == 1]
nonwhite <- clean[is_Caucasian == 0]

summ <- data.table(
  step = c("Starting genotyped taste-defined sample",
           "Recommended genomic-analysis exclusions (22010)",
           "Sex discordance, reported vs genetic (31 vs 22001)",
           "Sex-chromosome aneuploidy (22019)",
           "Heterozygosity/missingness-rate outliers (22027)",
           "Relatedness exclusion, high het/missing (22018=2)",
           "Relatedness/ancestry exclusion, mixed (22018=1)",
           "No genotype data (22001 missing)",
           "Any exclusion (union)",
           "QC-pass (all ancestries)",
           "QC-pass White British (22006 Caucasian)",
           "QC-pass non-White"),
  n = c(nrow(dt),
        dt[, sum(Rec_Exclusions)],
        dt[, sum(Discordant_Sex)],
        dt[, sum(Sex_Chr_aneuploidy)],
        dt[, sum(Hetero_missing_outliers)],
        dt[, sum(High_hetero_missing)],
        dt[, sum(Mixed_Ancestry)],
        dt[, sum(No_gene_participant)],
        dt[, sum(any_exclusion)],
        nrow(clean), nrow(white), nrow(nonwhite)))
print(summ)

# ---- write outputs ---------------------------------------------------------
cat("\n[4] writing outputs ...\n")
fwrite(dt[, c("eid", flag_cols, "any_exclusion"), with = FALSE],
       file.path(OUT_DIR, "qc_sample_flags.csv"))
fwrite(summ, file.path(OUT_DIR, "qc_step_summary.csv"))
fwrite(clean[,    .(FID = eid, IID = eid)], file.path(OUT_DIR, "qc_clean_eid.txt"),
       sep = " ", col.names = FALSE)
fwrite(white[,    .(FID = eid, IID = eid)], file.path(OUT_DIR, "qc_white_eid.txt"),
       sep = " ", col.names = FALSE)
fwrite(nonwhite[, .(FID = eid, IID = eid)], file.path(OUT_DIR, "qc_nonwhite_eid.txt"),
       sep = " ", col.names = FALSE)

# ---- reproduce the reported per-step counts --------------------------------
cat("\n[5] comparing per-step counts with the reported QC table ...\n")
chk <- c(start = nrow(dt), rec = dt[, sum(Rec_Exclusions)],
         sca = dt[, sum(Sex_Chr_aneuploidy)], het = dt[, sum(Hetero_missing_outliers)],
         relhi = dt[, sum(High_hetero_missing)], relmix = dt[, sum(Mixed_Ancestry)],
         nogeno = dt[, sum(No_gene_participant)], discsex = dt[, sum(Discordant_Sex)],
         excl_union = dt[, sum(any_exclusion)], clean = nrow(clean),
         white = nrow(white), nonwhite = nrow(nonwhite))
ok <- TRUE
for (k in names(REFERENCE)) {
  match_k <- isTRUE(chk[[k]] == REFERENCE[[k]])
  if (!match_k) ok <- FALSE
  cat(sprintf("    %-12s reported %7d | this run %7d | %s\n",
              k, REFERENCE[[k]], chk[[k]], ifelse(match_k, "OK", "DIFFERS")))
}
cat(if (ok) "    All per-step counts reproduce the reported QC table.\n" else
    "    [WARN] counts differ from the reported table - check the input extract.\n")

# NOTE: the primary GWAS additionally restricts the White British QC-pass set
# (158,867) to participants with a defined taste_2w_strict primary outcome,
# giving the analytic discovery sample of 158,370 (5,753 cases / 152,617
# controls). That phenotype thresholding is applied when assembling the SAIGE
# phenotype/covariate file (see 02_*), not in this individual-level QC step.

cat("\n=== GWAS sample QC complete ===\n")
