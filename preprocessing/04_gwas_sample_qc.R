# GWAS sample QC and phenotype export for the UK Biobank taste study
# Interactive script for RAP or a clean local review workspace.
#
# Sample-level QC and phenotype export for the genome-wide analysis.
#
# Inputs:
#   output/gwas_phenotypes/taste_gwas_phenotypes.csv   (02_define_phenotypes.R)
#   output/smell_features/smell_features.csv           (02_define_phenotypes.R)
#   output/base_table/base_table_covar.csv             (01_build_base_table.R)
#   input/gwas_qc/QCDATA.CSV
#   input/gwas_qc/qc2.CSV
#
# The exported phenotype file carries every column the SAIGE covariate models
# name, so the genome-wide analysis has no covariate input this package does not
# produce.


# Setup ==========
# Configure a project root without exposing local absolute paths.

suppressPackageStartupMessages({
  library(data.table)
})

project_dir <- Sys.getenv("PROJECT_DIR", unset = ".")
phenotype_dir <- file.path(project_dir, "output", "gwas_phenotypes")
input_qc <- file.path(project_dir, "input", "gwas_qc")
output_dir <- file.path(project_dir, "output", "gwas_sample_qc")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# Load Inputs ==========
# Definition logic:
# Load taste phenotype definitions from the companion script and UKB genotype
# QC fields. The old code used QCDATA.CSV and qc2.CSV, then dropped duplicate
# p22027 and p22019 columns from qc2 before merging.
#
# QC field use in the original code:
#   p22010 -> recommended genomic-analysis exclusion flag
#   p22027 -> heterozygosity or missingness outlier flag
#   p22018 -> original code used 1 as mixed ancestry and 2 as high missingness
#   p22001 -> genetic sex
#   p31    -> self-reported sex
#   p22019 -> sex-chromosome aneuploidy flag
#   p22006 -> genetic ethnic grouping; "Caucasian" defines the White subset
#   p22009_a1-p22009_a40 -> genetic principal components
#
# The genotyping batch (field 22000) is not read here: it is in the participant
# field extract, so 01_build_base_table.R carries it as the `batch_number`
# factor and it arrives with the other covariates below.

phenotypes <- fread(file.path(phenotype_dir, "taste_gwas_phenotypes.csv"))
qcdata <- fread(file.path(input_qc, "QCDATA.CSV"))
qcdata2 <- fread(file.path(input_qc, "qc2.CSV"))

duplicate_qc2 <- intersect(c("p22027", "p22019"), names(qcdata2))
if (length(duplicate_qc2) > 0) {
  qcdata2[, (duplicate_qc2) := NULL]
}

# Covariates for the genome-wide models. They are the same derived columns the
# panel scans use, read from the covariate table rather than re-derived here, so
# a covariate means the same thing in every analysis.
#   age                        age at the taste questionnaire (01: 2023 - birth year)
#   batch_number               genotype measurement batch, field 22000, a factor (01)
#   smoking, drink, BMI        as derived in 01
#   surg_taste_affecting_full  cumulative taste-affecting surgery (01)
#   smell_any                  binary smell change (02)
GWAS_COVARIATES <- c("age", "batch_number", "smoking", "drink", "BMI",
                     "surg_taste_affecting_full", "smell_any")

base_covar <- fread(file.path(project_dir, "output", "base_table", "base_table_covar.csv"))
smell_features <- fread(file.path(project_dir, "output", "smell_features", "smell_features.csv"),
                        select = c("eid", "smell_any"))

covar_from_base <- setdiff(GWAS_COVARIATES, "smell_any")
missing_covar <- setdiff(covar_from_base, names(base_covar))
if (length(missing_covar) > 0)
  stop("base_table_covar.csv is missing GWAS covariate(s): ",
       paste(missing_covar, collapse = ", "))

gwas_covariates <- merge(base_covar[, c("eid", covar_from_base), with = FALSE],
                         smell_features, by = "eid", all.x = TRUE)
cat(sprintf("GWAS covariates assembled for %s participants\n",
            format(nrow(gwas_covariates), big.mark = ",")))


# QC Functions ==========
# Definition logic:
# Build the sample-level QC flags. Participants are
# removed for recommended genomic-analysis exclusion, heterozygosity/missingness
# outlier status, high missingness, mixed ancestry flag, missing genetic sex,
# discordant reported/genetic sex, or sex-chromosome aneuploidy.

rename_pcs <- function(dt) {
  old <- paste0("p22009_a", 1:40)
  new <- paste0("PC", 1:40)
  present <- old %in% names(dt)
  if (any(present)) {
    setnames(dt, old[present], new[present])
  }
  dt
}

add_qc_flags <- function(dt) {
  dt[, `:=`(
    Rec_Exclusions = fifelse(p22010 == "poor heterozygosity/missingness", 1L, 0L, na = 0L),
    Hetero_missing_outliers = fifelse(p22027 == "Yes", 1L, 0L, na = 0L),
    High_hetero_missing = fifelse(p22018 == 2, 1L, 0L, na = 0L),
    Mixed_Ancestry = fifelse(p22018 == 1, 1L, 0L, na = 0L),
    No_gene_participant = fifelse(!(p22001 %in% c("Female", "Male")), 1L, 0L, na = 1L),
    Discordant_Sex = fifelse(p31 != p22001 & No_gene_participant == 0, 1L, 0L, na = 0L),
    Sex_Chr_aneuploidy = fifelse(p22019 == "Yes", 1L, 0L, na = 0L),
    is_Caucasian = fifelse(p22006 == "Caucasian", 1L, 0L, na = 0L),
    # Genetic sex in PLINK coding (1 male, 2 female). Named apart from the
    # self-reported `sex` (p31, 1 male / 0 female) the panel scans use: different
    # source field and different coding, so they must not share a name.
    sex_plink = fcase(p22001 == "Male", 1L, p22001 == "Female", 2L, default = 0L)
  )]
  dt
}

run_sample_qc <- function(phenotype_table, phenotype_col) {
  keep_ids <- phenotype_table[!is.na(get(phenotype_col)), .(eid)]
  merged <- merge(qcdata[eid %in% keep_ids$eid],
                  phenotype_table,
                  by = "eid",
                  all.x = TRUE,
                  sort = FALSE)
  merged <- merge(merged,
                  qcdata2[eid %in% merged$eid],
                  by = "eid",
                  all.x = TRUE,
                  sort = FALSE)
  add_qc_flags(merged)

  clean <- merged[
    Rec_Exclusions == 0 &
      Hetero_missing_outliers == 0 &
      High_hetero_missing == 0 &
      Mixed_Ancestry == 0 &
      No_gene_participant == 0 &
      Discordant_Sex == 0 &
      Sex_Chr_aneuploidy == 0
  ]
  rename_pcs(clean)
  clean
}


# Define GWAS Cohorts ==========
# Definition logic:
# The four main GWAS cohorts are retained. cohort6 reproduces the
# separate 4-week-or-impact sensitivity phenotype described in the README.

cohort_definitions <- data.table(
  cohort = c("cohort1_basic_strict",
             "cohort2_basic_relaxed",
             "cohort3_2w_strict",
             "cohort4_2w_relaxed",
             "cohort6_4w_strict"),
  phenotype = c("taste_basic_strict",
                "taste_basic_relaxed",
                "taste_2w_strict",
                "taste_2w_relaxed",
                "taste_4w_strict")
)


# Run Sample QC ==========
# Definition logic:
# For each phenotype, export pre-QC counts, post-QC counts, White and non-White
# subsets, and case-control counts. The White subset follows p22006 ==
# "Caucasian", matching the original reviewed scripts.

qc_results <- list()
qc_summary <- list()

for (i in seq_len(nrow(cohort_definitions))) {
  cohort_name <- cohort_definitions$cohort[i]
  pheno_col <- cohort_definitions$phenotype[i]

  clean <- run_sample_qc(phenotypes, pheno_col)
  white <- clean[is_Caucasian == 1]
  nonwhite <- clean[is_Caucasian == 0]

  qc_results[[cohort_name]] <- list(clean = clean, white = white, nonwhite = nonwhite)

  qc_summary[[cohort_name]] <- data.table(
    cohort = cohort_name,
    phenotype = pheno_col,
    n_after_sample_qc = nrow(clean),
    n_white = nrow(white),
    n_white_control = sum(white[[pheno_col]] == 0, na.rm = TRUE),
    n_white_case = sum(white[[pheno_col]] == 1, na.rm = TRUE),
    n_nonwhite = nrow(nonwhite),
    n_nonwhite_control = sum(nonwhite[[pheno_col]] == 0, na.rm = TRUE),
    n_nonwhite_case = sum(nonwhite[[pheno_col]] == 1, na.rm = TRUE)
  )
}

qc_summary <- rbindlist(qc_summary)
fwrite(qc_summary, file.path(output_dir, "gwas_sample_qc_summary.csv"))
print(qc_summary)


# Export GWAS Phenotype Files ==========
# Definition logic:
# Export files for SAIGE/PLINK-style GWAS workflows:
#   gwas_pheno: FID, IID, phenotype, genetic sex (PLINK coding), PC1-PC10, the
#               genotyping batch, and the covariates the SAIGE models name
#   keep_ids: FID/IID sample list
#   control_ids: FID/IID controls for control-only HWE filtering

export_gwas_files <- function(dt, phenotype_col, cohort_name, ancestry_label, use_pc = 10) {
  pc_cols <- paste0("PC", seq_len(use_pc))
  missing_pc <- setdiff(pc_cols, names(dt))
  if (length(missing_pc) > 0) {
    stop("Missing PC columns for ", cohort_name, ": ", paste(missing_pc, collapse = ", "))
  }

  out_prefix <- file.path(output_dir, paste(cohort_name, ancestry_label, sep = "_"))

  # A name already present on both sides would come back as .x/.y and silently
  # change which column the export picks, so refuse rather than guess.
  clash <- intersect(GWAS_COVARIATES, names(dt))
  if (length(clash) > 0)
    stop("covariate name(s) already present in the QC table: ",
         paste(clash, collapse = ", "))

  n_before <- nrow(dt)
  dt <- merge(dt, gwas_covariates, by = "eid", all.x = TRUE, sort = FALSE)
  if (nrow(dt) != n_before)
    stop("covariate join changed the row count for ", cohort_name, "_", ancestry_label,
         " (", n_before, " -> ", nrow(dt), "); check for duplicate eids")

  missing_cov <- setdiff(GWAS_COVARIATES, names(dt))
  if (length(missing_cov) > 0)
    stop("Missing GWAS covariate(s) for ", cohort_name, ": ",
         paste(missing_cov, collapse = ", "))

  gwas_pheno <- dt[, c(
    list(FID = eid, IID = eid, pheno = get(phenotype_col), sex_plink = sex_plink),
    .SD
  ), .SDcols = c(pc_cols, GWAS_COVARIATES)]

  # Report covariate completeness rather than letting SAIGE drop the rows silently
  for (cv in GWAS_COVARIATES) {
    n_na <- sum(is.na(gwas_pheno[[cv]]))
    if (n_na > 0)
      cat(sprintf("  %s_%s: %s of %s missing %s\n", cohort_name, ancestry_label,
                  format(n_na, big.mark = ","), format(nrow(gwas_pheno), big.mark = ","), cv))
  }

  fwrite(gwas_pheno,
         paste0(out_prefix, "_gwas_pheno.txt"),
         sep = "\t")
  # Export and dx upload to RAP  (the SAIGE phenotype/covariate file; copy the
  # cohort's file into $COHORT_DIR/SAIGE/input/ before running gwas/saige/)

  keep_ids <- dt[, .(FID = eid, IID = eid)]
  fwrite(keep_ids,
         paste0(out_prefix, "_keep_ids.txt"),
         sep = "\t",
         col.names = FALSE)

  control_ids <- dt[get(phenotype_col) == 0, .(FID = eid, IID = eid)]
  fwrite(control_ids,
         paste0(out_prefix, "_control_ids.txt"),
         sep = "\t",
         col.names = FALSE)
}

for (cohort_name in names(qc_results)) {
  pheno_col <- cohort_definitions[cohort == cohort_name, phenotype]
  export_gwas_files(qc_results[[cohort_name]]$white, pheno_col, cohort_name, "white")
  export_gwas_files(qc_results[[cohort_name]]$nonwhite, pheno_col, cohort_name, "nonwhite")
}
