# Genome-wide association sample QC: UK Biobank genotype QC exclusions, split by genetic ethnic
# grouping (Caucasian vs other), and SAIGE phenotype / keep / control files per outcome definition.
library(data.table)
setwd(Sys.getenv("PROJECT_DIR", "."))
OUT_DIR <- "output/gwas_sample_qc"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

phenotypes <- fread("output/gwas_phenotypes/taste_gwas_phenotypes.csv")
qcdata  <- fread("input/gwas_qc/QCDATA.CSV")
qcdata2 <- fread("input/gwas_qc/qc2.CSV")
qcdata2[, intersect(c("p22027", "p22019"), names(qcdata2)) := NULL]   # also in QCDATA

GWAS_COVARIATES <- c("age", "batch_number", "smoking", "drink", "BMI",
                     "surg_taste_affecting_full", "smell_any")
base_covar <- fread("output/base_table/base_table_covar.csv")
gwas_covariates <- merge(base_covar[, c("eid", setdiff(GWAS_COVARIATES, "smell_any")), with = FALSE],
                         fread("output/smell_features/smell_features.csv", select = c("eid", "smell_any")),
                         by = "eid", all.x = TRUE)

add_qc_flags <- function(dt) {
  dt[, `:=`(
    Rec_Exclusions          = fifelse(p22010 == "poor heterozygosity/missingness", 1L, 0L, na = 0L),
    Hetero_missing_outliers = fifelse(p22027 == "Yes", 1L, 0L, na = 0L),
    High_hetero_missing     = fifelse(p22018 == 2, 1L, 0L, na = 0L),
    Mixed_Ancestry          = fifelse(p22018 == 1, 1L, 0L, na = 0L),
    No_gene_participant     = fifelse(!(p22001 %in% c("Female", "Male")), 1L, 0L, na = 1L),
    Sex_Chr_aneuploidy      = fifelse(p22019 == "Yes", 1L, 0L, na = 0L),
    is_Caucasian            = fifelse(p22006 == "Caucasian", 1L, 0L, na = 0L),
    sex_plink               = fcase(p22001 == "Male", 1L, p22001 == "Female", 2L, default = 0L)
  )]
  dt[, Discordant_Sex := fifelse(p31 != p22001 & No_gene_participant == 0, 1L, 0L, na = 0L)]
  dt
}

run_sample_qc <- function(phenotype_col) {
  keep_ids <- phenotypes[!is.na(get(phenotype_col)), .(eid)]
  merged <- merge(qcdata[eid %in% keep_ids$eid], phenotypes, by = "eid", all.x = TRUE, sort = FALSE)
  merged <- merge(merged, qcdata2[eid %in% merged$eid], by = "eid", all.x = TRUE, sort = FALSE)
  add_qc_flags(merged)
  clean <- merged[Rec_Exclusions == 0 & Hetero_missing_outliers == 0 & High_hetero_missing == 0 &
                    Mixed_Ancestry == 0 & No_gene_participant == 0 & Discordant_Sex == 0 &
                    Sex_Chr_aneuploidy == 0]
  setnames(clean, paste0("p22009_a", 1:40), paste0("PC", 1:40), skip_absent = TRUE)
  clean
}

# SAIGE inputs: phenotype + sex + 10 PCs + covariates, sample keep list, control list
export_gwas_files <- function(dt, phenotype_col, ancestry_label) {
  out_prefix <- file.path(OUT_DIR, paste(phenotype_col, ancestry_label, sep = "_"))
  dt <- merge(dt, gwas_covariates, by = "eid", all.x = TRUE, sort = FALSE)
  gwas_pheno <- dt[, c(list(FID = eid, IID = eid, pheno = get(phenotype_col), sex_plink = sex_plink), .SD),
                   .SDcols = c(paste0("PC", 1:10), GWAS_COVARIATES)]
  fwrite(gwas_pheno, paste0(out_prefix, "_gwas_pheno.txt"), sep = "\t")
  fwrite(dt[, .(FID = eid, IID = eid)], paste0(out_prefix, "_keep_ids.txt"), sep = "\t", col.names = FALSE)
  fwrite(dt[get(phenotype_col) == 0, .(FID = eid, IID = eid)], paste0(out_prefix, "_control_ids.txt"),
         sep = "\t", col.names = FALSE)
}

# taste_2w_strict is the primary outcome; the others are the sensitivity definitions
for (pheno_col in c("taste_basic_strict", "taste_basic_relaxed", "taste_2w_strict",
                    "taste_2w_relaxed", "taste_4w_strict")) {
  clean <- run_sample_qc(pheno_col)
  export_gwas_files(clean[is_Caucasian == 1], pheno_col, "white")
  export_gwas_files(clean[is_Caucasian == 0], pheno_col, "nonwhite")
}
