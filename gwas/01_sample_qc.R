# Genotype sample QC for participants with a taste phenotype: UK Biobank central QC exclusions
# (fields 22010, 22027, 22018, 22001, 31, 22019), per-step exclusion counts, and the QC-passed ID lists
# (all ancestries, White British by field 22006, non-White).
library(data.table)
setwd(Sys.getenv("PROJECT_DIR", "."))
OUT_DIR <- "output/gwas_qc"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

qc1   <- fread("input/gwas_qc/qc_fields.csv")     # p22010, p22027, p22019, p22006, p31, p22001
qc2   <- fread("input/gwas_qc/qc_fields2.csv")    # p22018
pheno <- fread("input/gwas_qc/taste_pheno.csv")   # eid + taste phenotype

dt <- merge(pheno[, .(eid)], qc1, by = "eid", all.x = TRUE)
dt <- merge(dt, qc2[, .(eid, p22018)], by = "eid", all.x = TRUE)

dt[, Rec_Exclusions          := as.integer(p22010 == "poor heterozygosity/missingness")]
dt[, Hetero_missing_outliers := as.integer(p22027 == "Yes")]
dt[, High_hetero_missing     := as.integer(p22018 == 2)]   # relatedness, high heterozygosity/missingness
dt[, Mixed_Ancestry          := as.integer(p22018 == 1)]   # relatedness/ancestry exclusion
dt[, No_gene_participant     := as.integer(!(p22001 %in% c("Female", "Male")))]
dt[, Discordant_Sex          := as.integer(p31 != p22001 & No_gene_participant == 0)]
dt[, Sex_Chr_aneuploidy      := as.integer(p22019 == "Yes")]
dt[, is_Caucasian            := as.integer(p22006 == "Caucasian")]
flag_cols <- c("Rec_Exclusions", "Hetero_missing_outliers", "High_hetero_missing", "Mixed_Ancestry",
               "No_gene_participant", "Discordant_Sex", "Sex_Chr_aneuploidy", "is_Caucasian")
for (cc in flag_cols) dt[is.na(get(cc)), (cc) := 0L]   # QC fields are blank when not applicable
dt[, any_exclusion := as.integer(rowSums(.SD) > 0), .SDcols = setdiff(flag_cols, "is_Caucasian")]

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
  n = c(nrow(dt), dt[, sum(Rec_Exclusions)], dt[, sum(Discordant_Sex)], dt[, sum(Sex_Chr_aneuploidy)],
        dt[, sum(Hetero_missing_outliers)], dt[, sum(High_hetero_missing)], dt[, sum(Mixed_Ancestry)],
        dt[, sum(No_gene_participant)], dt[, sum(any_exclusion)], nrow(clean), nrow(white), nrow(nonwhite)))
fwrite(summ, file.path(OUT_DIR, "qc_step_summary.csv"))

for (x in list(list(clean, "qc_clean_eid.txt"), list(white, "qc_white_eid.txt"), list(nonwhite, "qc_nonwhite_eid.txt")))
  fwrite(x[[1]][, .(FID = eid, IID = eid)], file.path(OUT_DIR, x[[2]]), sep = " ", col.names = FALSE)
