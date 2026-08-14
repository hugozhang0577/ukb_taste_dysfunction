#!/usr/bin/env Rscript
# =============================================================================
# APOE diplotype — case/control association with taste dysfunction
# =============================================================================
# Infers the APOE diplotype (e2/e3/e4) from the rs429358 + rs7412 .raw dosages,
# merges the phenotype/covariate file, and tests three encodings against the
# primary taste outcome with the primary covariate set:
#   (1) e4 carrier (binary), (2) per-e4-allele dose (0/1/2 trend),
#   (3) 6-level diplotype (e3/e3 reference).
# Exports the forest-plot table + genotype-frequency tables (Fig 3d data).
#
# ALT-allele assumption, as used for the reported
# OR): in apoe_genotypes.raw, rs429358 is coded to the T allele and rs7412 to
# the C allele, so e4_count = 2 - rs429358_dose and e2_count = 2 - rs7412_dose.
# Verify against the .raw header (column suffix = ALT allele) before reuse;
# 01_extract_apoe_genotype.sh performs the same detection when exporting the
# multi-ancestry genotype calls.
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(broom)
})

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
APOE_DIR <- file.path(PROJECT_DIR, "gwas/cohort_primary/genotype/APOE_type")
PHENO_FILE <- Sys.getenv("PHENO_FILE", file.path(APOE_DIR, "gwas_pheno.txt"))
RES_DIR <- file.path(APOE_DIR, "Results")
dir.create(RES_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- genotype -> diplotype -------------------------------------------------
geno <- fread(file.path(APOE_DIR, "apoe_genotypes.raw"))
setnames(geno,
         grep("rs429358", names(geno), value = TRUE)[1], "rs429358_dose")
setnames(geno,
         grep("rs7412", names(geno), value = TRUE)[1], "rs7412_dose")

# rs429358_dose = T-allele dosage; rs7412_dose = C-allele dosage (see header note)
geno[, `:=`(e4_count = 2 - rs429358_dose,
            e2_count = 2 - rs7412_dose)]
geno[, e4_carrier := fifelse(e4_count > 0, 1L, 0L)]
geno[, APOE_genotype := fcase(
  rs429358_dose == 2 & rs7412_dose == 0, "e2/e2",
  rs429358_dose == 2 & rs7412_dose == 1, "e2/e3",
  rs429358_dose == 2 & rs7412_dose == 2, "e3/e3",
  rs429358_dose == 1 & rs7412_dose == 1, "e2/e4",
  rs429358_dose == 1 & rs7412_dose == 2, "e3/e4",
  rs429358_dose == 0 & rs7412_dose == 2, "e4/e4",
  default = NA_character_)]
print(geno[, .N, by = APOE_genotype][order(-N)])

# ---- merge phenotype/covariates --------------------------------------------
pheno <- fread(PHENO_FILE)
data <- merge(pheno, geno[, .(IID, e4_count, e2_count, e4_carrier, APOE_genotype)],
              by = "IID")
cat(sprintf("N=%d (case=%d, control=%d)\n",
            nrow(data), sum(data$pheno == 1, na.rm = TRUE), sum(data$pheno == 0, na.rm = TRUE)))
data[, APOE_genotype := factor(APOE_genotype,
        levels = c("e3/e3", "e2/e2", "e2/e3", "e2/e4", "e3/e4", "e4/e4"))]

# ---- association models (the primary covariate set) ------------------------
covar <- "age + sex_plink + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + smoking + drink + surg_taste_affecting_full"
fit <- function(term) glm(as.formula(paste("pheno ~", term, "+", covar)),
                          data = data, family = binomial)

carrier <- as.data.table(tidy(fit("e4_carrier"), conf.int = TRUE, exponentiate = TRUE))[term == "e4_carrier"]
dose    <- as.data.table(tidy(fit("e4_count"),   conf.int = TRUE, exponentiate = TRUE))[term == "e4_count"]
gtype   <- as.data.table(tidy(fit("APOE_genotype"), conf.int = TRUE, exponentiate = TRUE))[grepl("APOE_genotype", term)]
gtype[, Comparison := sub("APOE_genotype", "", term)]

cat(sprintf("\ne4 carrier : OR=%.3f (%.3f-%.3f) P=%.2e\n",
            carrier$estimate, carrier$conf.low, carrier$conf.high, carrier$p.value))
cat(sprintf("per e4 dose: OR=%.3f (%.3f-%.3f) P=%.2e\n",
            dose$estimate, dose$conf.low, dose$conf.high, dose$p.value))

# ---- assemble forest table -------------------------------------------------
mk <- function(comp, e) data.table(Comparison = comp, OR = e$estimate,
                                   CI_lower = e$conf.low, CI_upper = e$conf.high, P = e$p.value)
all_results <- rbindlist(list(
  data.table(Comparison = "e3/e3 (reference)", OR = 1, CI_lower = 1, CI_upper = 1, P = NA_real_),
  gtype[, .(Comparison, OR = estimate, CI_lower = conf.low, CI_upper = conf.high, P = p.value)],
  data.table(Comparison = "---", OR = NA, CI_lower = NA, CI_upper = NA, P = NA),
  mk("e4 carrier vs non-carrier", carrier),
  mk("Per e4 allele (trend)", dose)), use.names = TRUE)

# Export and dx upload to RAP  (the Results/*.csv below feed Fig 3d + supp tables)
fwrite(all_results, file.path(RES_DIR, "plot_data_association_results.csv"))

# ---- frequency + sample-size tables ----------------------------------------
gf <- data[!is.na(APOE_genotype) & !is.na(pheno),
           .(n = .N), by = .(pheno, APOE_genotype)]
gf[, `:=`(pct = 100 * n / sum(n), Group = fifelse(pheno == 1, "Case", "Control")), by = pheno]
fwrite(gf, file.path(RES_DIR, "plot_data_genotype_freq.csv"))

cf <- data[!is.na(e4_carrier) & !is.na(pheno),
           .(n_total = .N, n_carrier = sum(e4_carrier == 1),
             carrier_pct = 100 * mean(e4_carrier == 1)), by = pheno]
cf[, Group := fifelse(pheno == 1, "Case", "Control")]
fwrite(cf, file.path(RES_DIR, "plot_data_e4_carrier_freq.csv"))

fwrite(data[!is.na(APOE_genotype) & !is.na(pheno),
            .(N_total = .N, N_case = sum(pheno == 1), N_control = sum(pheno == 0),
              N_e4_carrier = sum(e4_carrier == 1, na.rm = TRUE))],
       file.path(RES_DIR, "sample_summary.csv"))

cat("\nAPOE association tables ->", RES_DIR, "\n")
