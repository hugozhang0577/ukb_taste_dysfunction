# APOE diplotypes from rs429358/rs7412 and their association with the primary outcome: logistic regression
# for each diplotype vs e3/e3, e4 carriage and e4 dose, adjusted for age, sex, ten genetic principal components,
# smoking, alcohol and taste-affecting surgery; diplotype and carrier frequencies by case status.
library(data.table)
library(broom)
setwd(Sys.getenv("PROJECT_DIR", "."))
RES_DIR <- "output/apoe_diplotype"
dir.create(RES_DIR, showWarnings = FALSE, recursive = TRUE)

geno <- fread("input/gwas_genotype/apoe_genotypes.raw")
setnames(geno, grep("rs429358", names(geno), value = TRUE)[1], "rs429358_dose")
setnames(geno, grep("rs7412", names(geno), value = TRUE)[1], "rs7412_dose")
# rs429358_dose counts the T allele and rs7412_dose the C allele (the counted allele is the column-name suffix)
geno[, `:=`(e4_count = 2 - rs429358_dose, e2_count = 2 - rs7412_dose)]
geno[, e4_carrier := fifelse(e4_count > 0, 1L, 0L)]
geno[, APOE_genotype := fcase(
  rs429358_dose == 2 & rs7412_dose == 0, "e2/e2",
  rs429358_dose == 2 & rs7412_dose == 1, "e2/e3",
  rs429358_dose == 2 & rs7412_dose == 2, "e3/e3",
  rs429358_dose == 1 & rs7412_dose == 1, "e2/e4",
  rs429358_dose == 1 & rs7412_dose == 2, "e3/e4",
  rs429358_dose == 0 & rs7412_dose == 2, "e4/e4",
  default = NA_character_)]

data <- merge(fread("output/gwas_sample_qc/taste_2w_strict_white_gwas_pheno.txt"),
              geno[, .(IID, e4_count, e2_count, e4_carrier, APOE_genotype)], by = "IID")
data[, APOE_genotype := factor(APOE_genotype, levels = c("e3/e3", "e2/e2", "e2/e3", "e2/e4", "e3/e4", "e4/e4"))]

covar <- "age + sex_plink + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + smoking + drink + surg_taste_affecting_full"
fit <- function(term) as.data.table(tidy(glm(as.formula(paste("pheno ~", term, "+", covar)), data = data, family = binomial),
                                         conf.int = TRUE, exponentiate = TRUE))
carrier <- fit("e4_carrier")[term == "e4_carrier"]
dose    <- fit("e4_count")[term == "e4_count"]
gtype   <- fit("APOE_genotype")[grepl("APOE_genotype", term)][, Comparison := sub("APOE_genotype", "", term)]

mk <- function(comp, e) data.table(Comparison = comp, OR = e$estimate, CI_lower = e$conf.low, CI_upper = e$conf.high, P = e$p.value)
fwrite(rbindlist(list(
  data.table(Comparison = "e3/e3 (reference)", OR = 1, CI_lower = 1, CI_upper = 1, P = NA_real_),
  gtype[, .(Comparison, OR = estimate, CI_lower = conf.low, CI_upper = conf.high, P = p.value)],
  data.table(Comparison = "---", OR = NA, CI_lower = NA, CI_upper = NA, P = NA),
  mk("e4 carrier vs non-carrier", carrier),
  mk("Per e4 allele (trend)", dose)), use.names = TRUE),
  file.path(RES_DIR, "plot_data_association_results.csv"))

gf <- data[!is.na(APOE_genotype) & !is.na(pheno), .(n = .N), by = .(pheno, APOE_genotype)]
gf[, `:=`(pct = 100 * n / sum(n), Group = fifelse(pheno == 1, "Case", "Control")), by = pheno]
fwrite(gf, file.path(RES_DIR, "plot_data_genotype_freq.csv"))

cf <- data[!is.na(e4_carrier) & !is.na(pheno), .(n_total = .N, n_carrier = sum(e4_carrier == 1),
                                                  carrier_pct = 100 * mean(e4_carrier == 1)), by = pheno]
cf[, Group := fifelse(pheno == 1, "Case", "Control")]
fwrite(cf, file.path(RES_DIR, "plot_data_e4_carrier_freq.csv"))

fwrite(data[!is.na(APOE_genotype) & !is.na(pheno), .(N_total = .N, N_case = sum(pheno == 1), N_control = sum(pheno == 0),
                                                     N_e4_carrier = sum(e4_carrier == 1, na.rm = TRUE))],
       file.path(RES_DIR, "sample_summary.csv"))
