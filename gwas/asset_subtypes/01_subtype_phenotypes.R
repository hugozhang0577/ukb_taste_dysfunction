# Phenotype files for the per-subtype GWAS: the cases of one subtype (pheno = 1) against the shared controls of
# the primary GWAS (pheno = 0); the other cases are left out. Covariates are those of the primary GWAS.
library(data.table)
source("unsupervised_subtyping/_subtype_map.R")   # SUBTYPE_MAP: cluster id -> subtype letter, per sex
setwd(Sys.getenv("PROJECT_DIR", "."))
OUT <- "output/gwas_subtypes"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

labels <- rbindlist(lapply(c("m", "f"), function(sx) {
  fl <- readRDS(sprintf("output/subtyping/clusters/cluster_assignments_g1_%s_k4.rds", sx))$final_labels
  data.table(IID = as.integer(names(fl)), subtype = unname(SUBTYPE_MAP[[sx]][as.character(as.integer(fl))]))
}))
gwas_pheno <- fread("output/gwas_sample_qc/taste_2w_strict_white_gwas_pheno.txt")

for (s in c("A", "B", "C", "D")) {
  d <- copy(gwas_pheno)[, pheno := fifelse(IID %in% labels[subtype == s, IID], 1L, fifelse(pheno == 0, 0L, NA_integer_))]
  fwrite(d[!is.na(pheno)], file.path(OUT, sprintf("pheno_subtype_%s.txt", s)), sep = "\t")
}
