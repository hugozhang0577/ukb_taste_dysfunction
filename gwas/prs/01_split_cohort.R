# Random 70/30 split of the primary GWAS sample, stratified by case status: the 70% discovery subsample
# for the PRS GWAS and the 30% held-out subsample in which the score is computed and evaluated.
library(data.table)
setwd(Sys.getenv("PROJECT_DIR", "."))
OUT <- "output/prs"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

ph <- fread("output/gwas_sample_qc/taste_2w_strict_white_gwas_pheno.txt")[!is.na(pheno)]
set.seed(20260413L)
ph[, .disc := FALSE]
for (lvl in sort(unique(ph$pheno))) {
  idx <- which(ph$pheno == lvl)
  ph[sample(idx, size = floor(0.70 * length(idx))), .disc := TRUE]
}
heldout <- ph[.disc == FALSE][, .disc := NULL]
fwrite(ph[.disc == TRUE][, .disc := NULL], file.path(OUT, "discovery_70_gwas_pheno.txt"), sep = "\t")
fwrite(heldout, file.path(OUT, "pheno_validation_subset.txt"), sep = "\t")

# PRSice-2 reads the sample IDs embedded in the plink2-written BGEN files, which have the form FID_IID
prsice_ids <- heldout[, .(FID = paste0(FID, "_", IID), IID = paste0(FID, "_", IID), pheno)]
fwrite(prsice_ids, file.path(OUT, "heldout_prsice_pheno.txt"), sep = "\t")
fwrite(prsice_ids[, .(FID, IID)], file.path(OUT, "heldout_keep.txt"), sep = "\t", col.names = FALSE)
