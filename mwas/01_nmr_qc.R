# Remove technical variation from the Nightingale NMR biomarkers with ukbnmr (algorithm
# version 3, outlier plates removed) and derive the extended biomarker ratios.
library(data.table)
library(ukbnmr)
setwd(Sys.getenv("PROJECT_DIR", "."))

OUT_DIR <- "output/metabolomics_qc"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

exported  <- fread("input/metabolomics/nmr_merged_all.csv", sep = ",", header = TRUE)
processed <- remove_technical_variation(x = exported, remove.outlier.plates = TRUE, version = 3L)

fwrite(processed$biomarkers, file.path(OUT_DIR, "nmr_biomarkers_qc.csv"))
fwrite(compute_extended_ratios(processed$biomarkers), file.path(OUT_DIR, "nmr_extended_ratios.csv"))
