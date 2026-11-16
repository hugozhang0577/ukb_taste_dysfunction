# Olink proteomics QC for one cohort (argument: group1, group2 or group3; default group1):
# baseline-visit NPX of cohort members, protein (> 30% missing) and sample (> 50% missing)
# filters, then imputation of below-LOD values (MinProb) and missing values (protein median).
library(data.table)
setwd(Sys.getenv("PROJECT_DIR", "."))

COHORT           <- c(commandArgs(trailingOnly = TRUE), "group1")[1]
PROTEIN_MISS_MAX <- 0.30
SAMPLE_MISS_MAX  <- 0.50
MINPROB_Q        <- 0.01   # MinProb: draws from N(1st percentile of observed, 0.3 x SD of observed)
MINPROB_SD_SCALE <- 0.3
OUT_DIR <- "output/pwas/olink_qc"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(42)

pheno_eids <- unique(fread(sprintf("input/analysis_ready/phenotype_%s.csv", COHORT), select = "eid")$eid)
prot <- fread("input/proteomics/olink_full_dataset.csv", drop = "V1", na.strings = c("", "NA", "NaN"))
prot <- prot[ins_index == 0 & eid %in% pheno_eids]

# Wide matrices of NPX and limit of detection; a participant-protein pair with no row is missing
npx <- dcast(prot[, .(eid, protein_id, result)], eid ~ protein_id, value.var = "result")
lod <- dcast(prot[, .(eid, protein_id, LOD)], eid ~ protein_id, value.var = "LOD",
             fun.aggregate = function(x) median(x, na.rm = TRUE))

# Protein filter, then sample filter on the retained proteins
# (the sample filter mainly removes assays from the half-panel batch 7)
proteins      <- setdiff(names(npx), "eid")
protein_miss  <- sapply(npx[, ..proteins], function(x) mean(is.na(x)))
keep_proteins <- names(protein_miss)[protein_miss <= PROTEIN_MISS_MAX]
npx <- npx[, c("eid", keep_proteins), with = FALSE]
lod <- lod[, c("eid", keep_proteins), with = FALSE]
keep_rows <- rowMeans(is.na(npx[, ..keep_proteins])) <= SAMPLE_MISS_MAX
npx <- npx[keep_rows]
lod <- lod[keep_rows]

# Imputation per protein: below-LOD values (left-censored) by MinProb, or LOD/2 if fewer than
# 5 values are above the LOD; missing values by the median after the below-LOD step
for (p in keep_proteins) {
  x <- npx[[p]]
  l <- lod[[p]]
  below <- !is.na(x) & !is.na(l) & x < l
  miss  <- is.na(x)
  if (any(below)) {
    obs <- x[!miss & !below]
    if (length(obs) >= 5) {
      x[below] <- rnorm(sum(below), mean = quantile(obs, MINPROB_Q), sd = sd(obs) * MINPROB_SD_SCALE)
    } else {
      x[below] <- median(l, na.rm = TRUE) / 2
    }
  }
  if (any(miss)) x[miss] <- median(x[!miss], na.rm = TRUE)
  set(npx, j = p, value = x)
}

fwrite(npx, file.path(OUT_DIR, sprintf("proteomics_%s.csv", COHORT)))
# Copied to input/analysis_ready/ for 02_run_pwas.R and later steps
