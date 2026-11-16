# Diagnosis-count covariate: number of distinct three-character ICD-10 codes
# recorded before the taste assessment (row sum of the exposure matrix).
library(data.table)
setwd(Sys.getenv("PROJECT_DIR", "."))

for (g in c("group1", "group2", "group3")) {
  mat <- readRDS(sprintf("output/dwas/derive/icd10_matrix_%s.rds", g))
  dx  <- integer(nrow(mat))
  for (cd in setdiff(names(mat), "eid")) dx <- dx + as.integer(mat[[cd]])
  fwrite(data.table(eid = mat$eid, dx_count_pre = dx),
         sprintf("output/dwas/derive/dx_count_pre_%s.csv", g))
}
# Copied to input/analysis_ready/ for 04_run_dwas.R
