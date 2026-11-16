# Per-cohort covariate tables: each cohort membership file (group1 discovery; group2, group3
# held-out; full and Olink subsets) left-joined to the base-table covariates.
library(data.table)
setwd(Sys.getenv("PROJECT_DIR", "."))
dir.create("output/analysis_ready", recursive = TRUE, showWarnings = FALSE)

covar <- fread("input/analysis_ready/base_table_covar.csv")
for (cohort in c("group1_full", "group2_full", "group3_full",
                 "group1_olink", "group2_olink", "group3_olink")) {
  d <- merge(fread(file.path("input/eids", paste0(cohort, ".csv"))), covar,
             by = "eid", all.x = TRUE, sort = FALSE)
  fwrite(d, file.path("output/analysis_ready", paste0(cohort, "_covariates.csv")))
}
