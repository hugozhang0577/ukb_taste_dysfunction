#> ship output/analysis_ready/{group}_covariates.csv -> input/analysis_ready/phenotype_{group}.csv
suppressPackageStartupMessages({
  library(data.table)
})

left_join_dt <- function(x, y, by = "eid") {
  merge(x, y, by = by, all.x = TRUE, sort = FALSE)
}

merge_covariate_tables <- function(cohort_file,
                                   covariate_files,
                                   output_file,
                                   by = "eid") {
  cohort <- fread(cohort_file)
  if (!by %in% names(cohort)) stop("Cohort file is missing key column: ", by)

  original_n <- nrow(cohort)
  added <- data.table(file = character(), n_rows = integer(), n_columns = integer())

  for (f in covariate_files) {
    cov <- fread(f)
    if (!by %in% names(cov)) stop("Covariate file is missing key column: ", f)
    before_cols <- names(cohort)
    cohort <- left_join_dt(cohort, cov, by = by)
    added <- rbind(added, data.table(
      file = basename(f),
      n_rows = nrow(cov),
      n_columns = length(setdiff(names(cohort), before_cols))
    ))
  }

  if (nrow(cohort) != original_n) {
    stop("Row count changed during left joins: ", original_n, " -> ", nrow(cohort))
  }

  fwrite(cohort, output_file)
  # Export and dx upload to RAP  (the per-cohort covariate tables the scans read)
  list(output_file = output_file, n_rows = nrow(cohort), n_columns = ncol(cohort), added = added)
}

if (sys.nframe() == 0) {
  project_dir <- Sys.getenv("PROJECT_DIR", unset = ".")
  cohort_dir <- file.path(project_dir, "input", "eids")
  out_dir <- file.path(project_dir, "output", "analysis_ready")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cohort_files <- file.path(cohort_dir, c(
    "group1_full.csv", "group2_full.csv", "group3_full.csv",
    "group1_olink.csv", "group2_olink.csv", "group3_olink.csv"
  ))

  covariate_files <- file.path(project_dir, "input", "analysis_ready",
                               "base_table_covar.csv")

  for (f in cohort_files) {
    if (!file.exists(f)) next
    out <- file.path(out_dir, sub("\\.csv$", "_covariates.csv", basename(f)))
    log <- merge_covariate_tables(f, covariate_files, out)
    print(log)
  }
}
