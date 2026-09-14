#!/usr/bin/env Rscript
#> out  output/evidence_tiering/exwas_continuous_SD.csv
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
suppressPackageStartupMessages(library(data.table))

OUT_DIR <- "output/evidence_tiering"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

DICTS <- c("input/analysis_ready/baseline_exwas_variable_dictionary.csv",
           "input/analysis_ready/followup_exwas_variable_dictionary.csv")
DATA  <- c("input/analysis_ready/exwas_baseline_group1.csv",
           "input/analysis_ready/exwas_followup_group1.csv")
for (f in c(DICTS, DATA)) if (!file.exists(f)) stop("input not found: ", f)

dict <- unique(rbindlist(lapply(DICTS, fread, select = c("var_name", "var_type"))),
               by = "var_name")
cat("dictionary entries:", nrow(dict), "\n")
print(dict[, .N, by = var_type][order(-N)])

dat <- lapply(DATA, fread)

# Pattern match, not equality: the dictionary carries three continuous
# labels, and a variable with no SD here silently stays on its raw scale.
cont <- dict[grepl("continuous", var_type, ignore.case = TRUE), var_name]
cat("continuous exposures declared:", length(cont), "\n")

sd_of <- function(d, v) {
  if (!v %in% names(d)) return(NA_real_)
  sd(suppressWarnings(as.numeric(d[[v]])), na.rm = TRUE)
}

lk <- rbindlist(lapply(cont, function(v) {
  s <- NA_real_
  for (d in dat) { s <- sd_of(d, v); if (!is.na(s)) break }
  data.table(variable = v, var_type = "continuous", sd = s)
}))

bad <- lk[is.na(sd) | sd <= 0]
if (nrow(bad)) {
  cat("dropped (absent or zero-variance):", nrow(bad), "\n")
  print(bad[, .(variable)])
}
lk <- lk[!is.na(sd) & sd > 0]
cat("continuous exposures with a usable SD:", nrow(lk), "\n")

fwrite(lk, file.path(OUT_DIR, "exwas_continuous_SD.csv"))
# Export and dx upload to RAP  (Figure 1 rescales continuous exposures with this)
cat("wrote", file.path(OUT_DIR, "exwas_continuous_SD.csv"), "\n")
