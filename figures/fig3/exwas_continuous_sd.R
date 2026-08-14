#!/usr/bin/env Rscript
# =============================================================================
# exwas_continuous_sd.R  (per-SD scale for the continuous exposures)
# =============================================================================
#
# Standard deviation of every continuous exposure in the discovery cohort, so
# that a continuous exposure's effect can be expressed per one standard
# deviation instead of per raw unit.
#
# Why this is needed. An exposure-wide scan puts variables measured on unrelated
# scales onto one axis. Per raw unit, a waist-to-hip ratio odds ratio is per unit
# of a quantity whose whole range is well under one, and a cell-count odds ratio
# is per single cell; the first comes out around 7.7 and the second in the
# thousands, and neither number means what it appears to mean. Rescaling to
# per-SD makes the effects comparable and stops the axis being set by whichever
# variable happens to have the smallest unit.
#
# This is applied at the figure, not upstream: the tier flags in
# tier_flagged_results.R are computed on each scan's own beta, so this lookup
# changes how an effect is displayed and never which rows passed.
#
# Binary and ordinal exposures are left alone. Their effects are already per
# level, and a "standard deviation" of a 0/1 variable is not a step anyone takes.
#
# Input:  input/analysis_ready/exwas_{baseline,followup}_group1.csv
#         input/analysis_ready/exwas_{baseline,followup}_dictionary.csv
# Output: output/evidence_tiering/exwas_continuous_SD.csv
# =============================================================================

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
suppressPackageStartupMessages(library(data.table))

OUT_DIR <- "output/evidence_tiering"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

DICTS <- c("input/analysis_ready/exwas_baseline_dictionary.csv",
           "input/analysis_ready/exwas_followup_dictionary.csv")
DATA  <- c("input/analysis_ready/exwas_baseline_group1.csv",
           "input/analysis_ready/exwas_followup_group1.csv")
for (f in c(DICTS, DATA)) if (!file.exists(f)) stop("input not found: ", f)

dict <- unique(rbindlist(lapply(DICTS, fread, select = c("var_name", "var_type"))),
               by = "var_name")
cat("dictionary entries:", nrow(dict), "\n")
print(dict[, .N, by = var_type][order(-N)])

dat <- lapply(DATA, fread)

cont <- dict[var_type == "continuous", var_name]
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

# A declared-continuous variable with no SD is either absent from both matrices
# or constant within the cohort. Either way it cannot be rescaled, so it is
# dropped here and named, rather than silently carried through as NA.
bad <- lk[is.na(sd) | sd <= 0]
if (nrow(bad)) {
  cat("dropped (absent or zero-variance):", nrow(bad), "\n")
  print(bad[, .(variable)])
}
lk <- lk[!is.na(sd) & sd > 0]
cat("continuous exposures with a usable SD:", nrow(lk), "\n")

fwrite(lk, file.path(OUT_DIR, "exwas_continuous_SD.csv"))
# Export and dx upload to RAP  (Figure 3 rescales continuous exposures with this)
cat("wrote", file.path(OUT_DIR, "exwas_continuous_SD.csv"), "\n")
