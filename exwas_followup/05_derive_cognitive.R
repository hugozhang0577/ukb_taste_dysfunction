#!/usr/bin/env Rscript
# Follow-up questionnaire exposures, block 5: cognitive function (Category 116). Timed and error
# scores are log-transformed; the online follow-up (instance 1) is used where available, otherwise
# the first online assessment (instance 0). Runs in the R session of 01.
library(data.table)

num <- function(x) suppressWarnings(as.numeric(x))
latest <- function(i1, i0) fifelse(!is.na(i1), num(i1), num(i0))
finite_or_na <- function(x) fifelse(is.finite(x), x, NA_real_)

cog <- raw[, .(eid,
               tmt_a = latest(p20156_i1, p20156_i0),                 # trail making A (s)
               tmt_b = latest(p20157_i1, p20157_i0),                 # trail making B (s)
               cog_symbol_digit = latest(p20159_i1, p20159_i0),
               cog_fluid_iq = latest(p20191_i1, p20191_i0),
               cog_matrix_reasoning = num(p20760_i1))]               # instance 1 only
cog[, pairs_errors := rowSums(raw[, .(p20132_i0_a0, p20132_i0_a1, p20132_i0_a2)])]  # instance 0 only

cog[, `:=`(cog_tmt_a = finite_or_na(log(tmt_a)),
           cog_tmt_b = finite_or_na(log(tmt_b)),
           cog_tmt_ba = fifelse(tmt_b - tmt_a > 0, log(tmt_b - tmt_a), NA_real_),
           cog_pairs_errors = log(pairs_errors + 1),
           cog_symbol_digit = pmin(pmax(cog_symbol_digit, 1), 40))]  # truncated to 1-40

save_section(cog, "Cognitive function", c(
  cog_tmt_a = "continuous", cog_tmt_b = "continuous", cog_tmt_ba = "derived_continuous",
  cog_symbol_digit = "continuous", cog_pairs_errors = "continuous",
  cog_fluid_iq = "continuous", cog_matrix_reasoning = "continuous"))
