#!/usr/bin/env Rscript
# Follow-up questionnaire exposures, block 6: work environment (Category 130). Each exposure is
# "ever" across the lifetime job history (up to 40 jobs); four chemical exposures are also summed.
# Runs in the R session of 01.
library(data.table)

WORK_FIELDS <- c(work_noisy = 22606, work_dusty = 22609, work_chemical_fumes = 22610,
                 work_cigarette_smoke = 22611, work_paints_thinners = 22613,
                 work_pesticides = 22614, work_diesel = 22615, work_breathing_problems = 22616,
                 work_shift_work_ever = 22620, work_night_shifts_ever = 22650)
EXPOSED   <- c("Sometimes", "Often", "Always", "2", "3", "4")
UNEXPOSED <- c("Rarely/never", "Rarely / never", "Never", "1")

# TRUE where any job slot carries a matching answer
any_job <- function(job_cols, levels, numeric_match) Reduce(`|`, lapply(job_cols, function(cn) {
  x <- raw[[cn]]
  if (is.numeric(x)) !is.na(x) & numeric_match(x) else as.character(x) %in% levels
}))

work <- raw[, .(eid)]
for (v in names(WORK_FIELDS)) {
  job_cols <- intersect(paste0("p", WORK_FIELDS[[v]], "_a", 0:39), names(raw))
  work[, (v) := fifelse(any_job(job_cols, EXPOSED, function(x) x >= 2), 1L,
                        fifelse(any_job(job_cols, UNEXPOSED, function(x) x == 1), 0L, NA_integer_))]
}
work[, work_chemical_score := fifelse(rowSums(!is.na(.SD)) > 0, as.integer(rowSums(.SD, na.rm = TRUE)),
                                      NA_integer_),
     .SDcols = c("work_chemical_fumes", "work_paints_thinners", "work_pesticides", "work_diesel")]

save_section(work, "Work environment",
             c(setNames(rep("binary", length(WORK_FIELDS)), names(WORK_FIELDS)),
               work_chemical_score = "continuous (0-4)"))
