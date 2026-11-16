#!/usr/bin/env Rscript
# Follow-up questionnaire exposures, assembly: merge the seven questionnaire blocks into one
# exposure matrix per cohort and build the variable dictionary that defines the scan and its
# FDR families. Both outputs are copied to input/analysis_ready/ for the scan (08).
library(data.table)

setwd(Sys.getenv("PROJECT_DIR", "."))
IN_DIR  <- "output/followup_exwas/derive/per_section"
OUT_DIR <- "output/followup_exwas/derive"

FAMILIES <- c("Food preferences", "Diet (24-hour recall)", "Digestive health", "Mental health",
              "Experience of pain", "Cognitive function", "Work environment")
stem <- file.path(IN_DIR, paste0("followup_exwas_",
                                 gsub("^_|_$", "", gsub("[^a-z0-9]+", "_", tolower(FAMILIES)))))

dict <- rbindlist(lapply(paste0(stem, "_variable_dict.csv"), fread))
# held out of the scan: its "interference with life" wording overlaps the outcome's
# functional-impact criterion
dict[, exclude_primary := var_name == "gi_life_interference"]
fwrite(dict, file.path(OUT_DIR, "followup_exwas_variable_dictionary.csv"))

for (g in c("group1", "group2", "group3")) {
  blocks <- lapply(paste0(stem, "_", g, ".csv"), fread)
  fwrite(Reduce(function(a, b) merge(a, b, by = "eid", all = TRUE), blocks),
         file.path(OUT_DIR, paste0("exwas_followup_", g, ".csv")))
}
