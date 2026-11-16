# Genome-wide polygenic scores: sum of the per-chromosome PRSice-2 scores at each P-value threshold,
# and the list of clumped variants contributing.
library(data.table)
setwd(Sys.getenv("PROJECT_DIR", "."))
IN <- "input/prs"
LEVELS <- c(5e-08, 1e-07, 1e-06, 1e-05, 1e-04, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1)

# PRSice writes a column only where a chromosome gains variants; its score at any other level equals
# the score at the nearest lower level written (0 if none)
m <- NULL
for (chr in 1:22) {
  d <- fread(file.path(IN, sprintf("prs_chr%d.all_score", chr)))
  if (is.null(m)) { m <- d[, .(FID, IID)]; for (l in LEVELS) m[, paste0("Pt_", l) := 0] }
  cols <- grep("^Pt_", names(d), value = TRUE)
  have <- as.numeric(sub("^Pt_", "", cols))
  for (l in LEVELS) {
    below <- which(have <= l * (1 + 1e-9))
    if (length(below)) set(m, j = paste0("Pt_", l), value = m[[paste0("Pt_", l)]] + d[[cols[below[which.max(have[below])]]]])
  }
}
m <- m[, c("FID", "IID", names(m)[-(1:2)][sapply(m[, -(1:2)], function(x) any(x != 0))]), with = FALSE]
fwrite(m, "output/prs/prs_combined_all_scores.txt", sep = "\t")

fwrite(rbindlist(lapply(sprintf(file.path(IN, "prs_chr%d.snp"), 1:22), fread), fill = TRUE),
       "output/prs/prs_snps_used.txt", sep = "\t")
