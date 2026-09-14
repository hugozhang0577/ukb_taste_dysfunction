#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

files <- sprintf("prs_chr%d.all_score", 1:22)
missing <- files[!file.exists(files)]
if (length(missing))
  stop("no score file for: ", paste(missing, collapse = ", "),
       "\n  re-run 04_score_by_chr.sh for those chromosomes")

merged <- fread(files[1], colClasses = list(character = c("FID", "IID")))
score_cols <- setdiff(names(merged), c("FID", "IID"))
cat(sprintf("chr1: %d samples, %d thresholds\n", nrow(merged), length(score_cols)))

for (f in files[-1]) {
  d <- fread(f, colClasses = list(character = c("FID", "IID")))
  if (!identical(sort(score_cols), sort(setdiff(names(d), c("FID", "IID")))))
    stop(f, " has a different set of thresholds than chr1")
  n0 <- nrow(merged)
  merged <- merge(merged, d, by = c("FID", "IID"), suffixes = c("", ".add"))
  if (nrow(merged) != n0)
    cat(sprintf("  [%s] samples %d -> %d after the join\n", f, n0, nrow(merged)))
  for (cc in score_cols) {
    set(merged, j = cc, value = merged[[cc]] + merged[[paste0(cc, ".add")]])
    set(merged, j = paste0(cc, ".add"), value = NULL)
  }
}
cat(sprintf("genome-wide: %d samples, %d thresholds\n", nrow(merged), length(score_cols)))

fwrite(merged, "prs_combined_all_scores.txt", sep = "\t")
# Export and dx upload to RAP  (the genome-wide scores 06_regression.R reads)

# The clumped variants actually contributing, for the count reported in the text.
snp_files <- sprintf("prs_chr%d.snp", 1:22)
snp_files <- snp_files[file.exists(snp_files)]
if (length(snp_files)) {
  snps <- rbindlist(lapply(snp_files, fread), fill = TRUE)
  fwrite(snps, "prs_snps_used.txt", sep = "\t")
  cat(sprintf("clumped variants retained: %d\n", nrow(snps)))
} else {
  cat("no .snp files found; skipping the variant list\n")
}
