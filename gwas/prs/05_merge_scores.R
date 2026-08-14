#!/usr/bin/env Rscript
# =============================================================================
# PRS — sum the per-chromosome scores into one genome-wide score
# =============================================================================
#
# A polygenic score is a sum over variants, and the variants were partitioned by
# chromosome, so the genome-wide score is the sum of the per-chromosome scores at
# the same threshold. Nothing is averaged and nothing is re-standardised here;
# standardisation happens once, in 06_regression.R, on the held-out subset.
#
# Every chromosome must be present. A score silently missing a chromosome is
# still a plausible-looking number, so a missing file stops the merge rather than
# quietly producing a partial score.
#
# The per-chromosome files must also agree on their samples. They are joined on
# FID/IID rather than assumed to be row-aligned: PRSice writes samples in target
# order, and two chromosomes whose bgen sample files differ would otherwise be
# added row-by-row across different people.
#
# Usage: Rscript 05_merge_scores.R
# Input:  prs_chr{1..22}.all_score
# Output: prs_combined_all_scores.txt
#         prs_snps_used.txt   (the clumped variants, for the count in the text)
# =============================================================================

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
