#!/usr/bin/env Rscript
# =============================================================================
# SAIGE — Manhattan and QQ plots, and the genomic-control lambda
# =============================================================================
#
# Diagnostics for the primary GWAS. It reads the merged summary statistics
# written by 03_merge_sumstats.sh, computes the genomic-control lambda (all
# variants, MAF>0.01, MAF>0.05), and writes a Manhattan plot, a QQ plot, the
# top-10 hits, and a row in all_runs_summary.csv.
#
# The summary-statistics file is taken from an environment variable, so a
# sensitivity run is this same script pointed at that run's merged file with a
# different RUN_TAG. Nothing about the diagnostics differs between runs, which
# is why there is no per-run branching here.
#
# lambda is reported in Results for the genome-wide analysis; this script only
# computes and records it. It is not reported for the panel scans, where an
# outcome-anchored panel makes a high value expected rather than informative.
#
# Memory: one run at a time, rm() + gc() between runs; scattermore for
# rasterised points; ASCII-only glyphs so cairo_pdf renders cleanly.
#
# Usage:
#   PROJECT_DIR=/path/to/project Rscript 04_manhattan_qq_lambda.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scattermore)
  library(ggrepel)
  library(scales)
})

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)

# ---- the run to diagnose ---------------------------------------------------
# Defaults to the primary scan. For a sensitivity run, set both variables:
#   SUMSTATS=/path/to/that/run/gwas_sumstats.tsv.gz RUN_TAG=sens_4w_strict
SUMSTATS <- Sys.getenv("SUMSTATS", unset = file.path(
  PROJECT_DIR, "gwas", "cohort_primary", "SAIGE", "step2", "gwas_sumstats.tsv.gz"))
RUN_TAG  <- Sys.getenv("RUN_TAG", unset = "primary")
if (!file.exists(SUMSTATS))
  stop("summary statistics not found -- run 03_merge_sumstats.sh first: ", SUMSTATS)

gwas_jobs <- list(list(tag = RUN_TAG, file = SUMSTATS))

out_root <- file.path(PROJECT_DIR, "output", "gwas_plots")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# ---- helpers ---------------------------------------------------------------
calc_lambda <- function(pvals) {
  chisq <- qchisq(1 - pvals, 1)
  median(chisq, na.rm = TRUE) / qchisq(0.5, 1)
}

read_saige <- function(path) {
  hdr <- names(fread(path, nrows = 0))
  need_cols <- c("CHR", "POS", "MarkerID", "Allele1", "Allele2",
                 "AF_Allele2", "BETA", "SE", "p.value", "N")
  dt <- fread(path, nThread = 4, select = intersect(need_cols, hdr))
  if (!"N" %in% names(dt)) dt[, N := NA_integer_]
  dt <- dt[!is.na(p.value) & p.value > 0 & p.value < 1]
  dt[, CHR := as.integer(CHR)]
  dt[p.value == 0, p.value := .Machine$double.xmin]
  dt
}

read_gwas <- function(path) read_saige(path)

make_manhattan <- function(gwas, tag, out_dir) {
  threshold_gw <- 5e-8; threshold_suggestive <- 1e-5
  gm <- gwas[!is.na(CHR) & CHR %in% 1:22][order(CHR, POS)]

  chr_info <- gm[, .(chr_len = max(POS)), by = CHR][order(CHR)]
  chr_info[, cumsum_len := cumsum(as.numeric(chr_len)) - chr_len]
  chr_info[, center := cumsum_len + chr_len / 2]

  gm <- merge(gm, chr_info[, .(CHR, cumsum_len)], by = "CHR", sort = FALSE)
  gm[, BP_cum := POS + cumsum_len]
  gm[, log10p := -log10(p.value)]
  gm[, snp_category := fifelse(p.value < threshold_gw, "genome_wide",
                       fifelse(p.value < threshold_suggestive, "suggestive", "background"))]

  n_gw  <- sum(gm$snp_category == "genome_wide")
  n_sug <- sum(gm$snp_category == "suggestive")
  cat("  GW:", n_gw, " Suggestive:", n_sug, "\n")

  bg_data         <- gm[snp_category == "background"][, chr_color := fifelse(CHR %% 2 == 1, "odd", "even")]
  suggestive_data <- gm[snp_category == "suggestive"]
  gw_data         <- gm[snp_category == "genome_wide"]

  top10 <- gm[order(p.value)][1:min(10, nrow(gm))]
  fwrite(top10[, .(MarkerID, CHR, POS, Allele1, Allele2, AF_Allele2, BETA, SE, p.value)],
         file.path(out_dir, "top10_snps.csv"))

  color_odd <- "#6BAED6"; color_even <- "#08519C"
  color_suggestive <- "#FF8C00"; color_gw <- "#E31A1C"
  y_max <- max(gm$log10p, na.rm = TRUE)

  # piecewise y-axis: compress the non-significant band, stretch the sig band
  break_low <- -log10(threshold_suggestive)   # 5
  y_top <- max(y_max * 1.08, break_low + 2)
  compress_frac <- 0.30
  squish_trans <- trans_new(
    name = "squish",
    transform = function(y) ifelse(
      y <= break_low, (y / break_low) * compress_frac,
      compress_frac + (y - break_low) / (y_top - break_low) * (1 - compress_frac)),
    inverse = function(z) ifelse(
      z <= compress_frac, z / compress_frac * break_low,
      break_low + (z - compress_frac) / (1 - compress_frac) * (y_top - break_low))
  )
  y_breaks <- c(0, 2, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50)
  y_breaks <- y_breaks[y_breaks <= y_top]

  p <- ggplot() +
    geom_scattermore(data = bg_data[chr_color == "odd"],  aes(BP_cum, log10p),
                     color = color_odd,  pointsize = 1.2, alpha = 0.6, pixels = c(2000, 800)) +
    geom_scattermore(data = bg_data[chr_color == "even"], aes(BP_cum, log10p),
                     color = color_even, pointsize = 1.2, alpha = 0.6, pixels = c(2000, 800)) +
    geom_point(data = suggestive_data, aes(BP_cum, log10p),
               color = color_suggestive, size = 2.5, alpha = 0.9) +
    geom_point(data = gw_data, aes(BP_cum, log10p),
               color = color_gw, size = 2.5, alpha = 0.9) +
    geom_hline(yintercept = -log10(threshold_suggestive),
               linetype = "dashed", color = color_suggestive, linewidth = 0.7) +
    geom_hline(yintercept = -log10(threshold_gw),
               linetype = "solid", color = color_gw, linewidth = 0.9) +
    geom_point(data = top10, aes(BP_cum, log10p),
               color = "black", shape = 1, size = 3, stroke = 0.6) +
    geom_text_repel(data = top10, aes(BP_cum, log10p, label = MarkerID),
                    size = 3.2, fontface = "bold",
                    min.segment.length = 0, box.padding = 0.5, point.padding = 0.3,
                    segment.color = "grey30", segment.size = 0.3,
                    max.overlaps = Inf, seed = 42) +
    scale_x_continuous(label = chr_info$CHR, breaks = chr_info$center, expand = c(0.01, 0.01)) +
    scale_y_continuous(trans = squish_trans, breaks = y_breaks, limits = c(0, y_top), expand = c(0, 0)) +
    labs(x = "Chromosome", y = expression(-log[10](italic(P))), title = tag) +
    theme_classic(base_size = 16) +
    theme(plot.title = element_text(face = "bold", size = 14, hjust = 0),
          axis.title = element_text(face = "bold", size = 18),
          axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          plot.margin = margin(15, 20, 10, 10))

  ggsave(file.path(out_dir, "Manhattan_plot.png"), p, width = 16, height = 5.5, dpi = 300, bg = "white")
  ggsave(file.path(out_dir, "Manhattan_plot.pdf"), p, width = 16, height = 5.5, device = cairo_pdf)
  rm(gm, bg_data, suggestive_data, gw_data, top10, p); gc(verbose = FALSE)
  invisible(list(n_gw = n_gw, n_sug = n_sug))
}

make_qq <- function(gwas, tag, out_dir, lambda_gc) {
  n <- nrow(gwas); sorted_p <- sort(gwas$p.value)
  qq_data <- data.table(expected = -log10((1:n) / (n + 1)), observed = -log10(sorted_p))
  ci_indices <- unique(c(seq(1, min(10000, n), by = 1),
                         seq(10001, min(100000, n), by = 10),
                         seq(100001, n, by = 100)))
  ci_indices <- ci_indices[ci_indices <= n]
  ci_lo <- -log10(qbeta(0.975, ci_indices, n - ci_indices + 1))
  ci_hi <- -log10(qbeta(0.025, ci_indices, n - ci_indices + 1))
  qq_data[, ci_lower := approx(ci_indices, ci_lo, 1:n, rule = 2)$y]
  qq_data[, ci_upper := approx(ci_indices, ci_hi, 1:n, rule = 2)$y]
  axis_max <- max(c(qq_data$expected, qq_data$observed), na.rm = TRUE) * 1.02

  p <- ggplot(qq_data, aes(expected, observed)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "grey80", alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1.2) +
    geom_scattermore(color = "#08519C", pointsize = 2, alpha = 0.6, pixels = c(1000, 1000)) +
    annotate("text", x = axis_max * 0.05, y = axis_max * 0.95,
             label = sprintf("lambda[GC] == %.3f", lambda_gc), parse = TRUE,
             hjust = 0, size = 6, fontface = "bold") +
    labs(x = expression(Expected ~ -log[10](italic(P))),
         y = expression(Observed ~ -log[10](italic(P))), title = tag) +
    scale_x_continuous(limits = c(0, axis_max), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, axis_max), expand = c(0, 0)) +
    coord_fixed(ratio = 1, xlim = c(0, axis_max), ylim = c(0, axis_max)) +
    theme_classic(base_size = 16) +
    theme(plot.title = element_text(face = "bold", size = 14, hjust = 0),
          axis.title = element_text(face = "bold", size = 18),
          axis.text = element_text(size = 14),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          plot.margin = margin(15, 15, 10, 10))

  ggsave(file.path(out_dir, "QQ_plot.png"), p, width = 8, height = 8, dpi = 300, bg = "white")
  ggsave(file.path(out_dir, "QQ_plot.pdf"), p, width = 8, height = 8, device = cairo_pdf)
  rm(qq_data, sorted_p, ci_lo, ci_hi, p); gc(verbose = FALSE)
}

# ---- main loop: one run at a time ------------------------------------------
all_summary <- list()
for (job in gwas_jobs) {
  tag <- job$tag; file <- job$file
  cat("\n", strrep("=", 60), "\n[", tag, "] ", file, "\n", strrep("=", 60), "\n", sep = "")

  if (!file.exists(file)) {
    cat("!! not found, skipping - run 03_merge_sumstats.sh for this cohort/model\n")
    next
  }

  out_dir <- file.path(out_root, tag)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  gwas <- read_gwas(file)
  cat("  variants:", format(nrow(gwas), big.mark = ","), "\n")

  lambda_all   <- calc_lambda(gwas$p.value)
  lambda_maf01 <- calc_lambda(gwas[AF_Allele2 > 0.01 & AF_Allele2 < 0.99, p.value])
  lambda_maf05 <- calc_lambda(gwas[AF_Allele2 > 0.05 & AF_Allele2 < 0.95, p.value])
  cat("  lambdaGC all/0.01/0.05:", round(lambda_all, 4), "/",
      round(lambda_maf01, 4), "/", round(lambda_maf05, 4), "\n")

  man_res <- make_manhattan(gwas, tag, out_dir)
  make_qq(gwas, tag, out_dir, lambda_all)

  all_summary[[tag]] <- data.frame(
    tag = tag, N_SNP = nrow(gwas),
    lambda_GC = round(lambda_all, 4),
    lambda_MAF01 = round(lambda_maf01, 4),
    lambda_MAF05 = round(lambda_maf05, 4),
    N_GW_sig = man_res$n_gw, N_suggestive = man_res$n_sug)

  rm(gwas, man_res); gc(verbose = FALSE)
  cat("  [done] ->", out_dir, "\n")
}

summary_all <- rbindlist(all_summary)
# Export and dx upload to RAP  (all_runs_summary.csv + per-run plots/top10 under output/gwas_plots/)
fwrite(summary_all, file.path(out_root, "all_runs_summary.csv"))
cat("\nAll runs complete. Summary ->", file.path(out_root, "all_runs_summary.csv"), "\n")
print(summary_all)
