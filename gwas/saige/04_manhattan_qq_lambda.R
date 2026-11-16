# Primary GWAS summary: genomic-control lambda (all variants, MAF > 0.01, MAF > 0.05), counts at P < 5e-8
# and P < 1e-5, the 10 lowest-P variants, and Manhattan and QQ plots.
library(data.table)
library(ggplot2)
library(scattermore)
library(ggrepel)
library(scales)
setwd(Sys.getenv("PROJECT_DIR", "."))
TAG     <- "primary"
OUT_DIR <- file.path("output/gwas_plots", TAG)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
P_GW <- 5e-8; P_SUG <- 1e-5

gwas <- fread("input/assoc_results/gwas_primary_sumstats.tsv.gz",
              select = c("CHR", "POS", "MarkerID", "Allele1", "Allele2", "AF_Allele2", "BETA", "SE", "p.value"))
gwas <- gwas[!is.na(p.value) & p.value > 0 & p.value < 1]
gwas[, CHR := as.integer(CHR)]
calc_lambda <- function(p) median(qchisq(1 - p, 1), na.rm = TRUE) / qchisq(0.5, 1)

# ---- Manhattan --------------------------------------------------------------
gm <- gwas[!is.na(CHR) & CHR %in% 1:22][order(CHR, POS)]
chr_info <- gm[, .(chr_len = max(POS)), by = CHR][order(CHR)]
chr_info[, cumsum_len := cumsum(as.numeric(chr_len)) - chr_len]
chr_info[, center := cumsum_len + chr_len / 2]
gm <- merge(gm, chr_info[, .(CHR, cumsum_len)], by = "CHR", sort = FALSE)
gm[, `:=`(BP_cum = POS + cumsum_len, log10p = -log10(p.value))]
gm[, snp_category := fifelse(p.value < P_GW, "genome_wide", fifelse(p.value < P_SUG, "suggestive", "background"))]
bg <- gm[snp_category == "background"][, chr_color := fifelse(CHR %% 2 == 1, "odd", "even")]

top10 <- gm[order(p.value)][1:min(10, nrow(gm))]
fwrite(top10[, .(MarkerID, CHR, POS, Allele1, Allele2, AF_Allele2, BETA, SE, p.value)], file.path(OUT_DIR, "top10_snps.csv"))

# piecewise y-axis: the band below P = 1e-5 takes 30% of the height
break_low <- -log10(P_SUG)
y_top <- max(max(gm$log10p) * 1.08, break_low + 2)
squish_trans <- trans_new("squish",
  transform = function(y) ifelse(y <= break_low, y / break_low * 0.3, 0.3 + (y - break_low) / (y_top - break_low) * 0.7),
  inverse   = function(z) ifelse(z <= 0.3, z / 0.3 * break_low, break_low + (z - 0.3) / 0.7 * (y_top - break_low)))
y_breaks <- c(0, 2, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50)
y_breaks <- y_breaks[y_breaks <= y_top]

p <- ggplot() +
  geom_scattermore(data = bg[chr_color == "odd"],  aes(BP_cum, log10p), color = "#6BAED6", pointsize = 1.2, alpha = 0.6, pixels = c(2000, 800)) +
  geom_scattermore(data = bg[chr_color == "even"], aes(BP_cum, log10p), color = "#08519C", pointsize = 1.2, alpha = 0.6, pixels = c(2000, 800)) +
  geom_point(data = gm[snp_category == "suggestive"],  aes(BP_cum, log10p), color = "#FF8C00", size = 2.5, alpha = 0.9) +
  geom_point(data = gm[snp_category == "genome_wide"], aes(BP_cum, log10p), color = "#E31A1C", size = 2.5, alpha = 0.9) +
  geom_hline(yintercept = -log10(P_SUG), linetype = "dashed", color = "#FF8C00", linewidth = 0.7) +
  geom_hline(yintercept = -log10(P_GW), linetype = "solid", color = "#E31A1C", linewidth = 0.9) +
  geom_point(data = top10, aes(BP_cum, log10p), color = "black", shape = 1, size = 3, stroke = 0.6) +
  geom_text_repel(data = top10, aes(BP_cum, log10p, label = MarkerID), size = 3.2, fontface = "bold",
                  min.segment.length = 0, box.padding = 0.5, point.padding = 0.3, segment.color = "grey30",
                  segment.size = 0.3, max.overlaps = Inf, seed = 42) +
  scale_x_continuous(label = chr_info$CHR, breaks = chr_info$center, expand = c(0.01, 0.01)) +
  scale_y_continuous(trans = squish_trans, breaks = y_breaks, limits = c(0, y_top), expand = c(0, 0)) +
  labs(x = "Chromosome", y = expression(-log[10](italic(P))), title = TAG) +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0), axis.title = element_text(face = "bold", size = 18),
        axis.text = element_text(size = 13), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.margin = margin(15, 20, 10, 10))
ggsave(file.path(OUT_DIR, "Manhattan_plot.png"), p, width = 16, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "Manhattan_plot.pdf"), p, width = 16, height = 5.5, device = cairo_pdf)

# ---- QQ with 95% pointwise envelope -----------------------------------------
lambda_all <- calc_lambda(gwas$p.value)
n <- nrow(gwas)
qq <- data.table(expected = -log10((1:n) / (n + 1)), observed = -log10(sort(gwas$p.value)))
idx <- unique(c(seq(1, min(10000, n)), seq(10001, min(100000, n), by = 10), seq(100001, n, by = 100)))
idx <- idx[idx <= n]
qq[, ci_lower := approx(idx, -log10(qbeta(0.975, idx, n - idx + 1)), 1:n, rule = 2)$y]
qq[, ci_upper := approx(idx, -log10(qbeta(0.025, idx, n - idx + 1)), 1:n, rule = 2)$y]
axis_max <- max(c(qq$expected, qq$observed), na.rm = TRUE) * 1.02
p <- ggplot(qq, aes(expected, observed)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "grey80", alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1.2) +
  geom_scattermore(color = "#08519C", pointsize = 2, alpha = 0.6, pixels = c(1000, 1000)) +
  annotate("text", x = axis_max * 0.05, y = axis_max * 0.95, label = sprintf("lambda[GC] == %.3f", lambda_all),
           parse = TRUE, hjust = 0, size = 6, fontface = "bold") +
  labs(x = expression(Expected ~ -log[10](italic(P))), y = expression(Observed ~ -log[10](italic(P))), title = TAG) +
  scale_x_continuous(limits = c(0, axis_max), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, axis_max), expand = c(0, 0)) +
  coord_fixed(ratio = 1, xlim = c(0, axis_max), ylim = c(0, axis_max)) +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0), axis.title = element_text(face = "bold", size = 18),
        axis.text = element_text(size = 14), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.margin = margin(15, 15, 10, 10))
ggsave(file.path(OUT_DIR, "QQ_plot.png"), p, width = 8, height = 8, dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "QQ_plot.pdf"), p, width = 8, height = 8, device = cairo_pdf)

fwrite(data.table(tag = TAG, N_SNP = n, lambda_GC = round(lambda_all, 4),
                  lambda_MAF01 = round(calc_lambda(gwas[AF_Allele2 > 0.01 & AF_Allele2 < 0.99, p.value]), 4),
                  lambda_MAF05 = round(calc_lambda(gwas[AF_Allele2 > 0.05 & AF_Allele2 < 0.95, p.value]), 4),
                  N_GW_sig = sum(gm$snp_category == "genome_wide"), N_suggestive = sum(gm$snp_category == "suggestive")),
       "output/gwas_plots/all_runs_summary.csv")
