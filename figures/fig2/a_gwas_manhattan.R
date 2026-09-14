#> in   input/assoc_results/gwas_primary_lead_snps.csv
#> in   input/assoc_results/gwas_primary_sumstats.tsv.gz
#> out  output/figures/fig2a_gwas_manhattan.pdf
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel)
})
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
source("figures/_common/style.R")

GWAS_F      <- "input/assoc_results/gwas_primary_sumstats.tsv.gz"
GWAS_HITS_F <- "input/assoc_results/gwas_primary_lead_snps.csv"

COL <- list(chr1 = "#2F5D7C", chr2 = "#9CB8C8", gws = "#B22222",
            sugg = "#9AA0A6", tas = "#5A2CA0")

PT_TITLE <- 7; PT_TEXT <- 6; PT_DENSE <- 5.5
gs <- function(pt) pt / .pt

theme_man <- function() {
  theme_classic(base_size = PT_TITLE, base_family = FIG_FONT) +
    theme(axis.line   = element_line(linewidth = 0.3, colour = "black"),
          axis.ticks  = element_line(linewidth = 0.25, colour = "black"),
          axis.text   = element_text(size = PT_TEXT, colour = "black"),
          axis.title  = element_text(size = PT_TITLE, colour = "black"),
          plot.title  = element_blank(), plot.subtitle = element_blank(),
          panel.grid  = element_blank(), plot.margin = margin(2, 3, 2, 3))
}

read_gwas <- function() {
  stopifnot(file.exists(GWAS_F))
  d <- fread(GWAS_F, select = c("CHR", "POS", "MarkerID", "BETA", "SE", "p.value"))
  d <- d[CHR %in% 1:22 & !is.na(POS) & !is.na(p.value) & p.value > 0]
  d[, neglog10p := -log10(pmax(p.value, 1e-300))]
  setorder(d, CHR, POS)
  d
}

add_cumulative_bp <- function(d) {
  chr_lens <- d[, .(maxbp = max(POS, na.rm = TRUE)), by = CHR][order(CHR)]
  chr_lens[, offset := c(0, cumsum(as.numeric(maxbp))[-.N])]
  out <- merge(d, chr_lens[, .(CHR, offset)], by = "CHR", all.x = TRUE)
  out[, BP_cum := POS + offset]
  axis_dt <- chr_lens[, .(
    CHR,
    center = offset + maxbp / 2
  )]
  list(data = out, axis = axis_dt, chr_lens = chr_lens)
}

thin_manhattan <- function(d, bin_size = 100000, quantile_keep = 0.72) {
  d[, bin := floor(POS / bin_size)]
  d[, keep_q := quantile(neglog10p, quantile_keep, na.rm = TRUE), by = .(CHR, bin)]
  out <- d[neglog10p >= keep_q | p.value < 1e-5]
  out[, c("bin", "keep_q") := NULL]
  out
}

panel_a_manhattan <- function(gwas) {
  x <- add_cumulative_bp(thin_manhattan(copy(gwas)))
  d <- x$data
  axis_dt <- x$axis
  chr_lens <- x$chr_lens
  d[, chr_group := factor(CHR %% 2)]

  hits <- if (file.exists(GWAS_HITS_F)) {
    h <- fread(GWAS_HITS_F)
    h <- h[p.value < 5e-8 & CHR %in% 1:22]
    h[, neglog10p := -log10(p.value)]
    h <- merge(h, unique(d[, .(CHR, offset)]), by = "CHR", all.x = TRUE)
    h[, BP_cum := POS + offset]
    h[order(p.value)][1:min(.N, 8)]
  } else data.table()

  tas_regions <- data.table(
    locus = c("TAS1R", "TAS2R chr7", "TAS2R chr12"),
    CHR = c(1, 7, 12),
    lo = c(1000000, 141500000, 10500000),
    hi = c(20000000, 142500000, 11500000)
  )
  tas_regions <- merge(tas_regions, chr_lens[, .(CHR, offset)], by = "CHR")
  tas_regions[, `:=`(xmin = offset + lo, xmax = offset + hi)]
  tas_top <- rbindlist(lapply(seq_len(nrow(tas_regions)), function(i) {
    r <- tas_regions[i]
    top <- gwas[CHR == r$CHR & POS >= r$lo & POS <= r$hi][order(p.value)][1]
    if (!nrow(top)) return(data.table())
    top[, `:=`(
      BP_cum = POS + r$offset,
      neglog10p = -log10(pmax(p.value, 1e-300)),
      label = sprintf("%s top\n%s", r$locus, MarkerID)
    )]
    top
  }), fill = TRUE)

  man_cloud <- if (requireNamespace("scattermore", quietly = TRUE)) {
    scattermore::geom_scattermore(aes(colour = chr_group),
                                  pointsize = 2.8, pixels = c(5200, 1150),
                                  alpha = 0.70)
  } else {
    geom_point(aes(colour = chr_group), size = 0.32, alpha = 0.70)
  }

  ggplot(d, aes(BP_cum, neglog10p)) +
    geom_rect(data = tas_regions,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = "#7E57C2", alpha = 0.10) +
    # eBM prefers solid contrasting rules over several dashed line styles
    geom_hline(yintercept = -log10(5e-8), colour = COL$gws, linewidth = 0.35) +
    geom_hline(yintercept = -log10(1e-5), colour = COL$sugg, linewidth = 0.3) +
    man_cloud +
    {if (nrow(hits))
      geom_point(data = hits, aes(BP_cum, neglog10p),
                 inherit.aes = FALSE, colour = COL$gws,
                 size = 1.35, alpha = 0.95) else NULL} +
    {if (nrow(hits))
      geom_text_repel(data = hits,
                      aes(BP_cum, neglog10p, label = MarkerID),
                      inherit.aes = FALSE, size = gs(PT_TEXT), family = "Arial",
                      box.padding = 0.25, point.padding = 0.1,
                      segment.size = 0.18, min.segment.length = 0,
                      max.overlaps = Inf) else NULL} +
    {if (nrow(tas_top))
      geom_point(data = tas_top, aes(BP_cum, neglog10p),
                 inherit.aes = FALSE, colour = COL$tas,
                 size = 1.25, shape = 23, fill = "white", stroke = 0.45) else NULL} +
    {if (nrow(tas_top))
      geom_text_repel(data = tas_top,
                      aes(BP_cum, neglog10p, label = label),
                      inherit.aes = FALSE, size = gs(PT_TEXT), family = "Arial",
                      colour = COL$tas, fontface = "bold",
                      lineheight = 0.9,
                      nudge_y = 1.4, ylim = c(5.2, 7.5), direction = "both",
                      force = 8, force_pull = 0.15,
                      box.padding = 0.5, point.padding = 0.2,
                      segment.size = 0.18, segment.colour = COL$tas,
                      min.segment.length = 0, seed = 42,
                      max.overlaps = Inf) else NULL} +
    scale_colour_manual(values = c("0" = COL$chr1, "1" = COL$chr2), guide = "none") +
    scale_x_continuous(breaks = axis_dt$center, labels = axis_dt$CHR,
                       expand = expansion(mult = c(0.004, 0.01))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = "Chromosome", y = "-log10 (P)") +
    theme_man() +
    theme(axis.text.x  = element_text(size = PT_DENSE, hjust = 0.5, vjust = 0.8),
          axis.ticks.x = element_line(linewidth = 0.18))
}

# Panel letters are added at assembly, so the tag slot is cleared here.
strip_tag <- function(p) p + labs(tag = NULL) +
  theme(plot.title = element_blank(), plot.subtitle = element_blank(),
        plot.tag = element_blank(), text = element_text(family = FIG_FONT))

pA <- panel_a_manhattan(read_gwas())
fig_save(strip_tag(pA), "fig2a_gwas_manhattan", 168, 40)
