# ===========================================================================
# Fig 2 - Multi-omics association summary
#
# Layout (one flat patchwork grid so every panel shares the same axis tracks):
#
#   row 1   A  GWAS Manhattan                                    full width
#   row 2   B  GWAS QQ  |  C  APOE genotype  |  D  Olink PWAS volcano
#   row 3   E  Olink GSEA enrichment                             full width
#   row 4   F  NMR MWAS forest, 57 metabolites x G1/G2/G3        full width
#
# eBioMedicine artwork rules applied here (info-for-authors p10-11 +
# artwork guidelines, April 2026):
#   - no titles inside the artwork; every panel heading lives in the legend
#   - panel identifiers are uppercase letters only, set well above body size
#     so they cannot be misread as the subtype labels A/B/C/D
#   - no box outline on any panel
#   - solid contrasting threshold lines instead of several dashed styles
#   - one type size throughout; 183 mm wide, single page, >=300 dpi
#
# Writes (all Fig 2 deliverables live in one dedicated folder):
#   output/figures/fig2/Figure2.{pdf,png,svg}
#   output/figures/fig2/Figure2_preview.png   (FIG2_PREVIEW=1)
#   output/figures/fig2/Figure2_panelE_gsea_plotdata.csv
# The PDF rasterises only panel A's point cloud (scattermore) to stay small
# while keeping every other element vector.
# ===========================================================================

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
})

# Dedicated Fig 2 folder: all active Fig 2 deliverables live here and nowhere
# else in figures_workspace.
OUT_DIR <- "output/figures/fig2"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
FIG_STEM <- "Figure2"   # submission-style basename ("Figure 2" in the manuscript)

GWAS_F       <- "input/assoc_results/gwas_primary_sumstats.tsv.gz"
GWAS_HITS_F  <- "input/assoc_results/gwas_primary_lead_snps.csv"
APOE_F       <- "input/assoc_results/apoe_diplotype_association.csv"
PWAS_F       <- "input/assoc_results/pwas_primary.csv"
OLINK_MAP_F  <- "input/reference/olink_panel_mapping.csv"
MWAS_F       <- "input/assoc_results/mwas_primary.csv"  # 20-covariate primary model
MANIFEST_F   <- "output/feature_manifest/master_feature_manifest_final.csv"
MWAS_G2_F    <- "input/assoc_results/mwas_group2_primary.csv"  # held-out G2
MWAS_G3_F    <- "input/assoc_results/mwas_group3_primary.csv"  # held-out G3
GSEA_F       <- "output/pwas_gsea/panel_restricted/tables/GSEA_all_significant.csv"

COL <- list(
  chr1 = "#2F5D7C",
  chr2 = "#9CB8C8",
  gws  = "#B22222",
  sugg = "#9AA0A6",
  risk = "#C43C39",
  prot = "#2B6CB0",
  ref  = "#2D2D2D",
  grey = "#8C8C8C",
  tas  = "#5A2CA0"
)

# --- Typography -------------------------------------------------------------
# eBM asks for a single type size across the whole figure. Two are used:
# PT_TITLE for axis titles and PT_TEXT for everything else. PT_DENSE is held
# back for the only two axes where PT_TEXT physically collides - the 22
# chromosome ticks in A and the 57 metabolite ticks in F.
PT_TITLE <- 7
PT_TEXT  <- 6
PT_DENSE <- 5.5
PT_TAG   <- 10          # panel letters, ~1.7x body text (eBM compound-figure rule)

gs <- function(pt) pt / .pt   # geom_text()/annotate() size, in mm, from points

NMR_CLASS_ORDER <- c(
  "VLDL subclasses", "LDL/IDL subclasses", "HDL subclasses",
  "Cholesterol", "Triglycerides", "Phospholipids", "Apolipoproteins",
  "Other lipids", "Particle size", "Fatty acids", "Amino acids",
  "Glycolysis metabolites", "Ketone bodies", "Inflammation", "Fluid balance"
)
NMR_CLASS_PAL <- c(
  "VLDL subclasses"        = "#7A4EA3",
  "LDL/IDL subclasses"     = "#D45A4C",
  "HDL subclasses"         = "#3778B8",
  "Cholesterol"            = "#D99A28",
  "Triglycerides"          = "#B85C2B",
  "Phospholipids"          = "#2F9C95",
  "Apolipoproteins"        = "#6F7782",
  "Other lipids"           = "#4E9F50",
  "Particle size"          = "#475569",
  "Fatty acids"            = "#E18636",
  "Amino acids"            = "#8E63B0",
  "Glycolysis metabolites" = "#2FAE9B",
  "Ketone bodies"          = "#C7A92B",
  "Inflammation"           = "#B44646",
  "Fluid balance"          = "#9AA0A6",
  "Other"                  = "#B8BDC6"
)

theme_fig2 <- function() {
  theme_classic(base_size = PT_TITLE, base_family = "Arial") +
    theme(
      axis.line          = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks         = element_line(linewidth = 0.25, colour = "black"),
      axis.text          = element_text(size = PT_TEXT, colour = "black"),
      axis.title         = element_text(size = PT_TITLE, colour = "black"),
      legend.title       = element_text(size = PT_TEXT, face = "bold"),
      legend.text        = element_text(size = PT_TEXT),
      legend.key.size    = unit(2.6, "mm"),
      legend.background  = element_blank(),
      legend.margin      = margin(0, 0, 0, 0),
      legend.box.spacing = unit(1, "mm"),
      # eBM: the artwork carries no titles - they belong at the head of the
      # figure legend. Blanking the elements keeps that true by construction.
      plot.title         = element_blank(),
      plot.subtitle      = element_blank(),
      plot.tag           = element_text(size = PT_TAG, face = "bold",
                                        family = "Arial"),
      panel.grid         = element_blank(),
      plot.margin        = margin(2, 3, 2, 3)
    )
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

  # The dense point cloud (~2.9M markers after thinning) is the only thing that
  # makes the vector PDF enormous. Rasterise just this layer with scattermore
  # (already a dependency; panel B uses it): the two-tone chromosome colouring
  # is preserved as a mapped aesthetic, everything else stays vector. Verified
  # that the discrete colour + scale_colour_manual survives rasterisation.
  man_cloud <- if (requireNamespace("scattermore", quietly = TRUE)) {
    scattermore::geom_scattermore(aes(colour = chr_group),
                                  pointsize = 2.1, pixels = c(3400, 780),
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
    labs(x = "Chromosome", y = expression(-log[10]~"(P)")) +
    theme_fig2() +
    theme(axis.text.x  = element_text(size = PT_DENSE, hjust = 0.5, vjust = 0.8),
          axis.ticks.x = element_line(linewidth = 0.18))
}

panel_b_qq <- function(gwas) {
  p <- gwas[!is.na(p.value) & p.value > 0 & p.value <= 1, p.value]
  p <- sort(p)
  n <- length(p)
  qq <- data.table(
    expected = -log10(seq_len(n) / (n + 1)),
    observed = -log10(p),
    i = seq_len(n)
  )

  # Smooth 95% null envelope: exact for the early tail, then interpolated.
  ci_idx <- unique(c(
    seq(1, min(10000, n), by = 1),
    seq(10001, min(100000, n), by = 10),
    seq(100001, n, by = 100)
  ))
  ci_idx <- ci_idx[ci_idx <= n]
  ci_lo <- -log10(qbeta(0.975, ci_idx, n - ci_idx + 1))
  ci_hi <- -log10(qbeta(0.025, ci_idx, n - ci_idx + 1))
  qq[, `:=`(
    ci_lo = approx(ci_idx, ci_lo, i, rule = 2)$y,
    ci_hi = approx(ci_idx, ci_hi, i, rule = 2)$y
  )]

  lambda <- median(qchisq(1 - p, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)
  lim <- max(qq$expected, qq$observed, na.rm = TRUE) * 1.02

  qq_points <- if (requireNamespace("scattermore", quietly = TRUE)) {
    scattermore::geom_scattermore(
      aes(expected, observed),
      pointsize = 3.2,
      pixels = c(900, 900),
      colour = "#B43A35"
    )
  } else {
    geom_point(size = 0.32, alpha = 0.55, colour = "#B43A35")
  }

  ggplot(qq, aes(expected, observed)) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
                fill = "#BFDDF4", alpha = 0.62, colour = NA) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.3, colour = "grey45") +
    qq_points +
    annotate("text", x = lim * 0.04, y = lim * 0.96,
             label = sprintf("lambda[GC] == %.2f", lambda),
             parse = TRUE, hjust = 0, vjust = 1,
             size = gs(PT_TEXT), family = "Arial") +
    coord_equal(xlim = c(0, lim), ylim = c(0, lim), expand = FALSE) +
    labs(x = expression(Expected~-log[10]~"(P)"),
         y = expression(Observed~-log[10]~"(P)")) +
    theme_fig2()   # no panel border: eBM asks for artwork without box outlines
}

# 2 significant figures, per eBM, except in the far tail where an exponent is
# clearer; NA (the reference row of a forest plot) carries no P value.
format_p <- function(p) {
  out <- rep(NA_character_, length(p))
  ok <- !is.na(p)
  out[ok] <- ifelse(
    p[ok] < 1e-4,
    sub("e-0", "e-", sprintf("%.1e", p[ok])),
    sprintf("%.2g", p[ok])
  )
  out
}

# Forest with a dedicated right-hand P column: the x axis stops at the data
# range and the text sits beyond it, so no axis line runs under the numbers.
forest_plot <- function(d, colour = COL$risk,
                        xlim = c(0.78, 1.75),
                        breaks = c(0.8, 1, 1.25, 1.5),
                        x_text = 1.95) {
  d <- copy(d)
  d <- d[!is.na(OR) & !is.na(CI_lo) & !is.na(CI_hi)]
  d[, label := factor(label, levels = rev(label))]
  d[, p_label := fifelse(is.na(pval), "ref", paste0("P=", format_p(pval)))]
  d[, OR_clip := pmin(pmax(OR, xlim[1]), xlim[2])]
  d[, CI_lo_clip := pmin(pmax(CI_lo, xlim[1]), xlim[2])]
  d[, CI_hi_clip := pmin(pmax(CI_hi, xlim[1]), xlim[2])]

  ggplot(d, aes(y = label)) +
    geom_vline(xintercept = 1, linewidth = 0.3, colour = "grey55") +
    geom_segment(aes(x = CI_lo_clip, xend = CI_hi_clip,
                     y = label, yend = label),
                 linewidth = 0.45, colour = colour) +
    geom_point(aes(x = OR_clip), size = 1.5, colour = colour) +
    geom_text(aes(x = x_text, label = p_label), hjust = 0,
              size = gs(PT_TEXT), family = "Arial", colour = "grey25") +
    scale_x_log10(breaks = breaks, labels = as.character(breaks)) +
    coord_cartesian(xlim = xlim, clip = "off") +
    labs(x = "Odds ratio (95% CI)", y = NULL) +
    theme_fig2() +
    theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          # room on the right for the P column drawn outside the panel
          plot.margin = margin(2, 30, 2, 2))
}

panel_c_apoe <- function() {
  stopifnot(file.exists(APOE_F))
  # encoding="UTF-8" matters: the Comparison column carries Greek epsilon, and
  # an unflagged string makes the gsub below silently miss on this machine.
  d <- fread(APOE_F, encoding = "UTF-8")
  setnames(d, c("Comparison", "CI_lower", "CI_upper", "P"),
           c("label", "CI_lo", "CI_hi", "pval"), skip_absent = TRUE)
  d <- d[label != "---"]
  d[, label := gsub("\u03b5", "e", label)]    # \u escape: encoding-safe source
  stopifnot(!any(grepl("[^ -~]", d$label)))   # every epsilon must have mapped

  # Short row labels; the comparison each one encodes is spelled out in the
  # figure legend, which keeps the panel readable at 60 mm.
  relabel <- c(
    "Per e4 allele (trend)"      = "Per e4 allele",
    "e4 carrier vs non-carrier"  = "e4 carrier",
    "e3/e3 (reference)"          = "e3/e3 (ref)"
  )
  d[label %in% names(relabel), label := relabel[label]]
  row_order <- c(
    "Per e4 allele", "e4 carrier",
    "e4/e4", "e3/e4", "e2/e4", "e2/e3", "e2/e2", "e3/e3 (ref)"
  )
  d[, ord := match(label, row_order)]
  stopifnot(!any(is.na(d$ord)))               # no row silently unordered
  setorder(d, ord)
  forest_plot(d[, .(label, OR, CI_lo, CI_hi, pval)], colour = "#8A5A2B")
}

compress_y <- function(y, cut, factor = 0.38) {
  ifelse(y <= cut, y * factor, cut * factor + (y - cut))
}

panel_d_olink <- function() {
  stopifnot(file.exists(PWAS_F))
  d <- fread(PWAS_F)
  setnames(d, c("protein", "or", "or_lower", "or_upper"),
           c("label", "OR", "CI_lo", "CI_hi"), skip_absent = TRUE)
  if (file.exists(OLINK_MAP_F)) {
    map <- fread(OLINK_MAP_F)
    if (all(c("assay", "panel") %in% names(map))) {
      d <- merge(d, map[, .(label = assay, panel)], by = "label", all.x = TRUE)
    }
  }
  d[is.na(panel), panel := "Unmapped"]
  d[, neglog10p := -log10(pmax(pval, 1e-300))]
  d[, signed_effect := beta]
  d[, sig_class := fcase(
    sig_bonf %in% TRUE, "Bonferroni",
    sig_fdr %in% TRUE, "FDR",
    default = "Nominal"
  )]
  d[, sig_class := factor(sig_class, levels = c("Nominal", "FDR", "Bonferroni"))]
  # colour encodes DIRECTION (matching panel E), not the significance tier;
  # significance is carried by the threshold lines and by point size/alpha
  d[, dir_class := fcase(
    sig_class == "Nominal", "Nominal",
    signed_effect >= 0, "Risk",
    default = "Protective")]
  d[, dir_class := factor(dir_class, levels = c("Nominal", "Protective", "Risk"))]

  x_lim <- max(abs(d$signed_effect), na.rm = TRUE) * 1.10
  bonf_y <- -log10(0.05 / nrow(d))
  fdr_y <- if (any(d$sig_fdr %in% TRUE)) {
    -log10(max(d[sig_fdr %in% TRUE]$pval, na.rm = TRUE))
  } else NA_real_
  cut_y <- if (!is.na(fdr_y)) fdr_y else bonf_y
  d[, y_plot := compress_y(neglog10p, cut_y)]
  lab <- d[sig_fdr %in% TRUE][order(pval)]
  lab[, y_plot := compress_y(neglog10p, cut_y)]

  raw_breaks <- unique(c(0, 1, 2, round(cut_y, 1),
                         pretty(d$neglog10p[d$neglog10p > cut_y], n = 4)))
  raw_breaks <- raw_breaks[raw_breaks >= 0 & raw_breaks <= max(d$neglog10p, na.rm = TRUE)]
  plot_breaks <- compress_y(raw_breaks, cut_y)

  ggplot(d, aes(signed_effect, y_plot)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf,
             ymax = compress_y(cut_y, cut_y),
             fill = "grey92", alpha = 0.55) +
    geom_hline(yintercept = compress_y(bonf_y, cut_y),
               colour = COL$gws, linewidth = 0.35) +
    {if (!is.na(fdr_y))
      geom_hline(yintercept = compress_y(fdr_y, cut_y),
                 colour = COL$sugg, linewidth = 0.3) else NULL} +
    geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.25) +
    geom_point(aes(colour = dir_class, size = sig_class, alpha = sig_class)) +
    # Threshold labels live in a reserved strip along the left edge; the repel
    # layer below is fenced out of that strip so protein names cannot drift in.
    annotate("text", x = -x_lim * 0.97, y = compress_y(bonf_y, cut_y) + 0.10,
             label = "Bonferroni", colour = COL$gws,
             hjust = 0, vjust = 0, size = gs(PT_TEXT), family = "Arial") +
    {if (!is.na(fdr_y))
      annotate("text", x = -x_lim * 0.97, y = compress_y(fdr_y, cut_y) + 0.08,
               label = "FDR", colour = "grey35",
               hjust = 0, vjust = 0, size = gs(PT_TEXT), family = "Arial") else NULL} +
    annotate("text", x = x_lim * 0.97, y = compress_y(cut_y * 0.45, cut_y),
             label = "Nominal", colour = "grey45",
             hjust = 1, vjust = 0.5, size = gs(PT_TEXT), family = "Arial") +
    geom_text_repel(
      data = lab,
      aes(label = label),
      size = gs(PT_TEXT), family = "Arial",
      box.padding = 0.30, point.padding = 0.08,
      segment.size = 0.18, min.segment.length = 0,
      xlim = c(-x_lim * 0.62, NA),   # keep clear of the threshold-label strip
      seed = 42, max.overlaps = Inf
    ) +
    scale_colour_manual(
      values = c(Nominal = "#A9A9A9", Protective = COL$prot, Risk = COL$risk),
      guide = "none"
    ) +
    scale_size_manual(values = c(Nominal = 0.75, FDR = 1.15, Bonferroni = 1.65),
                      guide = "none") +
    scale_alpha_manual(values = c(Nominal = 0.40, FDR = 0.82, Bonferroni = 0.95),
                       guide = "none") +
    scale_x_continuous(limits = c(-x_lim, x_lim),
                       expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(breaks = plot_breaks,
                       labels = round(raw_breaks, 1),
                       expand = expansion(mult = c(0, 0.10))) +
    labs(x = "Log-odds effect per SD NPX",
         y = expression(-log[10]~"(P)")) +
    theme_fig2() +
    theme(legend.position = "none")
}

clean_gsea_term <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("^GOBP_", "", x)
  x <- gsub("^REACTOME_", "", x)
  x <- gsub("^KEGG_", "", x)
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(tolower(x))
  x <- gsub("\\bIl\\b", "IL", x)
  x <- gsub("\\bTnf\\b", "TNF", x)
  x <- gsub("\\bNf Kb\\b", "NF-kB", x)
  x <- gsub("\\bDna\\b", "DNA", x)
  x <- gsub("\\bRna\\b", "RNA", x)
  x
}

wrap_label <- function(x, width = 30) {
  vapply(x, function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
}

shorten_gsea_term <- function(x) {
  repl <- c(
    "Gamma Carboxylation Hypusinylation Hydroxylation and Arylsulfatase Activation" =
      "\u03b3-carboxylation/hydroxylation/sulfatase activation",
    "Response to Elevated Platelet Cytosolic Ca2" =
      "Response to elevated platelet cytosolic Ca2+",
    "Complement and Coagulation Cascades" = "Complement and coagulation cascades",
    "Cell Adhesion Molecules Cams" = "Cell adhesion molecules (CAMs)",
    "Amino Acid Biosynthetic Process" = "Amino acid biosynthetic process",
    "Neuron Projection Guidance" = "Neuron projection guidance",
    "Axon Development" = "Axon development",
    "Innate Immune System" = "Innate immune system",
    "Xenobiotic Metabolism" = "Xenobiotic metabolism",
    "Apical Junction" = "Apical junction"
  )
  hit <- match(x, names(repl))
  ifelse(is.na(hit), x, repl[hit])
}

panel_e_gsea <- function() {
  stopifnot(file.exists(GSEA_F))
  gsea <- fread(GSEA_F)
  stopifnot(all(c("Database", "Pathway", "NES", "padj", "direction") %in% names(gsea)))

  gsea[, neg_log10_fdr := -log10(pmax(padj, .Machine$double.xmin))]
  gsea[, abs_nes := abs(NES)]
  gsea[, direction := ifelse(NES >= 0, "Activated", "Suppressed")]
  gsea[, term_label := shorten_gsea_term(clean_gsea_term(Pathway))]
  gsea[, source_label := fifelse(Database == "GO_BP", "GO BP", Database)]

  priority_db <- c("Hallmark", "Reactome", "KEGG", "GO BP", "Custom")
  gsea[, source_rank := match(source_label, priority_db)]
  gsea[is.na(source_rank), source_rank := length(priority_db) + 1]

  top_gsea <- gsea[order(padj, -abs_nes)][
    , .SD[1:min(.N, 3)], by = source_label
  ][order(source_rank, padj, -abs_nes)]
  top_gsea[, source_label := factor(source_label, levels = priority_db[priority_db %in% unique(source_label)])]
  top_gsea[, source_index := as.integer(source_label)]
  top_gsea[, y_pos := (.N:1) + (max(source_index) - source_index) * 0.7]
  top_gsea <- top_gsea[order(y_pos)]

  # --- panel geometry (data units; the panel spans the full figure width) ---
  X_BAND_LAB <- -3.30    # database name, left aligned
  X_TERM_R   <- -1.34    # pathway names, right aligned here
  X_AX_L     <- -1.18    # left edge of the NES field
  X_AX_R     <-  1.20
  NES_SCALE  <-  0.48    # NES units -> x units
  top_gsea[, x_plot := NES * NES_SCALE]

  band_df <- top_gsea[, .(
    ymin = min(y_pos) - 0.42,
    ymax = max(y_pos) + 0.42,
    ymid = mean(range(y_pos))
  ), by = source_label]

  # The panel is one wide field: pathway names on the left, NES on the right.
  # The x axis therefore belongs only under the NES part, so the axis rule is
  # drawn as a segment pinned to the panel edge while the ticks, tick labels
  # and title stay on ggplot's real axis (which cannot collide with them).
  X_LO <- X_BAND_LAB - 0.10
  nes_title_hjust <- (0 - X_LO) / (X_AX_R - X_LO)

  fwrite(top_gsea, file.path(OUT_DIR, paste0(FIG_STEM, "_panelE_gsea_plotdata.csv")))

  pal_db_bg <- c(
    "Hallmark" = "#C43C39",
    "Reactome" = "#2B6CB0",
    "KEGG" = "#2F855A",
    "GO BP" = "#805AD5",
    "Custom" = "#B7791F"
  )
  pal_dir <- c(Activated = COL$risk, Suppressed = COL$prot)

  ggplot(top_gsea, aes(x = x_plot, y = y_pos)) +
    geom_rect(
      data = band_df,
      aes(xmin = X_AX_L, xmax = X_AX_R, ymin = ymin, ymax = ymax, fill = source_label),
      inherit.aes = FALSE, alpha = 0.08, colour = NA
    ) +
    geom_text(data = band_df, aes(x = X_BAND_LAB, y = ymid, label = source_label),
              inherit.aes = FALSE, hjust = 0, fontface = "bold",
              size = gs(PT_TEXT), family = "Arial", colour = "grey20") +
    geom_text(aes(x = X_TERM_R, label = wrap_label(term_label, width = 52)),
              hjust = 1, size = gs(PT_TEXT), family = "Arial",
              lineheight = 0.9, colour = "grey20") +
    geom_segment(aes(x = X_AX_L - 0.11, xend = X_AX_L, yend = y_pos),
                 linewidth = 0.28, colour = "grey35", lineend = "round") +
    geom_vline(xintercept = X_AX_L, linewidth = 0.32, colour = "grey35") +
    geom_vline(xintercept = 0, linewidth = 0.32, colour = "grey65") +
    geom_segment(aes(x = 0, xend = x_plot, yend = y_pos, colour = direction),
                 linewidth = 0.82, alpha = 0.72, lineend = "round") +
    geom_point(aes(size = neg_log10_fdr, fill = direction),
               shape = 21, colour = "grey20", stroke = 0.22, alpha = 0.96) +
    annotate("segment", x = X_AX_L, xend = X_AX_R, y = -Inf, yend = -Inf,
             linewidth = 0.32, colour = "grey30") +
    scale_colour_manual(values = pal_dir, guide = "none") +
    scale_fill_manual(values = c(pal_dir, pal_db_bg),
                      breaks = c("Activated", "Suppressed"),
                      name = "Direction") +
    scale_size_continuous(name = expression(-log[10]~"(FDR)"),
                          range = c(1.3, 3.8), breaks = pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = NULL, labels = NULL,
                       expand = expansion(mult = c(0.03, 0.04))) +
    scale_x_continuous(limits = c(X_LO, X_AX_R),
                       breaks = (-2:2) * NES_SCALE, labels = -2:2,
                       expand = expansion(mult = c(0, 0))) +
    labs(x = "NES", y = NULL) +
    guides(fill = guide_legend(override.aes = list(size = 2.2), order = 1, nrow = 1),
           size = guide_legend(order = 2, nrow = 1)) +
    theme_fig2() +
    theme(
      legend.position  = "bottom",
      legend.box       = "horizontal",
      legend.spacing.x = unit(1.2, "mm"),
      # centre the axis title over the NES field, not over the whole panel
      axis.title.x = element_text(hjust = nes_title_hjust),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line    = element_blank(),
      plot.margin  = margin(2, 3, 1, 3)
    )
}

# NMR metabolome: 57 de-redundant metabolites x G1/G2/G3
# (Julkunen-style grouped point-range)
panel_f_nmr_forest <- function() {
  # The reported non-redundant NMR panel = the MWAS rows of the frozen feature
  # manifest (one representative per correlation cluster).
  man <- fread(MANIFEST_F)[source_analysis == "MWAS"]; set57 <- man$feature_id
  g1 <- fread(MWAS_F); g2 <- fread(MWAS_G2_F); g3 <- fread(MWAS_G3_F)
  pick <- function(d, coh) d[protein %in% set57, .(metabolite = protein, or, or_lower, or_upper, cohort = coh)]
  d <- rbindlist(list(pick(g1, "G1"), pick(g2, "G2"), pick(g3, "G3")))
  d[, metabolite := factor(metabolite, levels = sort(set57))]
  d[, cohort := factor(cohort, levels = c("G1", "G2", "G3"))]
  lmap <- fread("input/reference/nmr_label_map.csv")
  lut  <- setNames(lmap$label, lmap$code)   # official Nightingale abbreviation (else cleaned code)
  pal  <- c(G1 = "#E64B35", G2 = "#4DBBD5", G3 = "#3C5488")
  labs <- c(G1 = "G1 (discovery)", G2 = "G2 (other White)", G3 = "G3 (non-White)")
  ggplot(d, aes(metabolite, or, colour = cohort)) +
    geom_hline(yintercept = 1, linewidth = 0.4, colour = "black") +
    geom_hline(yintercept = c(0.8, 1.2), linewidth = 0.22, colour = "grey85") +
    geom_pointrange(aes(ymin = or_lower, ymax = or_upper),
                    position = position_dodge(width = 0.65),
                    size = 0.2, fatten = 1.3, linewidth = 0.36) +
    scale_colour_manual(values = pal, labels = labs, name = NULL) +
    scale_x_discrete(labels = function(x) lut[as.character(x)]) +
    scale_y_log10(breaks = c(0.6, 0.8, 1.0, 1.2, 1.5), limits = c(0.5, 1.9)) +
    labs(x = NULL, y = "Odds ratio (95% CI) per SD") +
    theme_fig2() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                     size = PT_DENSE),
          legend.position = "top",
          legend.direction = "horizontal",
          legend.justification = "left",
          legend.key.spacing.x = unit(2, "mm"))
}

# ---------------------------------------------------------------------------
# Assemble
# ---------------------------------------------------------------------------
message("Loading GWAS summary statistics")
gwas <- read_gwas()

pA <- panel_a_manhattan(gwas)   # GWAS Manhattan
pB <- panel_b_qq(gwas)          # GWAS QQ
pC <- panel_c_apoe()            # APOE genotype forest
pD <- panel_d_olink()           # Olink PWAS volcano
pE <- panel_e_gsea()            # Olink GSEA enrichment
pF <- panel_f_nmr_forest()      # NMR MWAS forest, G1/G2/G3

# One flat design rather than nested "/" blocks: with a single grid, patchwork
# puts every left-hand y axis in the same track, so the panel edges of A, B, E
# and F line up exactly. Columns are 1/100 of the figure width.
design <- c(
  area(1,  1, 1, 100),   # A
  area(2,  1, 2,  28),   # B
  area(2, 29, 2,  61),   # C
  area(2, 62, 2, 100),   # D
  area(3,  1, 3, 100),   # E
  area(4,  1, 4, 100)    # F
)

fig2 <- pA + pB + pC + pD + pE + pF +
  # row 2 is sized so that B's coord_equal square exactly fills its 28 columns;
  # shrinking it would letterbox the QQ panel
  plot_layout(design = design, heights = c(50, 52, 40, 66)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = PT_TAG, face = "bold", family = "Arial"))

FIG_W <- 183   # mm, eBM double-column
FIG_H <- 205   # mm, leaves room for a 10 pt legend on the same page

if (identical(Sys.getenv("FIG2_PREVIEW"), "1")) {
  ggsave(file.path(OUT_DIR, paste0(FIG_STEM, "_preview.png")), fig2,
         width = FIG_W, height = FIG_H, units = "mm", dpi = 220, bg = "white")
} else {
  pdf_path <- file.path(OUT_DIR, paste0(FIG_STEM, ".pdf"))
  ggsave(pdf_path, fig2,
         width = FIG_W, height = FIG_H, units = "mm", device = cairo_pdf)
  ggsave(file.path(OUT_DIR, paste0(FIG_STEM, ".png")), fig2,
         width = FIG_W, height = FIG_H, units = "mm", dpi = 600, bg = "white")
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(file.path(OUT_DIR, paste0(FIG_STEM, ".svg")), fig2,
           width = FIG_W, height = FIG_H, units = "mm", bg = "white")
  } else {
    message("Skipping SVG export because package 'svglite' is not installed.")
  }
  message(sprintf("PDF size: %.1f MB", file.info(pdf_path)$size / 1024^2))
}

message("Wrote ", FIG_STEM, " to ", OUT_DIR)
