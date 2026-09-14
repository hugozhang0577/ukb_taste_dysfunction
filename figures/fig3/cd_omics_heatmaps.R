#> in   output/subtyping/evidence_omics/omics_signature_nmr_plotdata.csv
#> in   output/subtyping/evidence_omics/omics_signature_olink_plotdata.csv
#> out  output/figures/fig3c_olink_heatmap.pdf
#> out  output/figures/fig3d_nmr_heatmap.pdf
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
source("figures/_common/style.R")
source("figures/_common/fig3_type.R")   # local 10% type bump
suppressPackageStartupMessages(library(scales))

ev <- "output/subtyping"
DIV_LOW <- "#2166AC"; DIV_HI <- "#B2182B"

## ---- data --------------------------------------------------------------------
nmr <- fread(file.path(ev, "evidence_omics/omics_signature_nmr_plotdata.csv"))
oli <- fread(file.path(ev, "evidence_omics/omics_signature_olink_plotdata.csv"))
cat(sprintf("[load] NMR rows=%d  Olink rows=%d\n", nrow(nmr), nrow(oli)))
agg <- function(d) d[, .(beta = mean(beta), psig = 100 * mean(q < 0.05)),
                     by = .(platform, group, subtype)]
om <- rbind(agg(nmr), agg(oli))

nmr_ord <- c("Cholesterol (totals)", "Triglycerides (totals)", "Apolipoproteins & particle size",
             "VLDL subclasses", "IDL subclasses", "LDL subclasses", "HDL subclasses",
             "Phospholipids & other lipids", "Fatty acids", "Amino acids", "Glycolysis")
NMR_SHORT <- c("Apolipoproteins & particle size" = "Apolipoproteins & size",
               "Phospholipids & other lipids"    = "Phospholipids & other",
               "Glycolysis, ketones & other"     = "Glycolysis & ketones")
oli_ord <- c("Cardiometabolic", "Cardiometabolic_II", "Inflammation", "Inflammation_II",
             "Neurology", "Neurology_II", "Oncology", "Oncology_II")
ord_of <- function(plat, curated) {
  have <- unique(as.character(om[platform == plat, group]))
  c(intersect(curated, have), setdiff(have, curated))
}
nmr_ord <- ord_of("NMR",   nmr_ord)
oli_ord <- ord_of("Olink", oli_ord)
stopifnot(nrow(om) == 4 * (length(nmr_ord) + length(oli_ord)))
cat(sprintf("[groups] NMR %d | Olink %d | cells %d\n",
            length(nmr_ord), length(oli_ord), nrow(om)))

## ---- one panel per platform ------------------------------------------------
CELL_W <- round(5.0*FIG3_K,2); CELL_H <- round(3.0*FIG3_K,2)     # mm; MUST match fig3e_disease_heatmap.R
OVERHEAD <- round(20.0*FIG3_K,1)                 # header + bottom legend + margin (mm), measured
W_MM <- 52.7                     # (168 - 2 gutters of 5) / 3
SZA  <- FIG3_SZ$anno / .pt

fix_cells <- function(p, n_rows, ncol, cw, ch) {
  g <- ggplotGrob(p)
  pan <- g$layout[grepl("^panel", g$layout$name), , drop = FALSE]
  for (l in unique(pan$l)) g$widths[l]  <- grid::unit(cw * ncol, "mm")
  for (tt in unique(pan$t)) g$heights[tt] <- grid::unit(ch * n_rows, "mm")
  g
}

mk <- function(plat, ord, lim, brk, legname) {
  d <- om[platform == plat]
  d[, group := factor(group, levels = rev(ord))]
  stopifnot(!anyNA(d$group))
  ggplot(d, aes(subtype, group)) +
    geom_tile(aes(fill = beta), colour = NA) +
    geom_point(aes(size = psig), shape = 21, colour = "grey15", fill = NA, stroke = .35) +
    scale_fill_gradient2(low = DIV_LOW, mid = "white", high = DIV_HI, midpoint = 0,
                         limits = lim, breaks = brk, oob = squish, name = legname,
                         guide = guide_colourbar(title.position = "top", order = 1,
                                                 barwidth = unit(9, "mm"),    # at this type size 11 mm pushes the 100 of %FDR-sig off the right edge
                                                 barheight = unit(1.8, "mm"))) +
    scale_size_area("% FDR", max_size = 1.1, breaks = c(25, 100), limits = c(0, 100),
                    guide = guide_legend(title.position = "top", order = 2, nrow = 2,
                                         override.aes = list(stroke = .5))) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    # the Olink panel names carry an internal `_II` suffix; render with a space
    scale_y_discrete(expand = c(0, 0), labels = function(x) {
      x <- ifelse(x %in% names(NMR_SHORT), NMR_SHORT[x], x)
      gsub("_", " ", x) }) +
    labs(x = NULL, y = NULL) +
    theme_fig(grid = "none") +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
          axis.text.x.top = element_text(size = FIG3_SZ$anno, face = "bold", colour = "black"),
          axis.text.y = element_text(size = FIG3_SZ$axtxt, colour = "black"),
          legend.position = "bottom", legend.box = "horizontal",
          legend.justification = "left", legend.box.just = "left",
          legend.key.width = unit(2.2, "mm"), legend.key.height = unit(2.2, "mm"),
          legend.margin = margin(1, 0, 0, 0, "mm"),
          legend.box.spacing = unit(1, "mm"),
          legend.spacing.x = unit(1, "mm"),
          plot.margin = margin(1, 1.5, 1, 1, "mm"))
}

pOli <- mk("Olink", oli_ord, c(-.3, .3), c(-.3, 0, .3), "Olink")
pNmr <- mk("NMR",   nmr_ord, c(-.8, .8), c(-.8, 0, .8), "NMR")

gOli <- fix_cells(pOli, length(oli_ord), 4, CELL_W, CELL_H)
gNmr <- fix_cells(pNmr, length(nmr_ord), 4, CELL_W, CELL_H)
h_of <- function(n) round(n * CELL_H + OVERHEAD, 1)
cat(sprintf("[height] Olink %d rows -> %.1f mm | NMR %d rows -> %.1f mm\n",
            length(oli_ord), h_of(length(oli_ord)), length(nmr_ord), h_of(length(nmr_ord))))
fig_save(gOli, "fig3c_olink_heatmap", W_MM, h_of(length(oli_ord)))
fig_save(gNmr, "fig3d_nmr_heatmap",   W_MM, h_of(length(nmr_ord)))

cat("\n[legend] Cell fill is the subtype mean adjusted effect against the other three\n")
cat("         pooled; circle area is the percentage of assays in that group significant\n")
cat("         at 5% FDR. NOTE the two panels use DIFFERENT fill limits and units:\n")
cat("         NMR betas are in SD, Olink betas in native NPX (log2 relative), and\n")
cat("         the ranges are not comparable between panels.\n")
