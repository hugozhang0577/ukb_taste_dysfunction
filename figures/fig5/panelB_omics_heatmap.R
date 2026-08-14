# =============================================================================
# Fig 5 Panel B  --  full-panel omics (NMR class x Olink panel) x subtype
# Rendered as a standalone image and composited by fig5/assemble.R.
# Theme, colours, geoms, scales and strip styling are shared with the other panels.
# ASCII glyphs only; cairo_pdf safe.  data.table only.
# =============================================================================
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(scales); library(grid)
                                 library(ggnewscale) })

ev  <- "output/subtyping"
out <- "output/figures/fig5/panels"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ---- shared aesthetics (identical across the fig5 panels) -------------------
LAYER   <- c(L1 = "#4E9080", L2a = "#2C6FAC", L2b = "#D4813A", L3 = "#9B59B6", L4 = "#C0392B")
DIV_LOW <- "#2166AC"; DIV_HI <- "#B2182B"

# Fonts are scale-compensated: this panel is composited at ~0.64x, so source pt are
# set larger so final sizes match the other panels (see fig5_master_assemble.R).
base_thm <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title    = element_text(face = "bold", size = 20, colour = "black"),
        axis.text     = element_text(colour = "black"),
        legend.title  = element_text(size = 13), legend.text = element_text(size = 12),
        legend.key.height = unit(4, "mm"), legend.key.width = unit(4, "mm"))

# ---- data --------------------------------------------------------------------
# Age-, sex- and centre-adjusted one-versus-rest effects, restricted to the
# biomarkers reaching FDR in at least one subtype (subtype_characterisation/
# 01_fullpanel_omics.R).
nmr <- fread(file.path(ev, "evidence_omics/omics_signature_nmr_plotdata.csv"))
oli <- fread(file.path(ev, "evidence_omics/omics_signature_olink_plotdata.csv"))
cat(sprintf("NMR rows=%d  Olink rows=%d\n", nrow(nmr), nrow(oli)))
agg <- function(d) d[, .(beta = mean(beta), psig = 100*mean(q < 0.05)),
                     by = .(platform, group, subtype)]
om <- rbind(agg(nmr), agg(oli))
nmr_ord <- c("Cholesterol (totals)","Triglycerides (totals)","Apolipoproteins & particle size",
             "VLDL subclasses","IDL subclasses","LDL subclasses","HDL subclasses",
             "Phospholipids & other lipids","Fatty acids","Amino acids",
             "Glycolysis, ketones & other")   # label as nmr_class() emits it
oli_ord <- c("Cardiometabolic","Cardiometabolic_II","Inflammation","Inflammation_II",
             "Neurology","Neurology_II","Oncology","Oncology_II")
om[, platform := factor(platform, c("NMR","Olink"))]
known <- c(nmr_ord, oli_ord)
extra <- setdiff(unique(as.character(om$group)), known)
om[, group := factor(group, levels = rev(c(known, extra)))]
n0 <- nrow(om); om <- om[!is.na(group)]
cat(sprintf("cells: %d -> %d after group-factor\n", n0, nrow(om)))

# ---- panel b ------------------------------------------------------------------
# The two platforms get SEPARATE fill scales, because their betas are in different
# units and one shared scale cannot serve both:
#   NMR   : the metabolite matrix is z-scored (all 327 measures SD 0.91-1.04), so
#           betas are in SD units. Group means span -0.743..+0.819.
#   Olink : the protein matrix is centred but NOT scaled (2,920 assays, SD
#           0.082-3.285, only 7% near 1), so betas are in native NPX, i.e. log2
#           relative units. Group means span only -0.213..+0.189.
# A single limit wide enough for the metabolome would flatten the proteome, and one
# narrow enough for the proteome would saturate the metabolome.
# NPX is normalised (same protein comparable across samples/plates) but not
# standardised (different proteins are not on a common variance scale), so an NPX
# beta reads as a log2 fold-change, not as an effect size in SD. Keeping the raw
# NPX scale is deliberate: Supplementary Results 6.5 compares the cis-pQTL CI to
# the subtype A-vs-C protein difference, and that comparison requires both in NPX.
pB <- ggplot(mapping = aes(subtype, group)) +
  geom_tile(data = om[platform == "NMR"], aes(fill = beta), colour = NA) +
  geom_point(data = om[platform == "NMR"], aes(size = psig),
             shape = 21, colour = "grey15", fill = NA, stroke = .75) +
  scale_fill_gradient2(low = DIV_LOW, mid = "white", high = DIV_HI, midpoint = 0,
                       limits = c(-.8, .8), breaks = c(-.8, -.4, 0, .4, .8),
                       oob = squish, name = "NMR\nmean adj.\neffect (SD)") +
  new_scale_fill() +
  geom_tile(data = om[platform == "Olink"], aes(fill = beta), colour = NA) +
  geom_point(data = om[platform == "Olink"], aes(size = psig),
             shape = 21, colour = "grey15", fill = NA, stroke = .75) +
  scale_fill_gradient2(low = DIV_LOW, mid = "white", high = DIV_HI, midpoint = 0,
                       # Group means top out at +/-0.213, so the +/-0.3 limit leaves the
                       # data occupying ~71% of the ramp. That headroom is deliberate:
                       # at a tighter limit every cell saturates and the whole block
                       # reads as uniformly dark, hiding the A-vs-C contrast.
                       limits = c(-.3, .3), breaks = c(-.3, -.15, 0, .15, .3),
                       oob = squish, name = "Olink\nmean adj.\neffect (NPX)") +
  scale_size_area("% FDR-sig", max_size = 3, breaks = c(25,50,75,100)) +
  scale_x_discrete(labels = c(A="A",B="B",C="C",D="D"), position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  facet_grid(platform ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(title = "B") +  # eBM: descriptive title moved to figure legend; panel letter only
  base_thm + theme(axis.title = element_blank(),
                   panel.grid = element_blank(),
                   panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
                   panel.spacing = unit(0, "pt"),
                   axis.text.x.top = element_text(size = 16, face = "bold", colour = "black"),
                   axis.text.y = element_text(size = 13, colour = "black"),
                   strip.text.y.left = element_text(angle = 90, face = "bold", size = 14, colour = LAYER["L2a"]),
                   strip.placement = "outside", legend.position = "right")

# ---- lock heatmap cell size (shared with Panel C: 11 mm wide x 7 mm tall) -----
# Fixes the panel BODY to an absolute size so B and C render identical cells;
# B has two facet rows (NMR 11 groups, Olink 8) -> set each to nrows x cell_h.
CELL_W <- 10; CELL_H <- 10                     # mm; MUST match Panel C
fix_cells <- function(p, panel_rows, ncol, cw, ch) {
  g <- ggplotGrob(p)
  pan <- g$layout[grepl("^panel", g$layout$name), , drop = FALSE]
  for (l in unique(pan$l)) g$widths[l]  <- grid::unit(cw * ncol, "mm")
  ts <- sort(unique(pan$t)); stopifnot(length(ts) == length(panel_rows))
  for (i in seq_along(ts)) g$heights[ts[i]] <- grid::unit(panel_rows[i] * ch, "mm")
  g
}
gB <- fix_cells(pB, panel_rows = c(11, 8), ncol = 4, cw = CELL_W, ch = CELL_H)  # NMR top, Olink bottom
ggsave(file.path(out, "Fig5_panelB_omics_heatmap.png"), gB,
       width = 185, height = 250, units = "mm", dpi = 300, bg = "white")   # generous canvas; master trims
ggsave(file.path(out, "Fig5_panelB_omics_heatmap.pdf"), gB,
       width = 185, height = 250, units = "mm", device = cairo_pdf)
cat("DONE panelB (cells", CELL_W, "x", CELL_H, "mm) ->", out, "\n")
