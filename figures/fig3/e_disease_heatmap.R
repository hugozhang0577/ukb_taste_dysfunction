#> in   output/subtyping/evidence_disease/phecode_enrichment_g1_by_subtype.csv
#> out  output/figures/fig3e_disease_heatmap.pdf
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
source("figures/_common/style.R")
source("figures/_common/fig3_type.R")   # local 10% type bump
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(scales); library(grid) })

ev  <- "output/subtyping"
out <- FIG_OUT
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ---- shared aesthetics -------------------------------------------------------
DIV_LOW <- "#2166AC"; DIV_HI <- "#B2182B"

# Fonts scale-compensated (panel composited at ~0.64x) to match the other panels.
base_thm <- theme_minimal(base_size = FIG3_SZ$axtxt, base_family = FIG_FONT) +
  theme(panel.grid.minor = element_blank(),
        plot.title    = element_blank(),          # panel letters are added at assembly
        axis.text     = element_text(colour = "black"),
        legend.title  = element_text(size = FIG3_SZ$lgd, face = "bold"),
        legend.text = element_text(size = FIG3_SZ$lgd),
        legend.key.height = unit(4, "mm"), legend.key.width = unit(4, "mm"))

# ---- data --------------------------------------------------------------------
phe <- fread(file.path(ev, "evidence_disease/phecode_enrichment_g1_by_subtype.csv"),
             colClasses = list(character = "phecode"))
cat(sprintf("phecode rows read: %d (%d distinct codes)\n", nrow(phe), uniqueN(phe$phecode)))
cur <- data.table(
  phecode = c("250.2","401.1","411.2","585","274","573",
              "535","574","531","562","550","721"),
  label   = c("Type 2 diabetes","Hypertension","Myocardial infarction","Renal failure",
              "Gout","Liver disease","Gastritis / duodenitis","Cholelithiasis",
              "Peptic ulcer","Diverticulosis","Abdominal hernia","Spondylosis"))
dc <- merge(phe[phecode %in% cur$phecode, .(phecode, subtype, OR, q)], cur, by = "phecode")
cat(sprintf("curated cells: %d (expect 48 = 12 x 4)\n", nrow(dc)))
dc[, signed := sign(log(OR)) * pmin(-log10(q), 12)]
dc[, star := fifelse(q < 0.001, "***", fifelse(q < 0.01, "**", fifelse(q < 0.05, "*", "")))]
ord <- dc[subtype == "C"][order(signed), label]
dc[, label := factor(label, ord)]

# ---- panel c (verbatim geoms/scales/theme; only leading letter -> C) --------
pC <- ggplot(dc, aes(subtype, label, fill = signed)) +
  geom_tile(colour = NA) +
  geom_text(aes(label = star), size = 1.9 * FIG3_K, fontface = "bold", vjust = .78,
            family = FIG_FONT, colour = "black") +   # stars sized for a 10 mm cell were 5.5
  scale_fill_gradient2(low = DIV_LOW, mid = "white", high = DIV_HI, midpoint = 0,
                       limits = c(-12,12), oob = squish,
                       # legend style matches panel b: title on top, bar 12 x 1.8 mm
                       name = "Signed -log10 q", breaks = c(-10,0,10),
                       labels = c("depleted","ns","enriched"),
                       guide = guide_colourbar(title.position = "top",
                                               barwidth = unit(15, "mm"),   # 12 mm crowds the three words; 20 mm pushes "enriched" out at this type size
                                               barheight = unit(1.8, "mm"))) +
  scale_x_discrete(labels = c(A="A",B="B",C="C",D="D"), position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title = NULL) +
  base_thm + theme(axis.title = element_blank(),
                   panel.grid = element_blank(),
                   panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
                   axis.text.x.top = element_text(size = FIG3_SZ$anno, face = "bold", colour = "black"),
                   axis.text.y = element_text(size = FIG3_SZ$axtxt, colour = "black"),
                   # legend below, as in fig4b: a right-hand column does not fit at 81.5 mm
                   legend.position = "bottom", legend.box = "horizontal",
                   legend.margin = margin(1, 0, 0, 0, "mm"),
                   legend.box.spacing = unit(1, "mm"))

# ---- lock heatmap cell size (shared with Panel B: 11 mm wide x 7 mm tall) -----
CELL_W <- round(5.0*FIG3_K,2); CELL_H <- round(3.0*FIG3_K,2)                   # mm; MUST match fig4cd_omics_heatmaps.R
fix_cells <- function(p, panel_rows, ncol, cw, ch) {
  g <- ggplotGrob(p)
  pan <- g$layout[grepl("^panel", g$layout$name), , drop = FALSE]
  for (l in unique(pan$l)) g$widths[l]  <- grid::unit(cw * ncol, "mm")
  ts <- sort(unique(pan$t)); stopifnot(length(ts) == length(panel_rows))
  for (i in seq_along(ts)) g$heights[ts[i]] <- grid::unit(panel_rows[i] * ch, "mm")
  g
}
gC <- fix_cells(pC, panel_rows = 12, ncol = 4, cw = CELL_W, ch = CELL_H)   # 12 diseases x 4 subtypes
fig_save(gC, "fig3e_disease_heatmap", w_mm = 52.7, h_mm = round(12 * CELL_H + 20.0 * FIG3_K, 1))
