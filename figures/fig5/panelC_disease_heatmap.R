# =============================================================================
# Fig 5 Panel C  --  disease phenome: curated PheCode enrichment x subtype
# Rendered as a standalone image and composited by fig5/assemble.R.
# Theme, colours, geoms, scales and strip styling are shared with the other panels.
# ASCII glyphs only; cairo_pdf safe.  data.table only.
# =============================================================================
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(scales); library(grid) })

ev  <- "output/subtyping"
out <- "output/figures/fig5/panels"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ---- shared aesthetics (identical across the fig5 panels) -------------------
LAYER   <- c(L1 = "#4E9080", L2a = "#2C6FAC", L2b = "#D4813A", L3 = "#9B59B6", L4 = "#C0392B")
DIV_LOW <- "#2166AC"; DIV_HI <- "#B2182B"

# Fonts scale-compensated (panel composited at ~0.64x) to match the other panels.
base_thm <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title    = element_text(face = "bold", size = 20, colour = "black"),
        axis.text     = element_text(colour = "black"),
        legend.title  = element_text(size = 13), legend.text = element_text(size = 12),
        legend.key.height = unit(4, "mm"), legend.key.width = unit(4, "mm"))

# ---- data --------------------------------------------------------------------
# Per-subtype PheCode enrichment, from
# subtyping/02_phecode_enrichment.R.
phe <- fread("output/subtyping/evidence_disease/phecode_enrichment_g1_by_subtype.csv",
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


# ---- panel c ------------------------------------------------------------------
pC <- ggplot(dc, aes(subtype, label, fill = signed)) +
  geom_tile(colour = NA) +
  geom_text(aes(label = star), size = 5.5, fontface = "bold", vjust = .78, colour = "black") +
  scale_fill_gradient2(low = DIV_LOW, mid = "white", high = DIV_HI, midpoint = 0,
                       limits = c(-12,12), oob = squish,
                       name = "signed\n-log10 q", breaks = c(-10,0,10),
                       labels = c("depleted","ns","enriched")) +
  scale_x_discrete(labels = c(A="A",B="B",C="C",D="D"), position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title = "C") +  # eBM: descriptive title moved to figure legend; panel letter only
  base_thm + theme(axis.title = element_blank(),
                   panel.grid = element_blank(),
                   panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
                   axis.text.x.top = element_text(size = 16, face = "bold", colour = "black"),
                   axis.text.y = element_text(size = 13, colour = "black"), legend.position = "right")

# ---- lock heatmap cell size (shared with Panel B: 11 mm wide x 7 mm tall) -----
CELL_W <- 10; CELL_H <- 10                     # mm; MUST match Panel B
fix_cells <- function(p, panel_rows, ncol, cw, ch) {
  g <- ggplotGrob(p)
  pan <- g$layout[grepl("^panel", g$layout$name), , drop = FALSE]
  for (l in unique(pan$l)) g$widths[l]  <- grid::unit(cw * ncol, "mm")
  ts <- sort(unique(pan$t)); stopifnot(length(ts) == length(panel_rows))
  for (i in seq_along(ts)) g$heights[ts[i]] <- grid::unit(panel_rows[i] * ch, "mm")
  g
}
gC <- fix_cells(pC, panel_rows = 12, ncol = 4, cw = CELL_W, ch = CELL_H)   # 12 diseases x 4 subtypes
ggsave(file.path(out, "Fig5_panelC_disease_heatmap.png"), gC,
       width = 170, height = 175, units = "mm", dpi = 300, bg = "white")   # generous canvas; master trims
ggsave(file.path(out, "Fig5_panelC_disease_heatmap.pdf"), gC,
       width = 170, height = 175, units = "mm", device = cairo_pdf)
cat("DONE panelC (cells", CELL_W, "x", CELL_H, "mm) ->", out, "\n")
