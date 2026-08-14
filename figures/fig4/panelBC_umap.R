# ===========================================================================
# Fig 4 — Panels B and C: MOFA factor UMAP, one square panel per sex
#
# The model is fitted separately in men and women, so the embedding is drawn
# per sex rather than as one faceted panel: a shared canvas would imply the two
# sexes share a factor space, which they do not.
#
# Subtype colours and labels follow the manuscript's user-facing names and match
# Fig 5, so a subtype keeps one colour across the whole paper.
#
# Outputs (square canvases; multi-panel assembly is done by assemble.R):
#   output/figures/fig4/fig4_panelB_umap_male.{pdf,png}
#   output/figures/fig4/fig4_panelC_umap_female.{pdf,png}
#
# Reads (produced by subtyping/):
#   output/subtyping/reports/fig_tables/umap_coords_{m,f}_k4.csv
# ===========================================================================

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

OUT_DIR <- "output/figures/fig4"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Subtype mapping (cluster id -> letter); single source of truth.
source("unsupervised_subtyping/_subtype_map.R")

# User-facing subtype names and palette, matching Figure 5.
SUBTYPE_LABEL <- c(A = "A. Aging frailty",
                B = "B. Psychosomatic",
                C = "C. Cardiometabolic",
                D = "D. Young idiopathic")
SUBTYPE_COLOR <- c("A. Aging frailty"    = "#4E79A7",
                "B. Psychosomatic"    = "#E15759",
                "C. Cardiometabolic"  = "#F28E2B",
                "D. Young idiopathic" = "#59A14F")


# ---------------------------------------------------------------------------
# Panels B & C — MOFA factor UMAP, one square panel per sex
# ---------------------------------------------------------------------------
umap_plot <- function(dt, sex_key, sex_label) {
  # Derive the user-facing label from the cluster id via _subtype_map rather than
  # from the file's own `subtype` column, so the figure and the map cannot drift.
  dt <- copy(dt)
  dt[, letter := SUBTYPE_MAP[[sex_key]][as.character(cluster)]]
  dt[, subtype := factor(SUBTYPE_LABEL[letter], levels = SUBTYPE_LABEL)]
  stopifnot(!any(is.na(dt$subtype)))
  ggplot(dt, aes(UMAP1, UMAP2, colour = subtype)) +
    geom_point(size = 0.7, alpha = 0.65) +
    scale_colour_manual(values = SUBTYPE_COLOR, guide = "none") +   # colour key lives in Panel A
    labs(subtitle = sex_label) +
    theme_bw(base_size = 12) +   # calibrated so on-page fonts match A / D / E after assembly
    theme(panel.grid.minor = element_blank(),
          plot.subtitle    = element_text(face = "bold", size = 14),
          axis.title       = element_text(size = 13),
          axis.text        = element_text(size = 11),
          aspect.ratio     = 1)
}

umap_m <- fread("output/subtyping/reports/fig_tables/umap_coords_m_k4.csv")
umap_f <- fread("output/subtyping/reports/fig_tables/umap_coords_f_k4.csv")

p_m <- umap_plot(umap_m, "m", sprintf("Male (n=%s)",   format(nrow(umap_m), big.mark = ",")))
p_f <- umap_plot(umap_f, "f", sprintf("Female (n=%s)", format(nrow(umap_f), big.mark = ",")))

# Uniform 5x5 square canvas for B/C/D (matched panel size in the assembled column)
cat("[panel B] male UMAP  n=", nrow(umap_m), "\n", sep = "")
ggsave(file.path(OUT_DIR, "fig4_panelB_umap_male.pdf"), p_m, width = 5, height = 5)
ggsave(file.path(OUT_DIR, "fig4_panelB_umap_male.png"), p_m, width = 5, height = 5, dpi = 300)

cat("[panel C] female UMAP  n=", nrow(umap_f), "\n", sep = "")
ggsave(file.path(OUT_DIR, "fig4_panelC_umap_female.pdf"), p_f, width = 5, height = 5)
ggsave(file.path(OUT_DIR, "fig4_panelC_umap_female.png"), p_f, width = 5, height = 5, dpi = 300)


cat("\n=== Done ===\n")
cat("Panel B: fig4_panelB_umap_male.{pdf,png}\n")
cat("Panel C: fig4_panelC_umap_female.{pdf,png}\n")
