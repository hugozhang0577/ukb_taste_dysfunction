# =============================================================================
# Fig 5 master assembler  --  3-row layout with a shared subtype legend
#   row 1  A               sensitivity (full width)
#   row 2  B | C            omics | disease heatmaps
#   row 3  D               merged genetics + severity (3 sub-panels horizontal)
#   foot   shared subtype colour legend (A/B/C/D), centred
# Panels are mixed base-graphics (A) + ggplot (B,C,D) -> composed as PNG tiles with
# magick. The subtype colour key is generated once here and shared by all panels
# (A keeps only "reassigned"; D's severity bar has no per-panel subtype legend).
# =============================================================================
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
suppressPackageStartupMessages({ library(magick); library(ggplot2) })

REBUILD <- FALSE   # TRUE = run each panel script first (panel A is slow)
# Panel PNGs live under output/; the panel scripts themselves sit next to this file.
out    <- "output/figures/fig5/panels"
final  <- "output/figures/fig5"
pandir <- Sys.getenv("CODE_DIR", unset = ".")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
SUB_PAL <- c(A = "#4E79A7", B = "#E15759", C = "#F28E2B", D = "#59A14F")
# Legend maps the A-D letters (used on the heatmaps + forest y-axes) to subtype names.
NAMES   <- c(A = "A  Ageing frailty", B = "B  Psychosomatic",
             C = "C  Cardiometabolic", D = "D  Young idiopathic")

if (REBUILD) for (s in c("panelA_reproducibility_radial", "panelB_omics_heatmap",
                         "panelC_disease_heatmap", "panelDE_genetics_covid"))
  { cat("rebuild", s, "\n"); source(file.path(pandir, paste0(s, ".R")), local = new.env()) }

# ---- shared subtype colour legend (generate, then trim to content) -----------
legdf <- data.frame(s = factor(c("A","B","C","D"), levels = c("A","B","C","D")))
legp <- ggplot(legdf, aes(s, 1, fill = s)) + geom_tile(alpha = 0) +
  scale_fill_manual(values = SUB_PAL, labels = NAMES, name = NULL,
                    guide = guide_legend(nrow = 1, override.aes = list(alpha = 1))) +
  theme_void() +
  theme(legend.position = "bottom", legend.text = element_text(size = 13),
        legend.key.size = unit(7, "mm"), legend.spacing.x = unit(4, "mm"))
ggsave(file.path(out, "_legend_strip.png"), legp, width = 200, height = 18,
       units = "mm", dpi = 300, bg = "white")
leg <- image_trim(image_read(file.path(out, "_legend_strip.png")))

# ---- read panels -------------------------------------------------------------
png_path <- function(stem) {
  f <- file.path(out, paste0(stem, ".png"))
  if (!file.exists(f)) stop("missing panel PNG (run panels/*.R first): ", f)
  image_read(f)
}
A  <- png_path("Fig5_panelA_radial_8panel")
B  <- png_path("Fig5_panelB_omics_heatmap")
C  <- png_path("Fig5_panelC_disease_heatmap")
Dh <- png_path("Fig5_panelD_hla")          # HLA forest alone -> tucked under C
Ecs<- png_path("Fig5_panelE_covidsev")     # COVID forest | severity -> bottom row

W   <- 2400L; GAP <- 26L; PAD <- 22L
to_w <- function(img, w) image_scale(img, as.character(w))
row_of <- function(imgs, ws) {
  parts <- Map(to_w, imgs, ws)
  h <- max(sapply(parts, function(im) image_info(im)$height))
  parts <- lapply(parts, function(im) image_extent(im, sprintf("%dx%d", image_info(im)$width, h),
                                                   gravity = "north", color = "white"))
  image_append(image_join(parts))
}
# vertical stack of two images at a fixed width (top over bottom, small gap)
col_stack <- function(top, bot, w) {
  parts <- list(to_w(top, w), image_blank(w, GAP, "white"), to_w(bot, w))
  image_append(image_join(parts), stack = TRUE)
}

# ---- rows --------------------------------------------------------------------
# B and C are rendered with IDENTICAL heatmap cells (11x7 mm, locked in the panel
# scripts) on generous white canvases; trim them so the cells align, keep them at
# NATIVE resolution through the B|[C/D] row, then scale the whole row once so the
# equal-cell property survives (no per-panel width forcing).
Btr <- image_trim(B); Ctr <- image_trim(C); Dtr <- image_trim(Dh)
Cw <- image_info(Ctr)$width; Bh <- image_info(Btr)$height; Ch <- image_info(Ctr)$height
# D (HLA forest) auto-fills the gap under C so the C/D column height matches B.
Dnat <- image_info(image_scale(Dtr, as.character(Cw)))$height    # D height at C's width
gapH <- max(Bh - Ch - GAP, Dnat)                                 # never squash below natural
Dfit <- image_extent(image_scale(Dtr, sprintf("%dx%d", Cw, gapH)),
                     sprintf("%dx%d", Cw, gapH), gravity = "center", color = "white")
rightcol <- image_append(image_join(list(Ctr, image_blank(Cw, GAP, "white"), Dfit)), stack = TRUE)
h2   <- max(Bh, image_info(rightcol)$height)
padN <- function(im) image_extent(im, sprintf("%dx%d", image_info(im)$width, h2),
                                  gravity = "north", color = "white")
row2 <- image_append(image_join(list(padN(Btr), image_blank(GAP, h2, "white"), padN(rightcol))))
rowA   <- to_w(A, W)
rowBCD <- to_w(row2, W)                             # B | [C over D], cells equal
rowE   <- to_w(Ecs, W)                              # E = COVID + severity (full width)
leg_s <- image_scale(leg, as.character(round(W * 0.82)))            # scale to fit width (no clip)
leg_c <- image_extent(leg_s, sprintf("%dx%d", W, image_info(leg_s)$height + 24),
                      gravity = "center", color = "white")          # centre under width W

gap <- image_blank(W, GAP, "white")
stack <- image_append(image_join(list(rowA, gap, rowBCD, gap, rowE, gap, leg_c)), stack = TRUE)
fig <- image_border(stack, "white", sprintf("%dx%d", PAD, PAD))

info <- image_info(fig); cat(sprintf("composite: %d x %d px\n", info$width, info$height))
image_write(fig, file.path(final, "Fig5_master_preview.png"), format = "png")
image_write(fig, file.path(final, "Fig5_master_preview.pdf"), format = "pdf")
cat("DONE ->", file.path(final, "Fig5_master_preview.{png,pdf}"), "\n")
