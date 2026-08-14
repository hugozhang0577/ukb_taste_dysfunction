# Fig 4 — raster composite of the four panels.
#
# Layout:  top    = A heatmap | (B male UMAP over C female UMAP)
#          bottom = D subtype-signature forest, spanning the top row's width
#
# This produces a single reviewable page from the four panel scripts. The
# submitted figure is assembled from the panels' vector (PDF) outputs, so the
# raster produced here is for checking layout and proportions, not for print.
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
suppressPackageStartupMessages(library(magick))
DIR <- "output/figures/fig4"
rd  <- function(f) image_trim(image_read(file.path(DIR, f)))
lab <- function(img, L) image_annotate(img, L, size = 60, weight = 700,
                                        color = "black", boxcolor = "white", location = "+10+6")

A <- lab(rd("fig4_feature_sample_heatmap.png"), "A")
B <- lab(rd("fig4_panelB_umap_male.png"),       "B")
C <- lab(rd("fig4_panelC_umap_female.png"),     "C")
D <- rd("fig4_panelD_signature_forest.png")          # subtype panels self-labelled A-D inside

H   <- 1650                                       # top-row height
Ac  <- image_scale(A, paste0("x", H))
bc  <- image_scale(image_append(c(B, C), stack = TRUE), paste0("x", H))
top <- image_append(c(Ac, bc))
Wt  <- image_info(top)$width

# panel letter D in a left margin band so it does not collide with the subtype
# 'A Aging frailty' title inside the strip
Dw   <- image_scale(D, as.character(Wt - 90))
Dw   <- image_border(Dw, "white", "45x0")
Dw   <- image_annotate(Dw, "D", size = 60, weight = 700, color = "black",
                       gravity = "NorthWest", location = "+4+2")
Dw   <- image_extent(Dw, geometry_size_pixels(Wt, image_info(Dw)$height), color = "white")

page <- image_append(c(top, Dw), stack = TRUE)
page <- image_border(page, "white", "24x24")
out  <- file.path(DIR, "Fig4_assembled.png")
image_write(page, out)
cat("wrote", out, " ", image_info(page)$width, "x", image_info(page)$height, "\n")
