#> out  output/figures/Figure2.pdf
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
suppressPackageStartupMessages({
  library(ggplot2); library(grid); library(patchwork)
})

grDevices::pdfFonts(Arial        = grDevices::pdfFonts()$Helvetica)
grDevices::postscriptFonts(Arial = grDevices::postscriptFonts()$Helvetica)

FIGDIR <- "figures"
OUT    <- "output/figures"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(FIGDIR, "_common/fig2_layout.R"))

CACHE_A <- file.path(OUT, "_cache_fig2a_grob.rds")

## ---- helpers ---------------------------------------------------------------
panel_from <- function(script, obj = "p") {
  e <- new.env(parent = globalenv())
  cat(sprintf("[build] %s\n", script))
  sys.source(file.path(FIGDIR, script), envir = e, toplevel.env = e)
  if (!exists(obj, envir = e, inherits = FALSE))
    stop(script, " did not leave an object called '", obj, "'")
  get(obj, envir = e)
}

as_grob <- function(p) {
  # a patchwork also inherits "ggplot", so test for patchwork FIRST
  if (inherits(p, "patchwork")) patchwork::patchworkGrob(p) else ggplot2::ggplotGrob(p)
}

## ---- panel a: cached, because the Manhattan reads 10.4M variants ------------
if (file.exists(CACHE_A)) {
  cat(sprintf("[cache] reusing %s (delete it to rebuild the Manhattan)\n", CACHE_A))
  gA <- readRDS(CACHE_A)
} else {
  cat("[build] fig2/a_gwas_manhattan.R  -- slow, ~3 min, caching the result\n")
  e <- new.env(parent = globalenv())
  sys.source(file.path(FIGDIR, "fig2/a_gwas_manhattan.R"), envir = e, toplevel.env = e)
  if (is.null(e$pA)) stop("panel a is NULL: the script ran with GWAS disabled")
  gA <- as_grob(e$strip_tag(e$pA))
  saveRDS(gA, CACHE_A)
}

## ---- panels b, c, d: fast, always rebuilt -----------------------------------
gB <- as_grob(panel_from("fig2/b_apoe_forest.R"))
gC <- as_grob(panel_from("fig2/c_protein_forest.R"))
gD <- as_grob(panel_from("fig2/d_metabolite_forest.R"))
grobs <- list(a = gA, b = gB, c = gC, d = gD)

## ---- draw ------------------------------------------------------------------
f_pdf <- file.path(OUT, "Figure2.pdf")
cairo_pdf(f_pdf, width = FIG2_W / 25.4, height = FIG2_H / 25.4, bg = "white")
grid.newpage()
for (i in seq_len(nrow(FIG2))) {
  # y in the layout table runs DOWNWARD from the top edge; grid runs upward
  pushViewport(viewport(
    x = unit(FIG2$x[i], "mm"), y = unit(FIG2_H - FIG2$y[i], "mm"),
    width = unit(FIG2$w[i], "mm"), height = unit(FIG2$h[i], "mm"),
    just = c("left", "top")))
  grid.draw(grobs[[FIG2$letter[i]]])
  popViewport()
}
for (i in seq_len(nrow(FIG2))) {
  grid.text(FIG2$letter[i],
            x = unit(FIG2$x[i] + 1.2, "mm"),
            y = unit(FIG2_H - FIG2$y[i] - 1.2, "mm"),
            just = c("left", "top"),
            gp = gpar(fontfamily = "Arial", fontface = "bold",
                      fontsize = FIG2_LETTER_PT, col = "black"))
}
invisible(dev.off())

## ---- verify the file we just wrote ------------------------------------------
sz <- file.info(f_pdf)$size
cat(sprintf("\n[assembled] %s  (%.1f x %.1f mm, %.1f MB)\n",
            f_pdf, FIG2_W, FIG2_H, sz / 1024^2))
cat("[note] Text is live text and lines are paths; only the Manhattan point cloud\n")
cat("       is raster, rasterised at 826 dpi inside its own panel by design.\n")
cat("[note] Delete out/_cache_fig2a_grob.rds after any change to the Manhattan.\n")
