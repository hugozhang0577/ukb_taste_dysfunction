suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(scales); library(grid)
})

# The journal font whitelist is Arial / Times / Frutiger / Sabon. A bare Linux
# worker has none of them; install msttcorefonts, or set FIG_FONT to one that is
# present and expect to restyle before submission.
FIG_FONT <- Sys.getenv("FIG_FONT", unset = "Arial")
if (requireNamespace("systemfonts", quietly = TRUE)) {
  .m <- try(systemfonts::match_fonts(FIG_FONT), silent = TRUE)
  if (inherits(.m, "try-error") || !nzchar(.m$path[1]))
    stop(FIG_FONT, " not found. Install it (Debian: ttf-mscorefonts-installer) ",
         "or set FIG_FONT to an installed whitelist font.")
}

# Page geometry (mm). W_TARGET leaves 5 mm of slack for vector assembly.
W_FULL <- 173.0; W_COL <- 84.5; H_MAX <- 207.0; W_TARGET <- 168.0

# Panels are placed 1:1, so these point sizes are the printed sizes.
SZ <- list(base = 7.0, axtit = 7.0, axtxt = 6.2, strip = 6.8,
           lgd = 6.2, anno = 6.0, tag = 8.5)
LW <- list(axis = 0.25, ref = 0.30, data = 0.55, err = 0.45, tick = 0.25)

theme_fig <- function(base_size = SZ$base, grid = c("y", "x", "both", "none")) {
  grid <- match.arg(grid)
  th <- theme_bw(base_size = base_size, base_family = FIG_FONT) +
    theme(
      text             = element_text(family = FIG_FONT, colour = "black"),
      plot.title       = element_blank(),
      plot.subtitle    = element_blank(),
      plot.background  = element_blank(),
      panel.background = element_blank(),
      panel.border     = element_blank(),
      axis.line        = element_line(colour = "black", linewidth = LW$axis),
      axis.ticks       = element_line(colour = "black", linewidth = LW$tick),
      axis.ticks.length = unit(1.1, "pt"),
      axis.title       = element_text(size = SZ$axtit, colour = "black"),
      axis.text        = element_text(size = SZ$axtxt, colour = "black"),
      strip.background = element_blank(),
      strip.text       = element_text(size = SZ$strip, face = "bold", colour = "black"),
      legend.background = element_blank(),
      legend.key       = element_blank(),
      legend.title     = element_text(size = SZ$lgd, face = "bold"),
      legend.text      = element_text(size = SZ$lgd),
      legend.key.size  = unit(3.2, "mm"),
      legend.margin    = margin(0, 0, 0, 0),
      plot.margin      = margin(1.5, 1.5, 1.5, 1.5, "mm")
    )
  gl <- element_line(colour = "grey92", linewidth = 0.18)
  th + switch(grid,
    y    = theme(panel.grid.major.y = gl, panel.grid.major.x = element_blank(),
                 panel.grid.minor = element_blank()),
    x    = theme(panel.grid.major.x = gl, panel.grid.major.y = element_blank(),
                 panel.grid.minor = element_blank()),
    both = theme(panel.grid.major = gl, panel.grid.minor = element_blank()),
    none = theme(panel.grid = element_blank()))
}

# A subtype keeps one colour across every figure.
SUBTYPE_COL <- c(A = "#4E79A7", B = "#E15759", C = "#F28E2B", D = "#59A14F")
COHORT_COL  <- c(G1 = "#D1495B", G2 = "#4CC3D9", G3 = "#2E4057")

# Exposure and diagnosis palettes are capped at 7 series for accessibility.
EXPO6_COL <- c(
  "Oral health"                   = "#4E79A7",
  "Symptoms & mental health"      = "#E15759",
  "Lifestyle & diet"              = "#59A14F",
  "Anthropometric & biochemistry" = "#F28E2B",
  "Sociodemographic"              = "#B07AA1",
  "Other"                         = "#9C9C9C")
DIS7_COL <- c(
  "Nervous & sensory"      = "#59A14F",
  "Mental & behavioural"   = "#B07AA1",
  "Digestive (incl. oral)" = "#4E79A7",
  "Respiratory & ENT"      = "#F28E2B",
  "Infectious"             = "#E15759",
  "Circulatory"            = "#8C564B",
  "Other"                  = "#9C9C9C")

GREY_NS  <- "#BDBDBD"
REF_LINE <- "#8B0000"

FIG_OUT <- "output/figures"

# One vector PDF per panel, text left editable; the PNG is a visual check only.
fig_save <- function(p, name, w_mm, h_mm, preview = TRUE, dpi_preview = 200) {
  stopifnot(w_mm <= W_FULL + 1e-6, h_mm <= H_MAX + 1e-6)
  dir.create(FIG_OUT, showWarnings = FALSE, recursive = TRUE)
  f_pdf <- file.path(FIG_OUT, paste0(name, ".pdf"))
  ggsave(f_pdf, p, width = w_mm, height = h_mm, units = "mm",
         device = grDevices::cairo_pdf, bg = "white")
  if (preview)
    ggsave(file.path(FIG_OUT, paste0(name, ".png")), p,
           width = w_mm, height = h_mm, units = "mm", dpi = dpi_preview, bg = "white")
  cat(sprintf("[save] %-34s %6.1f x %6.1f mm\n", name, w_mm, h_mm))
  invisible(f_pdf)
}

fmt_or <- function(x) sprintf("%.2f", x)
fmt_ci <- function(lo, hi) sprintf("%.2f to %.2f", lo, hi)
fmt_p  <- function(p) ifelse(p < 1e-4, "<0.0001", formatC(p, format = "g", digits = 2))
