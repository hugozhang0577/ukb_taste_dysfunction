#> out  output/figures/_fig3_layout.csv
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)

FIG3 <- data.frame(
  letter = c("a", "b", "c", "d", "e"),
  stem   = c("fig3a_feature_sample_heatmap", "fig3b_subtype_signature_forest",
             "fig3c_olink_heatmap", "fig3d_nmr_heatmap", "fig3e_disease_heatmap"),
  x      = c(0.0,   0.0,   0.0,  57.7, 115.4),
  y      = c(0.0,  93.5, 144.8, 144.8, 144.8),
  w      = c(168.0, 168.0, 52.7, 52.7, 52.7),
  h      = c(90.0,   47.8, 48.4, 58.3, 61.6),
  stringsAsFactors = FALSE)

FIG3_W <- 168.0
FIG3_H <- max(FIG3$y + FIG3$h)
FIG3_LETTER_PT <- 9

stopifnot(FIG3_H <= 207)
stopifnot(length(unique(FIG3$y[3:5])) == 1)
stopifnot(all(abs(diff(FIG3$x[3:5]) - (FIG3$w[3] + 5.0)) < 0.05))

LAYOUT_CSV <- "output/figures/_fig3_layout.csv"
dir.create(dirname(LAYOUT_CSV), showWarnings = FALSE, recursive = TRUE)
utils::write.csv(cbind(FIG3, canvas_w = FIG3_W, canvas_h = FIG3_H,
                       letter_pt = FIG3_LETTER_PT), LAYOUT_CSV, row.names = FALSE)
cat(sprintf("[fig3_layout] %.1f x %.1f mm | gutter %.1f mm -> %s\n",
            FIG3_W, FIG3_H, FIG3$x[4] - (FIG3$x[3] + FIG3$w[3]), LAYOUT_CSV))
