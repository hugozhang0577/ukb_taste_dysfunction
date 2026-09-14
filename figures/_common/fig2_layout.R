FIG2 <- data.frame(
  letter = c("a", "b", "c", "d"),
  stem   = c("fig2a_gwas_manhattan", "fig2b_apoe_forest",
             "fig2c_protein_forest", "fig2d_metabolite_forest"),
  x      = c(0.0,  0.0,  0.0,  86.5),
  y      = c(0.0, 44.0, 93.6,  44.0),
  w      = c(168.0, 81.5, 81.5, 81.5),
  h      = c(40.0,  45.6, 82.4, 132.0),
  stringsAsFactors = FALSE)

FIG2_W <- 168.0
FIG2_H <- max(FIG2$y + FIG2$h)
FIG2_LETTER_PT <- 9

stopifnot(FIG2_H <= 207)
stopifnot(abs((FIG2$y[3] + FIG2$h[3]) - (FIG2$y[4] + FIG2$h[4])) < 0.05)
stopifnot(abs((FIG2$x[4] - (FIG2$x[2] + FIG2$w[2])) - 5.0) < 0.05)

cat(sprintf("[fig2_layout] %.1f x %.1f mm | columns flush | gutter %.1f mm\n",
            FIG2_W, FIG2_H, FIG2$x[4] - (FIG2$x[2] + FIG2$w[2])))
