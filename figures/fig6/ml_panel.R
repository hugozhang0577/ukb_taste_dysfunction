# =============================================================================
# Fig 6  --  Machine-learning prediction & interpretability  (6 panels, 3x2)
#   A  ROC, 6 tiers, with smell (G1 OOF)        [square]
#   B  cross-cohort AUC forest, with smell (G1/G2/G3)
#   C  ROC, 6 tiers, smell ablated (G1 OOF)      [square]
#   D  cross-cohort AUC forest, smell ablated (G1/G2/G3)
#   E  SHAP importance bars (mean|SHAP|), with smell, coloured by direction; broken x-axis for smell
#   F  SHAP importance bars (mean|SHAP|), smell ablated
# Models labelled M1-M6 (defined in a companion table). ROC from saved OOF preds;
# forests + SHAP importance from eval CSVs (shap_*_mean_abs.csv).
# =============================================================================
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)

ev  <- "output/model_reports/eval"
out <- "output/figures/fig6"

MID  <- paste0("M", 1:6)
MIDS <- c("M1_TierA","M2_TierAB","M3_TierABC","M4_TierD_Olink","M5_TierD_NMR","M6_TierD_Full")
TAB10 <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b")   # model curves
COH     <- c(G1="#E64B35", G2="#4DBBD5", G3="#3C5488")                     # cohort (npg palette, matches fig2_nmr_forest)
COH_LAB <- c(G1="G1 (discovery)", G2="G2 (other-White)", G3="G3 (non-White)")
# ---- unified typography (eBM: no in-panel titles; bold panel letters only) --
FAM     <- "sans"     # single family across the whole figure
CEX_LET <- 1.55       # bold panel letter (A-F)
CEX_AX  <- 1.22       # axis tick labels
CEX_LAB <- 1.30       # axis titles + row labels
CEX_LEG <- 1.12       # legends / keys / value labels
# panel letter only; the descriptive title now lives in the figure legend (extra args ignored)
hdr <- function(L, ...) mtext(bquote(bold(.(L))), side = 3, line = 0.4, adj = 0,
                              cex = CEX_LET, family = FAM, xpd = NA)

PRETTY <- c(
  smell_any="Smell dysfunction (any)", smell_2w_strict_OR="Smell dysfunction, >=2 wk / impactful",
  smell_time="Smell dysfunction duration", body_fat_pct="Body fat %", neuroticism="Neuroticism",
  phq15_somatic_total="PHQ-15 somatic burden", phq9_total="PHQ-9 depression",
  health_rating="Self-rated health", pain_health_today="Pain limiting health",
  pain_chronic="Chronic pain", n_medications="No. of medications", cog_tmt_a="Trail-making A",
  pain_site_count="Pain-site count", pain_month_count="Pain sites (past month)",
  unenthusiasm="Anhedonia", hearing_diff="Hearing difficulty", abdom_pain_6mo="Abdominal pain (6 mo)",
  cog_tmt_b="Trail-making B", smoking="Non-current smoking", gad7_total="GAD-7 anxiety")
pretty <- function(x) ifelse(x %in% names(PRETTY), PRETTY[x], gsub("_", " ", x))

# ---- ROC from out-of-fold predictions --------------------------------------
mdir  <- "output/models/xgboost"
nsdir <- file.path(mdir, "sensitivity_no_smell")
roc_xy <- function(pred, y, npts = 600) {
  ok <- is.finite(pred) & !is.na(y); pred <- pred[ok]; y <- as.integer(y[ok])
  o <- order(pred, decreasing = TRUE); y <- y[o]
  tpr <- c(0, cumsum(y == 1)/sum(y == 1)); fpr <- c(0, cumsum(y == 0)/sum(y == 0))
  auc <- sum(diff(fpr)*(head(tpr,-1)+tail(tpr,-1))/2)
  idx <- unique(round(seq(1, length(fpr), length.out = npts)))
  list(fpr = fpr[idx], tpr = tpr[idx], auc = auc)
}
ROC_W <- lapply(MIDS, function(id){ d<-readRDS(file.path(mdir, paste0(id,"_oof_preds.rds"))); roc_xy(d$pred,d$y) })
ROC_N <- lapply(1:6, function(k){ d<-readRDS(file.path(nsdir, sprintf("M%da_NoSmell_oof_preds.rds",k))); roc_xy(d$pred,d$y) })

# ---- forest data -----------------------------------------------------------
es <- read.csv(file.path(ev, "xgb_eval_summary.csv"), stringsAsFactors=FALSE); es <- es[match(MIDS, es$model_id),]
xv <- read.csv(file.path(ev, "heldout_val_summary.csv"), stringsAsFactors=FALSE); xv <- xv[xv$variant=="platt",]
ns <- read.csv(file.path(ev, "sensitivity_no_smell_summary.csv"), stringsAsFactors=FALSE)
get_xv <- function(id, ch) { r <- xv[xv$model_id==id & xv$cohort==ch,]; if(nrow(r)==0) c(NA,NA,NA) else c(r$auc[1], r$auc_lo[1], r$auc_hi[1]) }
get_ns <- function(id, ch) { r <- ns[ns$original_id==id & ns$cohort==ch,]; if(nrow(r)==0) c(NA,NA,NA) else c(r$auc[1], r$auc_lo[1], r$auc_hi[1]) }
build_forest <- function(smell) {
  rows <- list()
  for (i in 1:6) { id <- MIDS[i]
    g1 <- if (smell) c(es$auc[i], es$auc_lo[i], es$auc_hi[i]) else get_ns(id,"G1")
    g2 <- if (smell) get_xv(id,"g2") else get_ns(id,"G2")
    g3 <- if (smell) get_xv(id,"g3") else get_ns(id,"G3")
    for (j in 1:3) { v <- list(g1,g2,g3)[[j]]
      rows[[length(rows)+1]] <- data.frame(mi=i, cohort=names(COH)[j], auc=v[1], lo=v[2], hi=v[3]) } }
  do.call(rbind, rows)
}
FOR_W <- build_forest(TRUE); FOR_N <- build_forest(FALSE)

# ============================== ROC panel (square) ==========================
panel_roc <- function(ROC, aucs, ci_lo, ci_hi, letter, title) {
  par(mar=c(4.4,4.6,3.0,1.0), pty="s")
  plot(NA, xlim=c(0,1), ylim=c(0,1), axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i", asp=1)
  rect(0,0,1,1, border="grey80"); abline(0,1, lty=3, col="grey65")
  axis(1, at=seq(0,1,0.2), cex.axis=CEX_AX, mgp=c(2.4,0.7,0)); axis(2, at=seq(0,1,0.2), las=1, cex.axis=CEX_AX, mgp=c(2.4,0.7,0))
  mtext("False-positive rate", side=1, line=2.6, cex=CEX_LAB); mtext("True-positive rate", side=2, line=2.8, cex=CEX_LAB)
  for (i in 6:1) lines(ROC[[i]]$fpr, ROC[[i]]$tpr, col=TAB10[i], lwd=3.0)
  leg <- sprintf("%s  %.3f (%.3f-%.3f)", MID, aucs, ci_lo, ci_hi)
  legend("bottomright", bty="n", cex=CEX_LEG, lwd=3.0, col=TAB10, legend=leg, seg.len=1.4,
         title=expression(bold("      AUROC (95% CI)")), title.adj=0)
  hdr(letter)
}

# ============================== forest panel ================================
panel_forest <- function(FD, letter, title, xlim) {
  par(mar=c(4.4,3.6,3.0,0.9), pty="m")
  plot(NA, xlim=xlim, ylim=c(0.5,6.5), axes=FALSE, xlab="", ylab="")
  axis(1, at=base::pretty(xlim,3), cex.axis=CEX_AX, mgp=c(2.3,0.6,0)); mtext("AUC (95% CI)", side=1, line=2.6, cex=CEX_LAB)
  abline(h=1:6, col="grey92"); off <- c(G1=0.24, G2=0, G3=-0.24)
  for (k in 1:nrow(FD)) { r <- FD[k,]; y <- (7-r$mi) + off[r$cohort]; cc <- COH[r$cohort]
    lo <- max(r$lo, xlim[1]); hi <- min(r$hi, xlim[2])
    segments(lo, y, hi, y, col=cc, lwd=2.8, lend=1)
    points(min(max(r$auc,xlim[1]),xlim[2]), y, pch=19, col=cc, cex=1.5) }
  for (i in 1:6) text(xlim[1]-diff(xlim)*0.04, 7-i, MID[i], adj=c(1,0.5), cex=CEX_LAB, xpd=NA, font=2)
  legend("topleft", bty="n", pch=19, col=COH, legend=names(COH), cex=CEX_LEG, horiz=TRUE, pt.cex=1.3)
  hdr(letter)
}

# ============================== feature importance ==========================
bee_offset <- function(x, halfwidth=0.36, cv=NULL) {
  n <- length(x); if (n==0) return(numeric(0))
  nbin <- max(25, round(sqrt(n)*1.5)); br <- seq(min(x), max(x), length.out=nbin+1)
  b <- findInterval(x, br, rightmost.closed=TRUE); maxc <- max(tabulate(b)); step <- 2*halfwidth/max(maxc,1)
  off <- numeric(n); for (bi in unique(b)) { ix <- which(b==bi); c <- length(ix)
    o <- if (!is.null(cv)) ix[order(cv[ix])] else ix; off[o] <- (seq_len(c)-(c+1)/2)*step }
  off
}

# ============================== SHAP importance bars ========================
# Horizontal mean|SHAP| bars, top-N features, coloured by direction. For the
# with-smell panel one feature (smell) is ~13x the rest, so a one-sided axis
# break keeps the smaller bars readable; the ablated panel needs no break.
panel_bar <- function(csv, letter, title, topn=12, brk=NA, vmax=NA, ticks=NULL) {
  d <- read.csv(csv, stringsAsFactors=FALSE)
  d <- d[order(-d$mean_abs_shap),, drop=FALSE][seq_len(min(topn, nrow(d))),, drop=FALSE]
  nf <- nrow(d); val <- d$mean_abs_shap
  bar_col <- ifelse(grepl("higher risk", d$direction), "#C0392B", "#3B6FB0")   # red up / blue down
  vmx <- if (is.na(vmax)) max(val) else vmax
  brk_on <- !is.na(brk) && vmx > brk
  if (brk_on) {                                  # piecewise x: linear 0..brk, compressed brk..vmx
    wLin <- 0.74; wCmp <- 0.26
    fx <- function(x) ifelse(x <= brk, x/brk*wLin, wLin + (x-brk)/(vmx-brk)*wCmp)
    xtop <- wLin + wCmp
    if (is.null(ticks)) ticks <- c(base::pretty(c(0,brk),3), round(vmx,1))
  } else {
    fx <- function(x) x/vmx; xtop <- 1
    if (is.null(ticks)) ticks <- base::pretty(c(0,vmx),4)
  }
  ticks <- unique(ticks[ticks >= 0 & ticks <= vmx])
  par(mar=c(4.4,18.5,3.0,3.4), pty="m")
  plot(NA, xlim=c(0, xtop*1.02), ylim=c(0.4, nf+0.6), axes=FALSE, xlab="", ylab="")
  bh <- 0.6
  for (i in seq_len(nf)) { y0 <- nf-i+1
    rect(0, y0-bh/2, fx(val[i]), y0+bh/2, col=adjustcolor(bar_col[i], 0.85), border=NA)
    text(-xtop*0.02, y0, pretty(d$feature[i]), adj=c(1,0.5), xpd=NA, cex=CEX_LAB, col="grey10")
    text(fx(val[i])+xtop*0.012, y0, sprintf("%.2f", val[i]), adj=c(0,0.5), cex=CEX_LEG, col="grey30", xpd=NA) }
  axis(1, at=fx(ticks), labels=ticks, cex.axis=CEX_AX, mgp=c(2.3,0.6,0))
  if (brk_on) { bx <- wLin; usr <- par("usr"); dy <- (usr[4]-usr[3])*0.015; dx <- xtop*0.006
    abline(v=bx, col="grey88", lwd=0.8)
    for (s in c(0, 2.6*dx)) segments(bx-dx+s, usr[3]-dy, bx+dx+s, usr[3]+dy, col="grey30", lwd=1.3, xpd=NA) }
  mtext("mean |SHAP| (log-odds)", side=1, line=2.6, cex=CEX_LAB)
  legend("bottomright", bty="n", cex=CEX_LEG*0.85, pch=15, col=c("#C0392B","#3B6FB0"), inset=c(0, 0.015),
         legend=c("increases predicted risk","decreases predicted risk"))
  hdr(letter)
}

# ============================== assemble ====================================
draw_fig <- function() {
  # top: A(square ROC) B(slim forest) C(square ROC) D(slim forest); bottom: E F importance
  # top row shortened (heights 0.84) so the square ROC panels fill their cells and
  # align with the forests rather than floating with centred whitespace.
  layout(matrix(c(1,2,3,4, 5,5,6,6), 2, 4, byrow=TRUE), widths=c(1,0.64,1,0.64), heights=c(0.84,1))
  panel_roc(ROC_W, es$auc, es$auc_lo, es$auc_hi, "A", "ROC - all features")
  panel_forest(FOR_W, "B", "Cross-cohort AUC", c(0.70,1.0))
  nsg1 <- sapply(MIDS, function(id) get_ns(id,"G1"))
  panel_roc(ROC_N, nsg1[1,], nsg1[2,], nsg1[3,], "C", "ROC - smell ablated")
  panel_forest(FOR_N, "D", "Cross-cohort (ablated)", c(0.50,0.78))
  panel_bar(file.path(ev,"shap_M1_TierA_mean_abs.csv"), "E", "SHAP importance - all features",
            brk=0.25, vmax=0.95, ticks=c(0,0.1,0.2,0.5,0.9))
  panel_bar(file.path(ev,"shap_M1a_NoSmell_mean_abs.csv"), "F", "SHAP importance - smell ablated")
}
# No in-figure title (eBM: title lives in the legend only; panel letters A-F carry the figure).
draw_all <- function() {
  par(oma=c(0.5,0.5,0.5,0.5), family=FAM); draw_fig()
}
png(file.path(out,"Fig6_ml_panel.png"), width=2300, height=1650, res=200); draw_all(); dev.off()   # 11.5 x 8.25 in
pdf(file.path(out,"Fig6_ml_panel.pdf"), width=11.5, height=8.25); draw_all(); dev.off()
cat("DONE ->", out, "\n")
