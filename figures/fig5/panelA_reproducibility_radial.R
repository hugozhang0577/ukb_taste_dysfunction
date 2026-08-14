# =============================================================================
# Fig 5 Panel A  --  subtype reproducibility, radial 8-panel minimal-label layout
#
# 2 rows (F / M) x 4 cols (Seed re-runs | Balanced weighting | Broadened 6-view | Broadened+balanced).
# On-figure labels (minimal addressing only):
#   - column word headers (top row only)        - row = sex (outer left)
#   - one kept% number per cell                 - A/B/C/D on top-left tree only
#   - one shared colour legend
# Backbone = within-subtype monophyletic clades on the main fit (pure by
# construction). Terminal pie = kept main subtype (colour) vs reassigned (grey).
# Subtype letters A-D are never reused as panel ids.
# =============================================================================
suppressWarnings(suppressMessages(library(ape)))

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)

base <- "output/subtyping"
cl   <- file.path(base, "clusters")
outdir <- "output/figures/fig5/panels"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

SUBTYPE_MAP <- list(m = c("1"="D","2"="A","3"="B","4"="C"),
                 f = c("1"="C","2"="B","3"="A","4"="D"))
PAL <- c(A="#4E79A7", B="#E15759", C="#F28E2B", D="#59A14F"); GREY <- "grey82"
TOTAL_TIPS <- 32

relabel_to_ref <- function(ref, x) {
  ref <- as.integer(ref); x <- as.integer(x); K <- 4
  conf <- matrix(0L, K, K); for (i in seq_along(ref)) conf[x[i], ref[i]] <- conf[x[i], ref[i]] + 1L
  map <- integer(K); used <- logical(K)
  for (idx in order(-conf)) { r <- ((idx-1)%%K)+1; cc <- ((idx-1)%/%K)+1
    if (map[r]==0 && !used[cc]) { map[r] <- cc; used[cc] <- TRUE } }
  for (r in which(map==0)) { cc <- which(!used)[1]; map[r] <- cc; used[cc] <- TRUE }
  out <- map[x]; names(out) <- names(x); out
}
floating_pie <- function(x, y, vals, radius, cols) {
  if (sum(vals) == 0) return(invisible()); vals <- vals/sum(vals); ang <- cumsum(c(0,vals))*2*pi + pi/2
  for (i in seq_along(vals)) { if (vals[i] <= 0) next
    th <- seq(ang[i], ang[i+1], length.out = max(2, ceiling((ang[i+1]-ang[i])/0.06)))
    polygon(c(x, x+radius*cos(th)), c(y, y+radius*sin(th)), col=cols[i], border="white", lwd=0.4) }
}

build_sex <- function(sex) {
  main <- readRDS(file.path(cl, sprintf("cluster_assignments_g1_%s_k4.rds", sex)))
  Z <- main$Z_active; if (is.null(rownames(Z))) rownames(Z) <- names(main$final_labels)
  ref <- as.integer(main$final_labels[rownames(Z)]); names(ref) <- rownames(Z); n <- nrow(Z)
  sizes <- table(factor(ref, levels=1:4)); Ks <- pmax(3, round(TOTAL_TIPS*as.numeric(sizes)/n))
  term_of <- setNames(rep(NA_integer_, n), rownames(Z)); SB <- list(); g <- 0
  for (s in 1:4) {
    idx <- rownames(Z)[ref==s]; Zs <- Z[idx,,drop=FALSE]
    hc <- hclust(dist(Zs),"ward.D2"); tk <- cutree(hc, k=Ks[s])
    cent <- t(sapply(1:Ks[s], function(j) colMeans(Zs[tk==j,,drop=FALSE])))
    st <- as.phylo(hclust(dist(cent),"ward.D2")); st$tip.label <- sprintf("K%d_t%d", s, 1:Ks[s])
    term_of[idx] <- g + tk; SB[[s]] <- list(tree=st, base=g, K=Ks[s]); g <- g + Ks[s]
  }
  gc <- t(sapply(1:4, function(s) colMeans(Z[ref==s,,drop=FALSE])))
  tree <- as.phylo(hclust(dist(gc),"ward.D2")); tree$tip.label <- paste0("K",1:4)
  for (s in 1:4) { pos <- which(tree$tip.label==paste0("K",s)); tree <- bind.tree(tree, SB[[s]]$tree, where=pos) }
  tree <- reorder(tree, "cladewise"); nT <- length(tree$tip.label)
  tip_sub <- as.integer(sub("K(\\d)_.*","\\1",tree$tip.label))
  tip_gid <- sapply(tree$tip.label, function(tl){ s<-as.integer(sub("K(\\d)_.*","\\1",tl)); j<-as.integer(sub(".*_t(\\d+)","\\1",tl)); SB[[s]]$base+j })
  tip_n <- sapply(tip_gid, function(gg) sum(term_of==gg, na.rm=TRUE))
  children <- split(tree$edge[,2], tree$edge[,1])
  desc_tips <- function(nd) if (nd<=nT) nd else unlist(lapply(children[[as.character(nd)]], desc_tips))
  edge_col <- character(nrow(tree$edge)); edge_w <- numeric(nrow(tree$edge))
  for (i in seq_len(nrow(tree$edge))) { tips <- desc_tips(tree$edge[i,2]); ss <- unique(tip_sub[tips])
    edge_col[i] <- if (length(ss)==1) PAL[SUBTYPE_MAP[[sex]][as.character(ss)]] else "grey75"; edge_w[i] <- sum(tip_n[tips]) }
  edge_lwd <- 0.8 + 2.4*(log1p(edge_w)/log1p(max(edge_w)))
  list(sex=sex, tree=tree, nT=nT, tip_sub=tip_sub, tip_gid=tip_gid, tip_n=tip_n,
       term_of=term_of, ref=ref, n=n, edge_col=edge_col, edge_lwd=edge_lwd)
}

twig_pie <- function(D, kind) {
  sex <- D$sex; ref <- D$ref; eid <- names(ref)
  if (kind == "seed") {
    fs <- list.files(cl, pattern=sprintf("cluster_assignments_g1_%s_reseed[0-9]+_k4\\.rds$", sex), full.names=TRUE)
    # Twenty re-seeds per sex, 10001-10020, as the Methods state. The summary file is
    # an append log, so it must be filtered by seed rather than read wholesale:
    # archived; female seeds 10021-10023 are valid but excluded to keep twenty per sex.
    # The retention arithmetic matches subtype_reproducibility/04_design_ari.R.
    sd_of <- as.integer(sub(".*reseed([0-9]+)_k4\\.rds$", "\\1", fs))
    fs <- fs[sd_of >= 10001 & sd_of <= 10020]
    stopifnot(length(fs) == 20)
    keptfrac <- rep(0, length(ref)); names(keptfrac) <- eid
    for (f in fs) { lb <- readRDS(f)$final_labels; rl <- relabel_to_ref(ref, lb[eid]); keptfrac <- keptfrac + (rl==ref) }
    keptfrac <- keptfrac/length(fs)
    cm <- sapply(seq_len(D$nT), function(t){ sel <- which(D$term_of==D$tip_gid[t]); c(sum(keptfrac[sel]), sum(1-keptfrac[sel])) })
    attr(cm,"kept") <- mean(keptfrac); return(cm)
  }
  lb <- readRDS(file.path(cl, sprintf("cluster_assignments_g1_%s_%s_k4.rds", sex, kind)))$final_labels
  rl <- relabel_to_ref(ref, lb[eid]); kept <- rl==ref
  cm <- sapply(seq_len(D$nT), function(t){ sel <- which(D$term_of==D$tip_gid[t]); c(sum(kept[sel]), sum(!kept[sel])) })
  attr(cm,"kept") <- mean(kept); cm
}

CELLS <- list(
  list(kind="seed", hdr="Seed re-runs"),
  list(kind="balanced",           hdr="Balanced weighting"),
  list(kind="broadened",          hdr="Broadened 6-view"),
  list(kind="broadened_balanced", hdr="Broadened + balanced"))

draw_one <- function(D, cell, header=NULL) {
  cm <- twig_pie(D, cell$kind); kept <- attr(cm,"kept")
  par(mar=c(0.2, 0.4, if (!is.null(header)) 2.8 else 0.2, 0.4), xpd=NA)
  #          bottom, left, top (header row) / top (no header), right;
  #          the no-header top margin is tightened to close the gap between the two sex rows
  pdf(NULL); plot.phylo(D$tree, type="fan", use.edge.length=FALSE, show.tip.label=FALSE,
                        edge.color=D$edge_col, edge.width=D$edge_lwd, open.angle=8)
  env <- get("last_plot.phylo", envir=ape::.PlotPhyloEnv); dev.off()
  xx <- env$xx[1:D$nT]; yy <- env$yy[1:D$nT]; span <- max(abs(c(env$xx, env$yy)))
  rad <- 0.085*span*sqrt(D$tip_n/max(D$tip_n)) + 0.03*span
  # A wider lim draws the tree smaller with more surrounding white space; a narrower one
  # fills the panel.
  lim <- c(-(span + max(rad) * 1.05), span + max(rad) * 1.05)
  plot.phylo(D$tree, type="fan", use.edge.length=FALSE, show.tip.label=FALSE,
             edge.color=D$edge_col, edge.width=D$edge_lwd, open.angle=8, x.lim=lim, y.lim=lim)
  for (t in seq_len(D$nT)) floating_pie(xx[t], yy[t], cm[,t], rad[t], c(PAL[SUBTYPE_MAP[[D$sex]][as.character(D$tip_sub[t])]], GREY))
  if (!is.null(header)) mtext(header, side=3, line=0.6, cex=2.2, font=2)  # font-harmonised
  # kept% in the centre "donut hole" on a white disc, large
  rc <- 0.36 * span; thc <- seq(0, 2*pi, length.out = 80)
  polygon(rc*cos(thc), rc*sin(thc), col = "white", border = NA)
  text(0, 0.12*span, sprintf("%.0f%%", 100*kept), cex = 2.4, font = 2,
       col = if (kept >= 0.7) "grey10" else "#B22222")
  # "retained" is a wide label; a smaller cex keeps it inside the
  # white centre circle rather than spilling onto the tree
  text(0, -0.21*span, "retained", cex = 1.5, font = 2, col = "grey20")
}

draw_fig <- function() {
  DF <- build_sex("f"); DM <- build_sex("m")
  # 4 tree columns + a right column (cell 9) for the vertical legend, spanning both rows
  layout(matrix(c(1,2,3,4,9, 5,6,7,8,9), 2, 5, byrow=TRUE), widths=c(1,1,1,1,0.42))
  # subtype colour key = SHARED bottom legend of the whole Fig 5; A keeps only "reassigned"
  for (i in seq_along(CELLS)) draw_one(DF, CELLS[[i]], header=CELLS[[i]]$hdr)
  for (i in seq_along(CELLS)) draw_one(DM, CELLS[[i]], header=NULL)
  par(mar=c(0,0.4,0,0.4)); plot.new()
  legend("center", bty="n", cex=1.9, pt.cex=3.2, pch=21, col="white", y.intersp=1.8,
         pt.bg=GREY, legend="Reassigned")
  mtext("Female", outer=TRUE, side=2, line=0.6, at=0.72, las=0, cex=2.1, font=2)
  mtext("Male",   outer=TRUE, side=2, line=0.6, at=0.27, las=0, cex=2.1, font=2)
  # eBM: descriptive title moved to figure legend; panel letter only
  mtext("A", outer=TRUE, side=3, line=0.6, adj=0, cex=2.7, font=2)
}

png(file.path(outdir, "Fig5_panelA_radial_8panel.png"), width=4300, height=2250, res=200)
par(oma=c(0,2.8,3.0,0.4)); draw_fig(); dev.off()
pdf(file.path(outdir, "Fig5_panelA_radial_8panel.pdf"), width=21.5, height=11.2)
par(oma=c(0,2.8,3.0,0.4)); draw_fig(); dev.off()
cat("DONE panelA ->", outdir, "\n")
