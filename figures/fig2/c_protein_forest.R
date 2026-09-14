#> in   input/assoc_results/pwas_primary.csv
#> out  output/figures/fig2c_protein_forest.pdf
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
source("figures/_common/style.R")
suppressPackageStartupMessages(library(patchwork))

## ---- 1. load + sanity gate -------------------------------------------------
g1 <- fread("input/assoc_results/pwas_primary.csv")
cat(sprintf("[load] PWAS main model: %d proteins\n", nrow(g1)))

N_FDR <- sum(g1$sig_fdr); N_BONF <- sum(g1$sig_bonf); N_TOT <- nrow(g1)
N_ANA <- unique(g1$n_total)
cat(sprintf("[check] FDR %d (17) | Bonferroni %d (4) | proteins %d (2920) | n %s (15098)\n",
            N_FDR, N_BONF, N_TOT, paste(N_ANA, collapse = "/")))
stopifnot(N_FDR == 17, N_BONF == 4, N_TOT == 2920, N_ANA == 15098)

## ---- 2. rows: every FDR-significant protein, no selection -------------------
NAMES <- c(
  a1bg    = "Alpha-1B-glycoprotein",
  f9      = "Coagulation factor IX",
  orm1    = "Orosomucoid 1",
  cfb     = "Complement factor B",
  apcs    = "Serum amyloid P",
  palm2   = "Paralemmin-2",
  inhbc   = "Inhibin beta C",
  ctsd    = "Cathepsin D",
  igsf9   = "Ig superfamily member 9",
  ces1    = "Carboxylesterase 1",
  lep     = "Leptin",
  oxt     = "Oxytocin-neurophysin 1",
  ccl19   = "C-C chemokine 19",
  cntn5   = "Contactin 5",
  wfikkn2 = "WFIKKN2",
  enpp6   = "Ectonucleotide PPase 6",
  pon3    = "Paraoxonase 3")

d <- g1[sig_fdr == TRUE, .(protein, or, or_lower, or_upper, pval, sig_bonf)]
stopifnot(nrow(d) == 17, all(d$protein %in% names(NAMES)))
d[, label := NAMES[protein]]
# Bonferroni-significant rows carry a trailing marker, explained in the legend
d[, label := fifelse(sig_bonf, paste0(label, " *"), label)]
d[, block := fifelse(or > 1, "Higher in cases", "Lower in cases")]
# order: raised first, each block by descending effect
setorder(d, -or)
d[, block := factor(block, levels = c("Higher in cases", "Lower in cases"))]
setorder(d, block, -or)
cat(sprintf("[rows] %d raised | %d reduced\n", sum(d$or > 1), sum(d$or < 1)))
stopifnot(max(nchar(d$label)) <= 25)

d[, or_txt := sprintf("%.2f (%.2f-%.2f)", or, or_lower, or_upper)]
d[, mark := factor(fifelse(or > 1, "up", "down"), levels = c("up", "down"))]

## ---- 3. y coordinates ------------------------------------------------------
d[, blk_id := as.integer(block)]
d[, y := rev(seq_len(.N)) + (max(blk_id) - blk_id) * 1.0]
d[, blk_lab := fifelse(!duplicated(block), as.character(block), NA_character_)]
Y_TOP <- max(d$y) + 2.5; Y_BOT <- 0.35

## ---- 4. layout -------------------------------------------------------------
PITCH <- 4.089
W_MM <- 81.5; H_MM <- round((Y_TOP - Y_BOT) * PITCH, 1)
XLIM <- c(0.40, 5.50)
XBRK <- c(0.5, 1, 2, 4)
SZA  <- SZ$anno / .pt

COLR  <- c(up = "#B23A48", down = "#4E79A7")

X_LAB <- 0.00; X_OR <- 1.05; X_R <- 1.72

p_left <- ggplot(d) +
  annotate("text", x = X_LAB, y = Y_TOP - 0.45, label = "Plasma protein",
           hjust = 0, size = SZA, family = FIG_FONT, fontface = "bold") +
  annotate("text", x = X_OR, y = Y_TOP - 0.45, label = "OR (95% CI)",
           hjust = 0, size = SZA, family = FIG_FONT, fontface = "bold") +
  annotate("segment", x = X_LAB, xend = X_R, y = Y_TOP - 1.05, yend = Y_TOP - 1.05,
           colour = "grey35", linewidth = 0.25) +
  geom_text(data = d[!is.na(blk_lab)], aes(x = X_LAB, y = y + 0.95, label = blk_lab),
            hjust = 0, size = SZA, family = FIG_FONT, fontface = "italic", colour = "grey30") +
  geom_text(aes(x = X_LAB + 0.045, y = y, label = label), hjust = 0,
            size = SZA, family = FIG_FONT, colour = "black") +
  geom_text(aes(x = X_OR, y = y, label = or_txt), hjust = 0,
            size = SZA, family = FIG_FONT, colour = "black") +
  scale_x_continuous(limits = c(-0.01, X_R), expand = c(0, 0)) +
  scale_y_continuous(limits = c(Y_BOT, Y_TOP), expand = c(0, 0)) +
  theme_void(base_family = FIG_FONT) +
  theme(plot.margin = margin(1.5, 0, 1.5, 1.5, "mm"))

p_mid <- ggplot(d, aes(x = or, y = y)) +
  annotate("segment", x = 1, xend = 1, y = Y_BOT, yend = Y_TOP - 1.05,
           colour = "grey45", linewidth = LW$ref) +
  geom_errorbar(aes(xmin = or_lower, xmax = or_upper, colour = mark),
                orientation = "y", width = 0.42, linewidth = LW$err) +
  geom_point(aes(colour = mark), shape = 16, size = 1.3) +
  scale_colour_manual(values = COLR, guide = "none") +
  scale_x_continuous(transform = "log", limits = XLIM, breaks = XBRK,
                     labels = function(x) format(x, drop0trailing = TRUE), expand = c(0, 0)) +
  scale_y_continuous(limits = c(Y_BOT, Y_TOP), expand = c(0, 0)) +
  labs(x = "Odds ratio per SD (0.4-5.5)") +
  theme_fig(grid = "none") +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        plot.margin = margin(1.5, 1.5, 1.5, 1.5, "mm"))

p <- p_left + p_mid + plot_layout(widths = c(X_R, 1.00))
fig_save(p, "fig2c_protein_forest", W_MM, H_MM)

cat("\n[rows]\n"); print(d[, .(label, or_txt, block)])

cat(sprintf("\n[legend] All %d of %d plasma proteins significant at 5%% FDR are shown;\n",
            N_FDR, N_TOT))
cat(sprintf("         * marks the %d that also pass Bonferroni.\n", N_BONF))
cat("[legend] NOTE the x-axis range (0.4-5.5) is far wider than the metabolite panel\n")
cat("         (0.80-1.26); effect sizes are NOT comparable between the two panels.\n")
cat("[note ]  Full 2,920-protein volcano -> appfig24_olink_volcano.\n")
