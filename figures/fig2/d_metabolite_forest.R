#> in   input/assoc_results/mwas_primary.csv
#> in   input/reference/nmr_metabolite_annotation.csv
#> out  output/figures/fig2d_metabolite_forest.pdf
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
source("figures/_common/style.R")
suppressPackageStartupMessages(library(patchwork))

## ---- 1. load + sanity gate -------------------------------------------------
g1  <- fread("input/assoc_results/mwas_primary.csv")
ann <- fread("input/reference/nmr_metabolite_annotation.csv")
cat(sprintf("[load] G1 %d measures | annotation %d rows\n", nrow(g1), nrow(ann)))

N_FDR <- sum(g1$sig_fdr); N_BONF <- sum(g1$sig_bonf); N_TOT <- nrow(g1)
cat(sprintf("[check] FDR %d (expect 222) | Bonferroni %d (expect 151) | total %d (expect 327)\n",
            N_FDR, N_BONF, N_TOT))
stopifnot(N_FDR == 222, N_BONF == 151, N_TOT == 327)

# the annotation table must describe the SAME model, or the categories are stale
chk <- merge(ann[, .(protein = metabolite, or_a = or)], g1[, .(protein, or)], by = "protein")
stopifnot(nrow(chk) == N_TOT, max(abs(chk$or_a - chk$or)) < 1e-9)
cat("[check] annotation matches N7 exactly (327/327, max abs dOR < 1e-9)\n")

## ---- 2. rows ---------------------------------------------------------------
ROWS <- data.table(
  protein = c(
    # -- apolipoproteins and lipids
    "ApoA1", "ApoB_by_ApoA1", "HDL_C", "VLDL_C", "Total_TG",
    # -- particle size
    "HDL_size", "VLDL_size", "LDL_size",
    # -- lipoprotein subclasses
    "XL_HDL_PL", "XL_VLDL_C", "M_HDL_TG", "S_HDL_L",
    # -- composition: TG-up / cholesterol-down axis, then esterification axis
    "XL_HDL_TG_pct", "L_HDL_C_pct", "HDL_TG_pct", "L_HDL_PL_pct",
    "XL_HDL_FC_by_CE", "L_VLDL_FC_pct_C", "XL_HDL_CE_pct_C", "M_VLDL_FC_pct_C",
    # -- fatty acids
    "MUFA_pct", "PUFA_pct", "Unsaturation", "MUFA",
    # -- amino acids and other metabolites
    "Gln", "GlycA", "Albumin", "Citrate",
    # -- negative controls
    "Total_C", "LDL_C", "Glucose"),
  label = c(
    "Apolipoprotein A1", "ApoB / ApoA1", "HDL cholesterol", "VLDL cholesterol",
    "Total triglycerides",
    "HDL particle size", "VLDL particle size", "LDL particle size",
    "Phospholipids, XL HDL", "Cholesterol, XL VLDL", "Triglycerides, M HDL",
    "Total lipids, S HDL",
    "Triglycerides, XL HDL %", "Cholesterol, L HDL %",
    "Triglycerides, all HDL %", "Phospholipids, L HDL %",
    "Free chol/ester, XL HDL", "Free chol %, L VLDL",
    "Chol ester %, XL HDL", "Free chol %, M VLDL",
    "Monounsaturated, %", "Polyunsaturated, %",
    "Degree of unsaturation", "Monounsaturated, total",
    "Glutamine", "Glycoprotein acetyls", "Albumin", "Citrate",
    "Total cholesterol", "LDL cholesterol", "Glucose"),
  block = c(rep("Apolipoproteins and lipids", 5),
            rep("Lipoprotein particle size", 3),
            rep("Lipoprotein subclasses", 4),
            rep("Particle composition", 8),
            rep("Fatty acids", 4),
            rep("Amino acids and other metabolites", 4),
            rep("Not associated", 3)))
stopifnot(nrow(ROWS) == 31)

d <- merge(ROWS, g1[, .(protein, or, or_lower, or_upper, pval, sig_fdr, sig_bonf)],
           by = "protein", all.x = TRUE, sort = FALSE)
if (anyNA(d$or)) stop("missing in G1: ", paste(d[is.na(or), protein], collapse = ", "))
cat(sprintf("[join] %d/%d rows matched\n", sum(!is.na(d$or)), nrow(ROWS)))

# the negative-control block must actually BE null, or the claim is wrong
stopifnot(!any(d[block == "Not associated", sig_fdr]))
cat("[check] negative-control block is null in every row\n")

d[, or_txt := sprintf("%.2f (%.2f-%.2f)", or, or_lower, or_upper)]
d[, mark := fcase(sig_fdr & or > 1, "up", sig_fdr & or < 1, "down", default = "ns")]
d[, mark := factor(mark, levels = c("up", "down", "ns"))]

## ---- 3. y coordinates: gap between blocks ----------------------------------
d[, blk_id := as.integer(factor(block, levels = unique(block)))]
d[, y := rev(seq_len(.N)) + (max(blk_id) - blk_id) * 1.0]
d[, blk_lab := fifelse(!duplicated(block), block, NA_character_)]
Y_TOP <- max(d$y) + 2.5; Y_BOT <- 0.35

## ---- 4. layout: one journal column wide, 3 columns of content ---------------
W_MM <- 81.5; H_MM <- 132
XLIM <- c(0.80, 1.26)
XBRK <- c(0.8, 0.9, 1.0, 1.1, 1.2)
SZA  <- SZ$anno / .pt

SHAPE <- c(up = 16, down = 16, ns = 21)
COLR  <- c(up = "#B23A48", down = "#4E79A7", ns = "grey50")
FILLC <- c(up = "#B23A48", down = "#4E79A7", ns = "white")

X_LAB <- 0.00; X_OR <- 1.05; X_R <- 1.72

p_left <- ggplot(d) +
  annotate("text", x = X_LAB, y = Y_TOP - 0.45, label = "Metabolic measure",
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
  geom_point(aes(shape = mark, colour = mark, fill = mark), size = 1.3, stroke = 0.45) +
  scale_shape_manual(values = SHAPE, guide = "none") +
  scale_colour_manual(values = COLR, guide = "none") +
  scale_fill_manual(values = FILLC, guide = "none") +
  scale_x_continuous(transform = "log", limits = XLIM, breaks = XBRK,
                     labels = function(x) format(x, drop0trailing = TRUE), expand = c(0, 0)) +
  scale_y_continuous(limits = c(Y_BOT, Y_TOP), expand = c(0, 0)) +
  labs(x = "Odds ratio per SD (0.80-1.26)") +
  theme_fig(grid = "none") +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        plot.margin = margin(1.5, 1.5, 1.5, 1.5, "mm"))

p <- p_left + p_mid + plot_layout(widths = c(X_R, 1.00))
fig_save(p, "fig2d_metabolite_forest", W_MM, H_MM)

cat("\n[rows]\n"); print(d[, .(label, or_txt, mark)])

cat(sprintf("\n[legend] %d of %d NMR measures were significant at 5%% FDR and %d at Bonferroni.\n",
            N_FDR, N_TOT, N_BONF))
cat("[legend] Rows are one per de-redundancy cluster, grouped by Nightingale category.\n")
cat("[legend] HDL cholesterol is collinear with HDL particle size, and total triglycerides\n")
cat("         with VLDL particle size; both pairs are shown and each pair is one cluster.\n")
cat("[legend] Open circles = not significant. Total-C, LDL-C and Glucose are retained as\n")
cat("         negative controls: the signal is in particle composition and inflammation,\n")
cat("         not in the standard lipid panel.\n")
cat("[legend] NOTE the x-axis range (0.80-1.26) differs from the Olink panel (0.4-5.5).\n")
cat("[note ]  Full 327-measure table with P values and held-out concordance -> appendix.\n")
