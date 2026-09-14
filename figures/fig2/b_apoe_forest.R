#> in   input/assoc_results/apoe_diplotype_association.csv
#> out  output/figures/fig2b_apoe_forest.pdf
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
source("figures/_common/style.R")
suppressPackageStartupMessages(library(patchwork))

APOE_F <- "input/assoc_results/apoe_diplotype_association.csv"
stopifnot(file.exists(APOE_F))

## ---- 1. load ---------------------------------------------------------------
d <- fread(APOE_F, encoding = "UTF-8")
d[, intersect(c("label", "OR_fmt", "CI_fmt", "P_fmt"), names(d)) := NULL]
setnames(d, c("Comparison", "CI_lower", "CI_upper", "P"),
         c("label", "or_lower", "or_upper", "pval"), skip_absent = TRUE)
setnames(d, "OR", "or", skip_absent = TRUE)
d <- d[label != "---"]
d[, label := gsub("\u03b5", "e", label)]     # \u escape keeps this source ASCII
stopifnot(!any(grepl("[^ -~]", d$label)))    # every epsilon must have mapped

relabel <- c("Per e4 allele (trend)"     = "Per e4 allele",
             "e4 carrier vs non-carrier" = "e4 carrier",
             "e3/e3 (reference)"         = "e3/e3 (ref)")
d[label %in% names(relabel), label := relabel[label]]

ROW_ORDER <- c("Per e4 allele", "e4 carrier",
               "e4/e4", "e3/e4", "e2/e4", "e2/e3", "e2/e2", "e3/e3 (ref)")
d[, ord := match(label, ROW_ORDER)]
stopifnot(!any(is.na(d$ord)), nrow(d) == 8)  # no row silently unordered or lost
setorder(d, ord)

d[, block := c(rep("Allele dose", 2), rep("Diplotype vs e3/e3", 6))]
d[, block := factor(block, levels = c("Allele dose", "Diplotype vs e3/e3"))]

# reference row carries no interval
d[, is_ref := label == "e3/e3 (ref)"]
d[, or_txt := fifelse(is_ref, "1.00 (ref)",
                      sprintf("%.2f (%.2f-%.2f)", or, or_lower, or_upper))]
# same encoding as panels c and d: filled = significant, open = not
d[, mark := fcase(is_ref,                        "ref",
                  pval < 0.05 & or > 1,          "up",
                  pval < 0.05 & or < 1,          "down",
                  default =                      "ns")]
d[, mark := factor(mark, levels = c("up", "down", "ns", "ref"))]
stopifnot(max(nchar(d$label)) <= 25)
cat(sprintf("[rows] %d | significant at P<0.05: %d\n",
            nrow(d), sum(d$mark %in% c("up", "down"))))

## ---- 2. y coordinates ------------------------------------------------------
d[, blk_id := as.integer(block)]
d[, y := rev(seq_len(.N)) + (max(blk_id) - blk_id) * 1.0]
d[, blk_lab := fifelse(!duplicated(block), as.character(block), NA_character_)]
Y_TOP <- max(d$y) + 2.5; Y_BOT <- 0.35

## ---- 3. layout -------------------------------------------------------------
PITCH <- 4.089
W_MM  <- 81.5
H_MM  <- round((Y_TOP - Y_BOT) * PITCH, 1)     # -> 45.6
cat(sprintf("[layout] %.2f y-units x %.3f mm = %.1f mm tall\n",
            Y_TOP - Y_BOT, PITCH, H_MM))

XLIM <- c(0.78, 1.75)
XBRK <- c(0.8, 1.0, 1.25, 1.5)
SZA  <- SZ$anno / .pt

SHAPE <- c(up = 16, down = 16, ns = 21, ref = 21)
COLR  <- c(up = "#B23A48", down = "#4E79A7", ns = "grey50", ref = "grey30")
FILLC <- c(up = "#B23A48", down = "#4E79A7", ns = "white", ref = "grey30")

# identical to fig2/c_protein_forest.R -- do not change one without the other
X_LAB <- 0.00; X_OR <- 1.05; X_R <- 1.72

p_left <- ggplot(d) +
  annotate("text", x = X_LAB, y = Y_TOP - 0.45, label = "APOE genotype",
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
  geom_errorbar(data = d[!is_ref & !is.na(or_lower)],
                aes(xmin = or_lower, xmax = or_upper, colour = mark),
                orientation = "y", width = 0.42, linewidth = LW$err) +
  geom_point(aes(shape = mark, colour = mark, fill = mark), size = 1.3, stroke = 0.45) +
  scale_shape_manual(values = SHAPE, guide = "none") +
  scale_colour_manual(values = COLR, guide = "none") +
  scale_fill_manual(values = FILLC, guide = "none") +
  scale_x_continuous(transform = "log", limits = XLIM, breaks = XBRK,
                     labels = function(x) format(x, drop0trailing = TRUE), expand = c(0, 0)) +
  scale_y_continuous(limits = c(Y_BOT, Y_TOP), expand = c(0, 0)) +
  labs(x = "Odds ratio (0.78-1.75)") +
  theme_fig(grid = "none") +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        plot.margin = margin(1.5, 1.5, 1.5, 1.5, "mm"))

p <- p_left + p_mid + plot_layout(widths = c(X_R, 1.00))
fig_save(p, "fig2b_apoe_forest", W_MM, H_MM)

cat("\n[rows]\n"); print(d[, .(label, or_txt, mark)])
cat("\n[legend] Odds ratios are per e4 allele (trend), for e4 carriage, and for each\n")
cat("         diplotype against e3/e3. Filled = P < 0.05, open = not significant.\n")
cat("[legend] P values are in the appendix table; the genome-wide P for rs429358 is\n")
cat("         labelled in panel a.\n")
cat("[legend] NOTE the x-axis range (0.78-1.75) differs from panels c and d.\n")
