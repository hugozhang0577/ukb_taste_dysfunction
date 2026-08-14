#!/usr/bin/env Rscript
# =============================================================================
# APOE diplotype — publication forest plot (left table + right forest)
# =============================================================================
# Draws the APOE association forest from 02_apoe_association.R output: the
# 6-level diplotype (e3/e3 reference) plus the e4-carrier and per-allele trend
# rows. This is the working figure behind main-text Fig 3d (the locked composite
# is assembled separately and uses ASCII-only glyphs for cairo_pdf).
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
APOE_DIR <- file.path(PROJECT_DIR, "gwas/cohort_primary/genotype/APOE_type")
FIG_DIR <- file.path(APOE_DIR, "Figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

all_results <- fread(file.path(APOE_DIR, "Results", "plot_data_association_results.csv"))

display_order <- c("e3/e3 (reference)", "e2/e2", "e2/e3", "e2/e4", "e3/e4", "e4/e4",
                   "e4 carrier vs non-carrier", "Per e4 allele (trend)")

fd <- all_results[Comparison != "---"]
fd[, is_reference := Comparison == "e3/e3 (reference)"]
fd[, color_group := fcase(is_reference | is.na(P), "reference",
                          OR > 1 & P < 0.05, "risk",
                          OR < 1 & P < 0.05, "protect",
                          default = "ns")]
fd[, Comparison := factor(Comparison, levels = display_order)]
setorder(fd, Comparison)
fd[, row_id := .N:1]
fd[, p_label := fifelse(is_reference, "-",
                fifelse(P < 0.001, sprintf("%.2e", P), sprintf("%.3f", P)))]
fd[, or_ci_label := fifelse(is_reference, "Reference",
                            sprintf("%.2f (%.2f-%.2f)", OR, CI_lower, CI_upper))]

valid <- fd[!is_reference & !is.na(CI_lower)]
x_max <- max(valid$CI_upper) * 1.1
y_max <- max(fd$row_id) + 1
table_width <- 2.2
col_gene <- -table_width + 0.05; col_pval <- -table_width + 1.2; col_or <- -table_width + 1.7
fs <- 4.2
cols <- c(risk = "#E64B35", protect = "#00A087", ns = "gray50", reference = "gray40")

p <- ggplot(fd) +
  annotate("segment", x = -table_width, xend = 0, y = y_max - 0.3, yend = y_max - 0.3, linewidth = 0.8) +
  annotate("segment", x = -table_width, xend = 0, y = y_max - 0.7, yend = y_max - 0.7, linewidth = 0.8) +
  annotate("segment", x = -table_width, xend = 0, y = 0.3, yend = 0.3, linewidth = 0.8) +
  annotate("text", x = col_pval, y = y_max - 0.5, label = "P value", hjust = 0.5, size = fs + 0.5, fontface = "bold") +
  annotate("text", x = col_or + 0.25, y = y_max - 0.5, label = "Odds Ratio", hjust = 0, size = fs + 0.5, fontface = "bold") +
  geom_text(aes(col_gene, row_id, label = Comparison), hjust = 0, size = fs) +
  geom_text(aes(col_pval, row_id, label = p_label), hjust = 0.5, size = fs) +
  geom_text(aes(col_or, row_id, label = or_ci_label), hjust = 0, size = fs) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray30", linewidth = 0.6) +
  geom_segment(data = fd[!is_reference == TRUE], aes(CI_lower, row_id, xend = CI_upper, yend = row_id), linewidth = 0.7) +
  geom_segment(data = fd[!is_reference == TRUE], aes(CI_lower, row_id - 0.2, xend = CI_lower, yend = row_id + 0.2), linewidth = 0.7) +
  geom_segment(data = fd[!is_reference == TRUE], aes(CI_upper, row_id - 0.2, xend = CI_upper, yend = row_id + 0.2), linewidth = 0.7) +
  geom_point(data = fd[!is_reference == TRUE], aes(OR, row_id, color = color_group), size = 4) +
  scale_color_manual(values = cols, guide = "none") +
  scale_x_continuous(limits = c(-table_width, x_max), breaks = seq(0.8, 1.8, 0.2),
                     expand = c(0, 0), oob = scales::oob_keep) +
  scale_y_continuous(limits = c(0, y_max), expand = c(0.01, 0.01)) +
  labs(x = "Odds Ratio", y = NULL) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 14) +
  theme(axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 8), hjust = 0.7),
        axis.text.x = element_text(size = 12, color = "black", margin = margin(t = 5)),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.ticks.x = element_line(color = "black", linewidth = 0.5),
        axis.ticks.length.x = unit(0.15, "cm"),
        plot.margin = margin(20, 25, 15, 15))

ggsave(file.path(FIG_DIR, "APOE_forest_plot.png"), p, width = 10, height = 6, dpi = 300, bg = "white")
ggsave(file.path(FIG_DIR, "APOE_forest_plot.pdf"), p, width = 10, height = 6)
cat("APOE forest plot ->", FIG_DIR, "\n")
