# =============================================================================
# Fig 5 Panels D and E  --  subtype genetics, COVID exposure and severity
# Renders two panels from one script, because they share the subtype-name axis:
#   D : HLA per-subtype forest (lead rs2071293)
#   E : COVID pre-taste OR, crude vs age+sex (null)  |  severity grouped bar
#       (taste/smell chronicity + functional impact + smell co-occurrence,
#        5 metrics x 4 subtypes, all proportions)
# The two forests share one subtype-name axis; names are shown on the HLA forest
# only. The ASSET Manhattan plot is a supplementary figure, not part of Figure 5.
# No in-graph descriptive titles (eBM); panel letter only; interpretation -> legend.
#
# HLA input: per-subtype SAIGE association at the ASSET lead variant rs2071293
# (chr6:32,062,687), effect allele A (SAIGE Allele2). ASSET joint p = 3.01e-7,
# positive subset {A,B,C} vs {D}. APOE and CETP are excluded: both are lipid loci
# whose genotype enters the subtyping through the lipid features, so a
# subtype-genotype association at either is circular by construction.
# =============================================================================
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
setDTthreads(1)

out <- "output/figures/fig5/panels"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

SUB_PAL <- c(A = "#4E79A7", B = "#E15759", C = "#F28E2B", D = "#59A14F")
NAMES   <- c(A = "Ageing frailty", B = "Psychosomatic",
             C = "Cardiometabolic", D = "Young idiopathic")
ord <- c("A", "B", "C", "D")
wilson <- function(k, n, z = 1.96) {
  if (n == 0 || is.na(n)) return(c(NA_real_, NA_real_))
  p <- k / n; d <- 1 + z^2 / n; ctr <- (p + z^2 / (2*n)) / d
  hw <- z * sqrt(p*(1-p)/n + z^2/(4*n^2)) / d; c(ctr - hw, ctr + hw) * 100
}
mk_bands <- function(xmin, xmax) {
  b <- data.table(yi = 1:4)
  b[, `:=`(ymin = yi - 0.5, ymax = yi + 0.5, xmin = xmin, xmax = xmax,
           fill = ifelse(yi %% 2 == 1, "grey92", "white"))]; b[]
}
thmD <- theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 9, colour = "black", margin = margin(t = 5)),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10, colour = "black", hjust = 0),
        axis.text.x = element_text(size = 8.5, colour = "black"),
        legend.position = "bottom", legend.title = element_text(size = 8),
        legend.text = element_text(size = 7.5), legend.key.size = unit(3.4, "mm"),
        plot.margin = margin(6, 10, 6, 6))

# Panel D (HLA forest) is composited at ~0.64x (in the B|[C/D] row) whereas Panel E
# is at ~1.0x; bump D's fonts ~1.5x so both render at matched final sizes.
thmD_D <- thmD + theme(axis.title.x = element_text(size = 14, colour = "black", margin = margin(t = 5)),
                       axis.text.y  = element_text(size = 15, colour = "black", hjust = 0),
                       axis.text.x  = element_text(size = 13, colour = "black"))

# Per-participant COVID timing relative to the taste questionnaire, joined to the
# subtype label. Built from output/subtyping/clusters/cluster_assignments_g1_{m,f}_k4.rds
# and the subtype map, so it must be rebuilt whenever the model is refitted.
# Subtype sizes at the reported fit: A 1837 / B 754 / C 1625 / D 1556 (n = 5772).
map <- fread("output/subtyping/evidence_covid/covid_subtype_labels.csv")
map <- map[subtype %in% ord]

# =============================================================================
# HLA forest: per-subtype association at the chr6 lead variant rs2071293
# (effect allele A). OR = exp(BETA), CI = exp(BETA +/- 1.96*SE), taken from the
# per-subtype genome-wide scans.
# =============================================================================
HLAF <- "output/subtyping/evidence_genetic/hla_persubtype_effects.csv"
if (!file.exists(HLAF))
  stop("missing HLA per-subtype table -- run subtype_characterisation/",
       "06_hla_variant_extract.sh then 07_hla_persubtype.R: ", HLAF)
DAT <- fread(HLAF)[stratum == "all", .(sub = subtype, or = OR, lo, hi, p)]
if (nrow(DAT) != 4L) stop("expected 4 pooled subtype rows, got ", nrow(DAT))
DAT[, sig := p < 0.05]
DAT[, sub := factor(sub, levels = ord)]
DAT[, fillcol := SUB_PAL[as.character(sub)]]   # all SOLID (open/filled reserved for COVID crude/adj)

# =============================================================================
# COVID forest data -- the one-versus-rest ORs the Methods and Figure 5E report,
# read rather than recomputed here so the figure and the text cannot diverge.
# Produced by subtype_characterisation/05_covid_association.R.
# =============================================================================
COVOR <- "output/subtyping/evidence_covid/covid_onevsrest_OR.csv"
if (!file.exists(COVOR))
  stop("missing COVID OR table -- run subtype_characterisation/05_covid_association.R: ", COVOR)
cov <- fread(COVOR)[model %in% c("sex_only", "age_sex")]
cov[, model := factor(fifelse(model == "sex_only", "Crude (sex)", "Age + sex adj"),
                      levels = c("Crude (sex)", "Age + sex adj"))]
setnames(cov, "subtype", "sub")
cov[, sub := factor(sub, ord)]
cov <- cov[, .(sub, OR, lo, hi, P, model)]
# Reference values as reported in covid_onevsrest_OR.csv
# (sex_only 1.600; age_sex 0.904). The +/-0.02 tolerance is deliberately tight:
# it is the gate that catches a stale input, so a wider window would defeat it.
stopifnot(abs(cov[sub == "D" & model == "Crude (sex)", OR] - 1.60) < 0.02,
          abs(cov[sub == "D" & model == "Age + sex adj", OR] - 0.90) < 0.02)

# =============================================================================
# severity metrics (5, proportions) -> vertical grouped bar
# =============================================================================
jb <- fread("input/analysis_ready/base_table_full.csv",
            select = c("eid","taste_time","taste_extent",
                       "smell_change","smell_time","smell_extent"), showProgress = FALSE)
dT <- merge(map[, .(eid, subtype)], jb, by = "eid", all.x = TRUE)
for (c0 in c("taste_time","taste_extent","smell_change","smell_time","smell_extent"))
  dT[get(c0) < 0, (c0) := NA]
m1 <- dT[!is.na(taste_time),                     .(k=sum(taste_time==3),   n=.N), by=subtype][, metric:="Taste\n>12wk"]
m2 <- dT[smell_change==1 & !is.na(smell_time),   .(k=sum(smell_time==3),   n=.N), by=subtype][, metric:="Smell\n>12wk"]
m3 <- dT[!is.na(taste_extent),                   .(k=sum(taste_extent==1), n=.N), by=subtype][, metric:="Taste\nimpact"]
m4 <- dT[smell_change==1 & !is.na(smell_extent), .(k=sum(smell_extent==1), n=.N), by=subtype][, metric:="Smell\nimpact"]
m5 <- dT[!is.na(smell_change),                   .(k=sum(smell_change==1), n=.N), by=subtype][, metric:="Smell\nco-occur"]
sev <- rbindlist(list(m1,m2,m3,m4,m5))
sev[, pct := 100*k/n]
sev[, c("lo","hi") := as.data.table(t(mapply(wilson, k, n)))]
sev[, sub := factor(subtype, ord)]
sev[, metric := factor(metric, levels = c("Taste\n>12wk","Smell\n>12wk","Taste\nimpact","Smell\nimpact","Smell\nco-occur"))]
stopifnot(abs(sev[metric=="Smell\nimpact" & subtype=="B", pct] - 16.5) < 1.5)  # as reported

# =============================================================================
# PLOTS
# =============================================================================
## HLA forest (subtype NAMES shown here; shared with COVID forest) --------------
XB_HLA <- c(0.85, 1.0, 1.2)
p_hla <- ggplot(DAT, aes(or, sub)) +
  geom_rect(data = mk_bands(0.3, 3), inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill)) +
  scale_fill_identity() +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50", linewidth = .5) +
  geom_errorbarh(aes(xmin = lo, xmax = hi, colour = sub), height = .14, linewidth = .7) +
  geom_point(aes(colour = sub, fill = fillcol), shape = 21, size = 3.4, stroke = 1.2) +
  annotate("text", x = 0.955, y = 1, label = "opposite (n.s.)", vjust = 2.1, size = 2.4,
           fontface = 3, colour = "grey40") +
  scale_colour_manual(values = SUB_PAL, guide = "none") +
  scale_x_log10(breaks = XB_HLA, labels = sprintf("%.2f", XB_HLA)) +
  scale_y_discrete(limits = rev(ord)) +   # subtype letters A-D; names in shared legend
  coord_cartesian(xlim = c(0.80, 1.42)) +
  labs(x = "HLA per-allele OR (rs2071293)", y = NULL) + thmD_D

## COVID forest (y labels HIDDEN -> shares HLA's subtype names) -----------------
p_cov <- ggplot(cov, aes(OR, sub, colour = sub, shape = model, alpha = model)) +
  geom_rect(data = mk_bands(0.4, 2.6), inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill)) +
  scale_fill_identity() +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50", linewidth = .5) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = .6,
                 position = position_dodge(width = .58)) +
  geom_point(size = 2.4, position = position_dodge(width = .58)) +
  annotate("text", x = 1.15, y = 0.62, label = "D: age artefact (null)",
           size = 2.5, hjust = 0.5, colour = "grey25") +
  scale_shape_manual(values = c("Crude (sex)" = 1, "Age + sex adj" = 16), name = NULL) +
  scale_alpha_manual(values = c("Crude (sex)" = .55, "Age + sex adj" = 1), guide = "none") +
  scale_colour_manual(values = SUB_PAL, guide = "none") +
  scale_x_log10(breaks = c(0.5, 0.7, 1, 1.5, 2)) +
  scale_y_discrete(limits = rev(ord)) +   # show subtype letters A-D (matches Panel D forest)
  coord_cartesian(xlim = c(0.5, 1.95)) +
  labs(x = "COVID pre-taste OR (one-vs-rest)", y = NULL) + thmD

## severity VERTICAL grouped bar -----------------------------------------------
# Drawn as two blocks on separate y scales. On a shared 0-100 axis the three
# chronicity/co-occurrence metrics (62-82%) fill the panel while the two
# functional-impact metrics (9-23%) compress into the bottom sixth, which hides the
# only contrast that discriminates the subtypes. Wilson intervals have capped ends
# so the B-vs-D overlap on smell impact is judgeable rather than implied. Each
# block carries its own y-axis range in the title, because bar heights are
# comparable within a block but NOT between blocks.
sev[, block := ifelse(grepl("impact", metric), "impact", "flat")]

sev_bar <- function(dat, ymax, ybreaks, ylab)
  ggplot(dat, aes(metric, pct, fill = sub, group = sub)) +
    geom_col(width = .78, colour = "grey25", linewidth = .2,
             position = position_dodge(width = .82)) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = .3, linewidth = .35,
                  colour = "grey30", position = position_dodge(width = .82)) +
    scale_fill_manual(values = SUB_PAL, guide = "none") +   # subtype key = shared bottom legend
    scale_y_continuous(limits = c(0, ymax), breaks = ybreaks,
                       expand = expansion(c(0, 0.02))) +
    labs(x = NULL, y = ylab) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          axis.title.y = element_text(size = 8.5, colour = "black"),
          axis.text.x = element_text(size = 8.5, colour = "black", lineheight = 0.9),
          axis.text.y = element_text(size = 8.5, colour = "black"),
          legend.position = "none", plot.margin = margin(4, 6, 4, 6))

p_sev <- (sev_bar(sev[block == "flat"],   100, seq(0, 100, 25), "% of cases (0-100)") |
          sev_bar(sev[block == "impact"],  30, seq(0,  30, 10), "% of cases (0-30)")) +
         plot_layout(widths = c(3, 2))

# =============================================================================
# ASSEMBLE — two files, so the master can place them independently:
#   D = HLA forest alone             -> narrow column, tucked under Panel C
#   E = COVID forest | severity bars -> full-width bottom row
# =============================================================================
titer <- function(t, sz = 16) plot_annotation(title = t,
  theme = theme(plot.title = element_text(face = "bold", size = sz, colour = "black")))

Dpanel <- p_hla + titer("D", 20)   # composited ~0.64x -> larger letter
ggsave(file.path(out, "Fig5_panelD_hla.png"), Dpanel,
       width = 150, height = 92, units = "mm", dpi = 300, bg = "white")   # wide enough for x-title at larger font
ggsave(file.path(out, "Fig5_panelD_hla.pdf"), Dpanel,
       width = 150, height = 92, units = "mm", device = cairo_pdf)

Epanel <- (p_cov | p_sev) + plot_layout(widths = c(1.05, 1.3)) + titer("E", 12)   # composited ~1.0x -> smaller letter
# height chosen so the full-width E row renders at the same height as Panel D's forest
# block in the composite (measured D_final 651 px vs W=2400 -> aspect h/w ~= 0.27).
ggsave(file.path(out, "Fig5_panelE_covidsev.png"), Epanel,
       width = 215, height = 62, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(out, "Fig5_panelE_covidsev.pdf"), Epanel,
       width = 215, height = 62, units = "mm", device = cairo_pdf)

cat("LEGEND NOTE: HLA forest = lead rs2071293 (chr6:32Mb, C4/TNXB/RCCX), effect",
    "allele A; ASSET subset {A,B,C} positive (OR ~1.13-1.18) vs {D} opposite/null (OR 0.96, n.s.),",
    "joint P=3.01e-7 (suggestive); APOE/CETP dropped (circular lipid loci). COVID = crude OR 1.60 -> age+sex-adj 0.90",
    "(null; age artefact). Severity bars are drawn in TWO BLOCKS ON SEPARATE Y SCALES",
    "(chronicity & co-occurrence 0-100%; functional impact 0-30%) because the two sets differ",
    "by an order of magnitude -- bar heights are comparable WITHIN a block, not between blocks.",
    "Chronicity and smell co-occurrence are flat across subtypes; functional impact is modestly",
    "higher in B, with B and D overlapping on smell impact. Error bars = Wilson 95% CI. Denominators:",
    "smell chronic/impact among smell-positive; others all cases. ASSET Manhattan in supplement.\n")
cat("DONE panelDE ->", out, "\n")
