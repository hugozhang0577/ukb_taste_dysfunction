#> in   output/ml_ready/group1_full.rds
#> in   output/subtyping/clusters/cluster_assignments_g1_{}_k4.rds
#> out  output/figures/fig3b_subtype_signature_forest.pdf
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
source("figures/_common/style.R")
source("figures/_common/fig3_type.R")   # local 10% type bump
suppressPackageStartupMessages(library(patchwork))

SUB <- "output/subtyping"
ML  <- "output/ml_ready"

# k-means cluster ids are arbitrary; the letter mapping has one source.
source("unsupervised_subtyping/_subtype_map.R")
SUB_TTL  <- c(A = "A  Aging frailty",     B = "B  Psychosomatic",
              C = "C  Cardiometabolic",   D = "D  Young idiopathic")

CORE    <- c("age", "BMI")                     # drawn as triangles
ANCHORS <- c("age", "BMI", "grip_max")         # pinned to the top rows, in this order
PANEL_VARS <- list(
  A = c("age", "BMI", "waist_circ", "grip_max", "hdl", "shbg", "urate", "GlycA"),
  B = c("age", "BMI", "phq9_total", "gad7_total", "neuroticism",
        "phq15_somatic_total", "pain_interference_mean", "dx_count_psych"),
  C = c("age", "BMI", "GlycA", "hba1c", "triglycerides", "hdl", "urate", "cystatin_c"),
  D = c("age", "BMI", "cog_symbol_digit", "grip_max", "igf1", "cystatin_c",
        "GlycA", "hba1c"))

LAB <- c(age = "Age", BMI = "BMI", waist_circ = "Waist circ.",
         grip_max = "Grip strength", hdl = "HDL-C", shbg = "SHBG",
         urate = "Urate", GlycA = "GlycA",
         phq9_total = "Depression", gad7_total = "Anxiety",
         neuroticism = "Neuroticism", phq15_somatic_total = "Somatic symptoms",
         pain_interference_mean = "Pain interference",
         dx_count_psych = "Psychiatric dx",
         hba1c = "HbA1c", triglycerides = "Triglycerides",
         cystatin_c = "Cystatin C", cog_symbol_digit = "Symbol-digit",
         igf1 = "IGF-1")

## ---- 1. within-sex z, pooled across sexes ----------------------------------
dat <- as.data.table(readRDS(file.path(ML, "group1_full.rds")))
cases <- dat[pheno_2w_strict == 1]
cat(sprintf("[load] ml_ready %d rows | cases %d\n", nrow(dat), nrow(cases)))
stopifnot(nrow(cases) == 5772)          # locked cohort size; fails loudly if inputs drift

need <- unique(unlist(PANEL_VARS))
miss <- setdiff(need, names(cases))
if (length(miss)) stop("absent from ml_ready: ", paste(miss, collapse = ", "))

zl <- list()
for (sx in c("m", "f")) {
  cl  <- readRDS(file.path(SUB, "clusters", sprintf("cluster_assignments_g1_%s_k4.rds", sx)))
  sub <- data.table(eid = as.integer(rownames(cl$Z_active)),
                    subtype = SUBTYPE_MAP[[sx]][as.character(cl$final_labels)])
  cs  <- cases[sex == (if (sx == "m") 1L else 0L)][sub, on = "eid", nomatch = 0]
  cat(sprintf("[%s] %d cases | %s\n", sx, nrow(cs),
              paste(names(table(cs$subtype)), table(cs$subtype), sep = "=", collapse = " ")))
  for (fid in need) {
    x <- suppressWarnings(as.numeric(cs[[fid]]))
    if (sum(!is.na(x)) < 50) { cat("   sparse, skipped:", sx, fid, "\n"); next }
    zl[[paste(sx, fid)]] <- data.table(
      feature_id = fid, subtype = cs$subtype,
      z = (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  }
}
Z <- rbindlist(zl)[!is.na(z) & !is.na(subtype)]
pooled <- Z[, .(m = mean(z), se = sd(z) / sqrt(.N), n = .N), by = .(feature_id, subtype)]
pooled[, `:=`(lo = m - 1.96 * se, hi = m + 1.96 * se)]

fin <- rbindlist(lapply(names(PANEL_VARS), function(s) {
  d <- pooled[subtype == s & feature_id %in% PANEL_VARS[[s]]]
  gap <- setdiff(PANEL_VARS[[s]], d$feature_id)
  if (length(gap)) stop("subtype ", s, " lost rows: ", paste(gap, collapse = ", "))
  d }))
fin[, `:=`(label = LAB[feature_id], core = feature_id %in% CORE)]
stopifnot(nrow(fin) == 32, !anyNA(fin$label))
stopifnot(max(nchar(fin$label)) <= 17)          # 42 mm cell budget in the 1x4 strip
cat(sprintf("[rows] %d (4 subtypes x 8) | z range %.2f to %.2f\n",
            nrow(fin), min(fin$lo), max(fin$hi)))

## ---- 2. layout -------------------------------------------------------------
W_MM <- 168; PITCH <- round(3.42 * FIG3_K, 2)   # pitch scales with the type size
XLIM <- c(-2.60, 2.60)
XBRK <- c(-1, 0, 1)
SZA  <- FIG3_SZ$anno / .pt

mk <- function(s) {
  d   <- fin[subtype == s]
  col <- SUBTYPE_COL[[s]]
  anch <- intersect(ANCHORS, d$feature_id)
  alab <- d[match(anch, feature_id), label]
  tlab <- d[!(feature_id %in% anch)][order(m), label]
  d[, label := factor(label, levels = c(tlab, rev(alab)))]   # age, BMI, grip on top
  ggplot(d, aes(m, label)) +
    geom_vline(xintercept = 0, linewidth = LW$ref, colour = "grey45") +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0,
                  linewidth = LW$err, colour = col) +
    geom_point(aes(shape = core), size = 1.3, colour = col) +
    geom_text(aes(label = sprintf("%+.2f", m)),
              hjust = ifelse(d$m >= 0, -0.30, 1.30),
              size = SZA * 0.92, family = FIG_FONT, colour = "grey35") +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17), guide = "none") +
    scale_x_continuous(limits = XLIM, breaks = XBRK, expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    ggtitle(SUB_TTL[[s]]) +
    theme_fig(grid = "none") +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_text(hjust = 1),
          plot.title.position = "plot",
          plot.title = element_text(size = FIG3_SZ$anno, face = "bold", colour = col,
                                    family = FIG_FONT, hjust = 0, margin = margin(b = 0.6)),
          plot.margin = margin(0.6, 1.5, 0.6, 1, "mm"))
}

plots <- lapply(c("A", "B", "C", "D"), mk)
# one shared axis title, drawn once instead of four times
xlab <- wrap_elements(full = grid::textGrob(
  "Standardized mean difference vs cohort (z)",
  gp = grid::gpar(fontfamily = FIG_FONT, fontsize = FIG3_SZ$axtit, col = "black")))

N_ROW  <- 8
# 11.6 mm = the fixed furniture of each cell (title, own x axis, margin), measured
H_CELL <- N_ROW * PITCH + 11.6 * FIG3_K
H_MM   <- round(H_CELL + 4.5 * FIG3_K, 1)        # one row of cells + the shared axis title
stopifnot(H_MM <= 207)
cat(sprintf("[layout] cell %.1f mm x 1 + axis title -> %.1f mm\n", H_CELL, H_MM))

fig <- ((plots[[1]] | plots[[2]] | plots[[3]] | plots[[4]]) / xlab) +
  plot_layout(heights = c(H_CELL, 4.5))
fig_save(fig, "fig3b_subtype_signature_forest", W_MM, H_MM)

cat("\n[rows plotted]\n")
for (s in c("A", "B", "C", "D")) {
  cat(sprintf("-- %s --\n", SUB_TTL[[s]]))
  print(fin[subtype == s][order(-abs(m)), .(label, z = round(m, 2), n)])
}
cat("\n[legend] Standardized mean difference of each variable against the whole-cohort\n")
cat("         mean, computed within sex and pooled; bars are 95% CI. Triangles mark the\n")
cat("         two variables shared by all four subtypes (age, BMI).\n")
cat("[legend] Descriptive: these variables entered the clustering, so no P values are\n")
cat("         shown. Independent characterisation is in panels c and d.\n")
