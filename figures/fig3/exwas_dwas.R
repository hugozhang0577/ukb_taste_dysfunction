# =============================================================================
# Fig 3 — Exposure-wide and disease-wide associations
#   Panel A  Baseline + online follow-up ExWAS volcano  - colour by source
#   Panel B  Signal by source (significant / tested bar)
#   Panel C  First-occurrence disease (DWAS) volcano    - colour by ICD chapter
#   Panel D  Signal by ICD chapter (significant / tested bar)
# x = log(OR)  (protective < 0 < risk),  y = -log10(P), nonlinear (expanded near 0).
# Significance = per-source FDR (A/B) or global FDR (C); colour only when significant.
# G2 direction-concordant hits ringed in black (incl. outlier triangles).
# Effects are per-SD standardised; the four panels are assembled with patchwork.
# Single source: output/evidence_tiering/tier_flagged_results.csv (same as Supp S9-S11).
# ASCII glyphs only; cairo_pdf safe.
# =============================================================================
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})
HAS_REPEL <- requireNamespace("ggrepel", quietly = TRUE)

# =============================================================================
#  >>>>>  TUNABLE PARAMETERS  —  editing this block is normally enough  <<<<<
#         The function bodies and the data reads below rarely need touching.
# =============================================================================

## --- 1. Label counts ---------------------------------------------------------
N_LABEL_P  <- 4    # per volcano, per side (protective / risk): label the N smallest-P points
N_LABEL_OR <- 2    # plus the N largest-|effect| points, whether or not they are significant.
                   #   These always show their OR, so an extreme effect on a non-significant
                   #   point cannot be mistaken for a significant one.

## --- 2. Figure size and panel proportions ------------------------------------
FIG_W_MM    <- 180          # overall width (mm)
FIG_H_MM    <- 215          # overall height (mm)
PANEL_RATIO <- c(2.25, 0.90) # width ratio per row = volcano : signal bar chart

## --- 3. Axis compression -----------------------------------------------------
##   One list per panel. Fields:
##     XLIM     : full x-axis (log OR) range c(left, right); only points outside are clipped.
##     X_BREAK  : piecewise x compression. c(bl, br) keeps [bl, br] linear, so the bulk of
##                the points do not move, and compresses the two outer intervals so that a
##                few distant points can be shown at their real value rather than clipped.
##                Set NULL to disable (out-of-range points are then clipped to the edge and
##                drawn as triangles). To compress one side only, set bl = XLIM[1] or
##                br = XLIM[2].
##     Y_BREAK  : y-axis (-log10 P) compression start; above it the axis compresses but the
##                true value is retained.
##     EFF_FLOOR: position of the grey effect reference lines, at +/-EFF_FLOOR.
PANEL_A <- list(XLIM = c(-0.7, 1.0), X_BREAK = NULL, Y_BREAK = 40, EFF_FLOOR = 0.2)
PANEL_C <- list(XLIM = c(-0.7, 1.0), X_BREAK = NULL, Y_BREAK = 12, EFF_FLOOR = 0.2)
# Example: to pull panel C's points below logOR -0.7 into view rather than clipping them:
#     PANEL_C <- list(XLIM = c(-1.15, 1.0), X_BREAK = c(-0.7, 1.0), Y_BREAK = 12, EFF_FLOOR = 0.2)
# =============================================================================

# The exposure-family and ICD-chapter labels are stored as display strings in
# the evidence-tiering table, so the panels label themselves with no lookup.

PAL_A <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD","#8C564B","#E377C2",
           "#7F7F7F","#BCBD22","#17BECF","#AEC7E8","#FFBB78","#98DF8A","#FF9896",
           "#C5B0D5","#C49C94","#F0B27A","#393B79","#8C6D31")
# ICD-10 chapter palette + phecode range mapper (fallback for unmapped phecodes)
PAL_ICD <- c(
  "A-B Infectious"="#D62728","C-D Neoplasms"="#7F3C8D","D50-89 Blood"="#8C564B",
  "E Endocrine"="#F28E2B","F Mental"="#B07AA1","G Nervous"="#59A14F","H Eye/Ear"="#1F77B4",
  "I Circulatory"="#E15759","J Respiratory"="#76B7B2","K Digestive"="#4E79A7",
  "L Skin"="#FF9DA7","M Musculoskel."="#9C755F","N Genitourinary"="#EDC948",
  "O Pregnancy"="#E377C2","P Perinatal"="#AEC7E8","Q Congenital"="#17BECF",
  "R Symptoms"="#BCBD22","S-T Injury"="#BAB0AC","V-Y External"="#7F7F7F",
  "Z Health status"="#C7C7C7","Other"="#B0B0B0","Unmapped"="#D0D0D0")
phecode_to_icd_chapter <- function(phe_var) {
  pc <- suppressWarnings(as.numeric(sub("^phe", "", phe_var)))
  fcase(
    is.na(pc),              "Unmapped",
    pc < 140,               "A-B Infectious",
    pc >= 140 & pc < 240,   "C-D Neoplasms",
    pc >= 280 & pc < 290,   "D50-89 Blood",
    pc >= 240 & pc < 280,   "E Endocrine",
    pc >= 290 & pc < 320,   "F Mental",
    pc >= 320 & pc < 380,   "G Nervous",
    pc >= 380 & pc < 390,   "H Eye/Ear",
    pc >= 390 & pc < 460,   "I Circulatory",
    pc >= 460 & pc < 520,   "J Respiratory",
    pc >= 520 & pc < 580,   "K Digestive",
    pc >= 680 & pc < 710,   "L Skin",
    pc >= 710 & pc < 740,   "M Musculoskel.",
    pc >= 580 & pc < 630,   "N Genitourinary",
    pc >= 780 & pc < 800,   "R Symptoms",
    pc >= 800,              "S-T Injury",
    default                 = "Other")
}
# accurate WHO ICD-10 map -> chapter (fig3 label style + O/P/Q/Z); falls back to range above
.ICDMAP <- NULL
phecode_icd_primary <- function(v) {
  if (is.null(.ICDMAP)) {
    m <- fread("input/reference/phecode_map_v1_2_icd10.csv")
    m[, icd3 := substr(ICD10, 1, 3)]
    .ICDMAP <<- m[, .(icd_primary = icd3[1]), by = PHECODE][, .(variable = paste0("phe", PHECODE), icd_primary)]
  }
  .ICDMAP[data.table(variable = v), on = "variable", x.icd_primary]
}
icd3_to_chapter <- function(icd3) {
  L <- substr(icd3, 1, 1); n <- suppressWarnings(as.integer(substr(icd3, 2, 3)))
  fcase(
    is.na(icd3) | icd3 == "",                    NA_character_,
    L %in% c("A","B"),                           "A-B Infectious",
    L == "C" | (L == "D" & !is.na(n) & n <= 48), "C-D Neoplasms",
    L == "D" & !is.na(n) & n >= 50,              "D50-89 Blood",
    L == "E","E Endocrine", L == "F","F Mental", L == "G","G Nervous", L == "H","H Eye/Ear",
    L == "I","I Circulatory", L == "J","J Respiratory", L == "K","K Digestive",
    L == "L","L Skin", L == "M","M Musculoskel.", L == "N","N Genitourinary",
    L == "O","O Pregnancy", L == "P","P Perinatal", L == "Q","Q Congenital",
    L == "R","R Symptoms", L %in% c("S","T"),"S-T Injury",
    L %in% c("V","W","X","Y"),"V-Y External", L == "Z","Z Health status",
    default = NA_character_)
}
orfmt <- function(x) ifelse(x>=100, sprintf("%.0f",x), ifelse(x>=10, sprintf("%.1f",x), sprintf("%.2f",x)))

# ---- load --------------------------------------------------------------------
tier <- fread("output/evidence_tiering/tier_flagged_results.csv")
ex <- tier[layer %in% c("ExWAS-baseline","ExWAS-followup","DWAS")]
ex[, logOR := ifelse(!is.na(OR) & OR > 0, log(OR), beta)]
ex[, nlp   := -log10(pmax(pval, 1e-300))]
ex[, tier_lab := fcase(pass_bonf == TRUE, "Bonferroni",
                       pass_fdr  == TRUE, "FDR",
                       default = "Non-sig")]
ex[, replicated := pass_g2_replication == TRUE]
ex[, concord    := g2_sign_concord == TRUE]
BONF_LINE <- -log10(0.05 / nrow(ex))   # global Bonferroni threshold (-log10P)
cat(sprintf("rows: questionnaire=%d  disease=%d | Bonf=%d FDR=%d G2-replicated=%d\n",
            nrow(ex[layer %in% c("ExWAS-baseline","ExWAS-followup")]),
            nrow(ex[layer == "DWAS"]),
            sum(ex$tier_lab=="Bonferroni"), sum(ex$tier_lab=="FDR"), sum(ex$replicated)))

# resolve pretty labels (tier file variable_label = raw name)
.lk <- rbindlist(list(
  fread("input/reference/exwas_baseline_labels.csv", select=c("variable","pretty_label")),
  fread("input/reference/exwas_followup_labels.csv", select=c("variable","pretty_label"))), use.names=TRUE)
.cc <- fread("output/dwas/results/dwas_fdr_significant_table.csv")
.cc[, variable := paste0("phe", PheCode)]
.lk <- unique(rbind(.lk, .cc[, .(variable, pretty_label = Description)]), by="variable")
ex <- merge(ex, .lk, by="variable", all.x=TRUE, sort=FALSE)
ex[is.na(pretty_label) | pretty_label=="", pretty_label := variable]

# per-SD standardisation for continuous exposures (comparable; removes raw-unit
# artefacts e.g. WHR/HLS-retic). Binary/ordinal/disease kept per-level/per-category.
.sdlk <- fread("output/evidence_tiering/exwas_continuous_SD.csv")
ex <- merge(ex, .sdlk, by = "variable", all.x = TRUE, sort = FALSE)
ex[, logOR_raw := logOR]
ex[var_type == "continuous" & !is.na(sd), logOR := logOR_raw * sd]
ex[, nlp := -log10(pmax(pval, 1e-300))]   # unchanged by standardisation
cat(sprintf("per-SD standardised continuous vars: %d\n", ex[var_type=="continuous" & !is.na(sd), .N]))

# piece-wise x-axis transform: linear inside [bl,br], compressed in the two outer
# bands so a few far points show their TRUE value without shifting the bulk.
.xtrans <- function(lo, hi, bl, br, fL = 0.12, fR = 0.12) {
  if (bl <= lo) fL <- 0                            # no left compression band
  if (br >= hi) fR <- 0                            # no right compression band
  wM <- 1 - fL - fR
  dL <- (bl - lo) + 1e-9; dR <- (hi - br) + 1e-9; dM <- (br - bl)
  tf <- function(z) { z <- pmax(pmin(z, hi), lo)
    ifelse(z < bl, (z - lo)/dL * fL,
      ifelse(z > br, (1 - fR) + (z - br)/dR * fR, fL + (z - bl)/dM * wM)) }
  iv <- function(u)
    ifelse(u < fL, lo + u/(fL + 1e-9) * dL,
      ifelse(u > 1 - fR, br + (u - (1 - fR))/(fR + 1e-9) * dR, bl + (u - fL)/wM * dM))
  scales::trans_new("xbrk", transform = tf, inverse = iv)
}

# ---- polished volcano builder copied from A/C prototypes ---------------------
A_CORE_LABELS <- c(
  "phq15_somatic_total", "phq9_total", "pain_health_today", "pain_month_count",
  "health_rating", "phq15_nausea", "gad7_total", "n_medications", "phq9_binary",
  "hdl", "crp", "hba1c", "getting_up", "cog_fluid_iq", "cog_symbol_digit",
  "income", "vitamin_d"
)
A_SHORT_LABEL <- c(
  phq15_somatic_total = "PHQ-15 somatic symptoms",
  phq15_fatigue = "PHQ-15 fatigue",
  phq9_total = "PHQ-9 depression score",
  pain_health_today = "Self-rated health today",
  pain_month_count = "Monthly pain-site count",
  health_rating = "Self-rated health",
  phq15_nausea = "PHQ-15 nausea",
  gad7_total = "GAD-7 anxiety score",
  n_medications = "Number of medications",
  phq9_binary = "Clinical depression",
  hdl = "HDL cholesterol",
  crp = "C-reactive protein",
  hba1c = "HbA1c",
  getting_up = "Ease of getting up",
  cog_fluid_iq = "Fluid intelligence",
  cog_symbol_digit = "Symbol-digit substitution",
  income = "Household income",
  vitamin_d = "Vitamin D"
)
A_POINT_NUDGE <- data.table(
  variable = c("n_medications", "phq2_score"),
  dx = c(-0.012, 0.012),
  dy = c(1.85, -1.45)
)

C_CORE_LABELS <- c(
  "phe332", "phe471", "phe369", "phe369.5", "phe79", "phe296", "phe296.2",
  "phe495", "phe530", "phe530.1", "phe530.11", "phe465.2", "phe53",
  "phe411.3", "phe475"
)
C_SHORT_LABEL <- c(
  "phe332" = "Parkinson's disease",
  "phe471" = "Nasal polyps",
  "phe369" = "Infection of the eye",
  "phe369.5" = "Infectious conjunctivitis",
  "phe79" = "Viral infection",
  "phe296" = "Mood disorders",
  "phe296.2" = "Depression",
  "phe495" = "Asthma",
  "phe530" = "Diseases of esophagus",
  "phe530.1" = "Esophagitis / GERD",
  "phe530.11" = "GERD",
  "phe465.2" = "Acute pharyngitis",
  "phe53" = "Herpes zoster",
  "phe411.3" = "Angina pectoris",
  "phe475" = "Chronic sinusitis"
)
C_POINT_NUDGE <- data.table(
  variable = c("phe369", "phe369.5", "phe296", "phe296.2", "phe530.11"),
  dx = c(-0.010, 0.010, 0.010, 0.014, 0.012),
  dy = c(0.60, -0.50, 0.55, -0.45, 0.35)
)

symmetric_x_focus_trans <- function(limit = 1.0, focus = 0.35, focus_width = 0.64) {
  tail_width <- (1 - focus_width) / 2
  trans <- function(z) {
    z <- pmax(pmin(z, limit), -limit)
    ifelse(z < -focus,
           (z + limit) / (limit - focus) * tail_width,
           ifelse(z <= focus,
                  tail_width + (z + focus) / (2 * focus) * focus_width,
                  tail_width + focus_width + (z - focus) / (limit - focus) * tail_width))
  }
  inv <- function(u) {
    ifelse(u < tail_width,
           -limit + u / tail_width * (limit - focus),
           ifelse(u <= tail_width + focus_width,
                  -focus + (u - tail_width) / focus_width * (2 * focus),
                  focus + (u - tail_width - focus_width) / tail_width * (limit - focus)))
  }
  scales::trans_new("symmetric_x_focus", transform = trans, inverse = inv)
}

make_volcano_polished <- function(d, grp_pal, title, cfg, core_labels, short_label,
                                  point_nudge, y_dense, y_focus, y_break,
                                  y_weights, jitter_cfg, point_cfg,
                                  threshold_band = NULL,
                                  top_expand = 0.04, repel_force = 3.0,
                                  repel_box_pad = 0.45) {
  xl <- c(-1.0, 1.0)
  eff_floor <- cfg$EFF_FLOOR
  d <- copy(d)
  d[, x := pmax(pmin(logOR, xl[2]), xl[1])]
  d[, outlier := logOR > xl[2] | logOR < xl[1]]
  d[, global_bonf := pval < (0.05 / nrow(ex))]
  d[, sig_level := fcase(global_bonf == TRUE, "Global Bonferroni",
                         pass_fdr == TRUE, "FDR",
                         default = "Non-sig")]
  d[, or_val := exp(logOR)]

  ns <- d[sig_level == "Non-sig"]
  sig_fdr <- d[sig_level == "FDR"]
  sig_bonf <- d[sig_level == "Global Bonferroni"]

  mx <- max(d$nlp, na.rm = TRUE)
  w0 <- y_weights[1]; w1 <- y_weights[2]; w2 <- y_weights[3]; w3 <- y_weights[4]
  visual_y <- function(z) {
    z <- pmax(z, 0)
    y_max <- max(mx, y_break + 1)
    ifelse(z <= y_dense, (z / y_dense) * w0,
           ifelse(z <= y_focus, w0 + (z - y_dense) / (y_focus - y_dense) * w1,
                  ifelse(z <= y_break, w0 + w1 + (z - y_focus) / (y_break - y_focus) * w2,
                         w0 + w1 + w2 + (z - y_break) / (y_max - y_break) * w3)))
  }
  inverse_visual_y <- function(t) {
    t <- pmax(t, 0)
    y_max <- max(mx, y_break + 1)
    ifelse(t <= w0, t / w0 * y_dense,
           ifelse(t <= w0 + w1, y_dense + (t - w0) / w1 * (y_focus - y_dense),
                  ifelse(t <= w0 + w1 + w2, y_focus + (t - w0 - w1) / w2 * (y_break - y_focus),
                         y_break + (t - w0 - w1 - w2) / w3 * (y_max - y_break))))
  }

  density_jitter <- function(dt, x_width, y_width, x_bin, y_bin, seed) {
    dt <- copy(dt)
    if (!nrow(dt)) return(dt[, `:=`(xj = numeric(), yj = numeric())])
    dt[, `:=`(x_bin_local = floor((x - xl[1]) / x_bin),
              y_bin_local = floor(nlp / y_bin))]
    dt[, local_n := .N, by = .(x_bin_local, y_bin_local)]
    max_n <- max(dt$local_n, na.rm = TRUE)
    dt[, density_scale := if (max_n > 1) sqrt(local_n) / sqrt(max_n) else 0]
    dt[, `:=`(jw_local = x_width * (0.25 + 0.75 * density_scale),
              jh_local = y_width * (0.20 + 0.80 * density_scale))]
    set.seed(seed)
    dt[, `:=`(xj = x + runif(.N, -jw_local, jw_local),
              yj = pmax(nlp + runif(.N, -jh_local, jh_local), 0))]
    dt[, c("x_bin_local", "y_bin_local", "local_n",
           "density_scale", "jw_local", "jh_local") := NULL]
    dt
  }
  apply_point_nudge <- function(dt) {
    if (!nrow(dt)) return(dt)
    dt <- merge(dt, point_nudge, by = "variable", all.x = TRUE, sort = FALSE)
    dt[is.na(dx), dx := 0]
    dt[is.na(dy), dy := 0]
    dt[, `:=`(xj = xj + dx, yj = pmax(yj + dy, 0))]
    dt[, c("dx", "dy") := NULL]
    dt
  }
  enforce_bonf_boundary <- function(dt, passed_global_bonf, pad = 0.035) {
    dt <- copy(dt)
    if (!nrow(dt)) return(dt)
    if (passed_global_bonf) {
      dt[, yj := pmax(yj, BONF_LINE + pad)]
    } else {
      dt[, yj := pmin(yj, BONF_LINE - pad)]
      dt[, yj := pmax(yj, 0)]
    }
    dt
  }
  separate_overlaps_y <- function(dt, x_bin, y_bin, step, max_shift) {
    dt <- copy(dt)
    if (!nrow(dt)) return(dt)
    dt[, y_visual := visual_y(yj)]
    dt[, `:=`(stack_xbin = floor((xj - xl[1]) / x_bin),
              stack_ybin = floor(y_visual / y_bin))]
    dt[, stack_n := .N, by = .(stack_xbin, stack_ybin)]
    dt[, stack_rank := seq_len(.N), by = .(stack_xbin, stack_ybin)]
    dt[stack_n > 1, y_stack_shift := (stack_rank - (stack_n + 1) / 2) * step]
    dt[stack_n <= 1, y_stack_shift := 0]
    dt[, y_stack_shift := pmax(pmin(y_stack_shift, max_shift), -max_shift)]
    dt[, yj := pmax(inverse_visual_y(y_visual + y_stack_shift), 0)]
    dt[, c("stack_xbin", "stack_ybin", "stack_n",
           "stack_rank", "y_stack_shift", "y_visual") := NULL]
    dt
  }
  separate_threshold_band <- function(fdr_dt, bonf_dt, band, x_bin = 0.075,
                                      step_visual = 0.070, max_shift_visual = 0.32,
                                      pad = 0.045) {
    lane_one_side <- function(dt, side = c("below", "above")) {
      side <- match.arg(side)
      dt <- copy(dt)
      if (!nrow(dt)) return(dt)
      idx <- which(dt$yj >= band[1] & dt$yj <= band[2])
      if (!length(idx)) return(dt)
      sub <- copy(dt[idx])
      sub[, y_visual := visual_y(yj)]
      sub[, x_lane := floor((xj - xl[1]) / x_bin)]
      if (side == "above") {
        setorder(sub, x_lane, y_visual)
        sub[, lane_rank := seq_len(.N), by = x_lane]
        sub[, y_visual_new := pmax(
          y_visual, visual_y(BONF_LINE + pad) + (lane_rank - 1) * step_visual
        ), by = x_lane]
        sub[, y_visual_new := pmin(y_visual_new, y_visual + max_shift_visual)]
        sub[, yj := pmax(inverse_visual_y(y_visual_new), BONF_LINE + pad)]
      } else {
        setorder(sub, x_lane, -y_visual)
        sub[, lane_rank := seq_len(.N), by = x_lane]
        sub[, y_visual_new := pmin(
          y_visual, visual_y(BONF_LINE - pad) - (lane_rank - 1) * step_visual
        ), by = x_lane]
        sub[, y_visual_new := pmax(y_visual_new, y_visual - max_shift_visual)]
        sub[, yj := pmin(inverse_visual_y(y_visual_new), BONF_LINE - pad)]
        sub[, yj := pmax(yj, 0)]
      }
      sub[, c("y_visual", "x_lane", "lane_rank", "y_visual_new") := NULL]
      dt[idx] <- sub
      dt
    }
    list(fdr = lane_one_side(fdr_dt, "below"),
         bonf = lane_one_side(bonf_dt, "above"))
  }

  ns <- density_jitter(ns, jitter_cfg$ns$x_width, jitter_cfg$ns$y_width,
                       jitter_cfg$ns$x_bin, jitter_cfg$ns$y_bin, jitter_cfg$ns$seed)
  sig_fdr <- density_jitter(sig_fdr, jitter_cfg$fdr$x_width, jitter_cfg$fdr$y_width,
                            jitter_cfg$fdr$x_bin, jitter_cfg$fdr$y_bin, jitter_cfg$fdr$seed)
  sig_bonf <- density_jitter(sig_bonf, jitter_cfg$bonf$x_width, jitter_cfg$bonf$y_width,
                             jitter_cfg$bonf$x_bin, jitter_cfg$bonf$y_bin, jitter_cfg$bonf$seed)
  sig_fdr <- apply_point_nudge(sig_fdr)
  sig_bonf <- apply_point_nudge(sig_bonf)
  sig_fdr <- separate_overlaps_y(sig_fdr, jitter_cfg$fdr$stack_x_bin,
                                 jitter_cfg$fdr$stack_y_bin, jitter_cfg$fdr$step,
                                 jitter_cfg$fdr$max_shift)
  sig_bonf <- separate_overlaps_y(sig_bonf, jitter_cfg$bonf$stack_x_bin,
                                  jitter_cfg$bonf$stack_y_bin, jitter_cfg$bonf$step,
                                  jitter_cfg$bonf$max_shift)
  sig_fdr <- enforce_bonf_boundary(sig_fdr, FALSE)
  sig_bonf <- enforce_bonf_boundary(sig_bonf, TRUE)
  ns <- enforce_bonf_boundary(ns, FALSE)
  if (!is.null(threshold_band)) {
    sep <- separate_threshold_band(sig_fdr, sig_bonf, threshold_band)
    sig_fdr <- sep$fdr
    sig_bonf <- sep$bonf
  }

  labs <- d[variable %in% core_labels]
  labs[, lab_txt := fifelse(variable %in% names(short_label),
                            short_label[variable], pretty_label)]
  labs[outlier == TRUE, lab_txt := sprintf("%s (OR=%s)", lab_txt, orfmt(or_val))]
  labs[, side := ifelse(logOR < 0, -1, 1)]
  plotted_xy <- rbind(
    sig_fdr[, .(variable, x_lab = xj, y_lab = yj)],
    sig_bonf[, .(variable, x_lab = xj, y_lab = yj)],
    ns[, .(variable, x_lab = xj, y_lab = yj)],
    use.names = TRUE, fill = TRUE
  )
  plotted_xy <- unique(plotted_xy, by = "variable")
  labs <- merge(labs, plotted_xy, by = "variable", all.x = TRUE, sort = FALSE)
  labs[is.na(x_lab), x_lab := x]
  labs[is.na(y_lab), y_lab := nlp]

  ytr <- scales::trans_new("y_piecewise", transform = visual_y, inverse = inverse_visual_y)
  y_breaks <- sort(unique(c(0, 2, 5, 10, 20, 40, 80, round(mx))))
  y_breaks <- y_breaks[y_breaks <= mx * 1.02]

  p <- ggplot() +
    geom_vline(xintercept = c(-eff_floor, eff_floor), linetype = 2,
               colour = "grey78", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.35) +
    geom_hline(yintercept = BONF_LINE, linetype = 2,
               colour = "#C0392B", linewidth = 0.35) +
    annotate("text", x = xl[1], y = BONF_LINE,
             label = "global Bonferroni reference",
             hjust = 0, vjust = -0.4, size = 2.2, colour = "#C0392B") +
    geom_point(data = ns, aes(xj, yj), colour = "grey70",
               size = point_cfg$ns_size, alpha = 0.62) +
    geom_point(data = sig_fdr[outlier == FALSE], aes(xj, yj, fill = grp),
               shape = 21, colour = "grey55", stroke = point_cfg$fdr_stroke,
               size = point_cfg$fdr_size, alpha = 0.93) +
    geom_point(data = sig_bonf[outlier == FALSE], aes(xj, yj, fill = grp),
               shape = 21, colour = "grey55", stroke = point_cfg$bonf_stroke,
               size = point_cfg$bonf_size, alpha = 0.93) +
    geom_point(data = sig_fdr[outlier == FALSE & concord == TRUE], aes(xj, yj),
               shape = 1, colour = "black", size = point_cfg$fdr_size + 0.10,
               stroke = point_cfg$fdr_ring) +
    geom_point(data = sig_bonf[outlier == FALSE & concord == TRUE], aes(xj, yj),
               shape = 1, colour = "black", size = point_cfg$bonf_size + 0.10,
               stroke = point_cfg$bonf_ring) +
    geom_point(data = sig_fdr[outlier == TRUE], aes(x, nlp, fill = grp),
               shape = 24, colour = "grey55", size = point_cfg$fdr_size + 0.55,
               stroke = point_cfg$fdr_stroke + 0.06, alpha = 0.65) +
    geom_point(data = sig_bonf[outlier == TRUE], aes(x, nlp, fill = grp),
               shape = 24, colour = "grey55", size = point_cfg$bonf_size + 0.50,
               stroke = point_cfg$bonf_stroke + 0.12) +
    geom_point(data = sig_fdr[outlier == TRUE & concord == TRUE], aes(x, nlp),
               shape = 24, colour = "black", fill = NA,
               size = point_cfg$fdr_size + 0.63, stroke = 0.30) +
    geom_point(data = sig_bonf[outlier == TRUE & concord == TRUE], aes(x, nlp),
               shape = 24, colour = "black", fill = NA,
               size = point_cfg$bonf_size + 0.58, stroke = 0.42) +
    scale_fill_manual(values = grp_pal, guide = "none") +
    scale_x_continuous(
      "Log odds ratio (protective < 0 < risk)",
      trans = symmetric_x_focus_trans(),
      breaks = c(-1.0, -0.7, -0.35, -0.2, 0, 0.2, 0.35, 0.7, 1.0),
      labels = function(x) sprintf("%.2g", x)
    ) +
    scale_y_continuous(
      bquote(-log[10]~italic(P)~" (nonlinear scale; expanded below "~.(y_dense)~")"),
      trans = ytr,
      breaks = y_breaks,
      expand = expansion(c(0, top_expand))
    ) +
    coord_cartesian(xlim = xl, clip = "off") +
    labs(title = title) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 10.8),
      axis.title = element_text(size = 8.5),
      axis.text = element_text(size = 7.2),
      legend.position = "none",
      plot.margin = margin(5, 14, 5, 6)
    )

  if (nrow(labs) && HAS_REPEL) {
    p <- p + ggrepel::geom_text_repel(
      data = labs,
      aes(x_lab, y_lab, label = lab_txt),
      size = 2.25, colour = "black", fontface = "plain",
      max.overlaps = Inf, min.segment.length = 0,
      point.padding = 0.22, segment.size = 0.16,
      segment.colour = "grey65", box.padding = repel_box_pad,
      force = repel_force, force_pull = 0.12,
      max.iter = 40000, max.time = 3,
      nudge_x = labs$side * 0.13,
      hjust = ifelse(labs$side < 0, 1, 0),
      seed = 12
    )
  } else if (nrow(labs)) {
    p <- p + geom_text(data = labs, aes(x_lab, y_lab, label = lab_txt), size = 2)
  }
  p
}

# ---- signal proportion bar: significant / tested by class --------------------
make_sigbar <- function(d, pal, title, axis_lab = "Significant / tested (%)") {
  d <- copy(d)
  d[, sig := tier_lab != "Non-sig"]
  agg <- d[, .(n_sig = sum(sig), n_tot = .N), by = grp]
  agg[, prop_sig := fifelse(n_tot > 0, n_sig / n_tot, 0)]
  agg[, lab := sprintf("%d/%d", n_sig, n_tot)]
  agg <- agg[order(n_sig, prop_sig)]
  lv <- as.character(agg$grp)
  agg[, grp := factor(as.character(grp), levels = lv)]
  tick_x <- -0.028

  ggplot(agg, aes(y = grp)) +
    geom_col(aes(x = 1), width = .66, fill = "grey89",
             colour = "white", linewidth = .18) +
    geom_col(aes(x = prop_sig, fill = grp), width = .66,
             colour = "white", linewidth = .18) +
    geom_point(aes(x = tick_x, colour = grp),
               size = 2.2, shape = 16) +
    geom_text(aes(x = 0.985, label = lab), hjust = 1,
              size = 1.75, colour = "grey25") +
    scale_fill_manual(values = pal, guide = "none") +
    scale_colour_manual(values = pal, guide = "none") +
    scale_x_continuous(
      axis_lab,
      limits = c(-0.055, 1.02),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0", "25", "50", "75", "100"),
      expand = expansion(c(0, 0))
    ) +
    labs(title = title) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "grey90", linewidth = .25),
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 7),
      axis.text.y = element_text(size = 6.6, colour = "grey30"),
      axis.text.x = element_text(size = 6.6),
      plot.title = element_text(face = "bold", size = 10.8),
      plot.margin = margin(3, 2, 3, 4),
      axis.ticks.x = element_line(colour = "grey70")
    )
}

# ---- panels A and B : baseline + follow-up questionnaire scans ----------------
# Selected on `layer`, the same column the load step filtered on, so the two
# panels cannot end up selecting on a vocabulary the tiering step has moved on
# from.
da <- ex[layer %in% c("ExWAS-baseline","ExWAS-followup")]
da[, grp := domain]
da[is.na(grp), grp := "Other (handedness)"]
grp_a <- sort(unique(da$grp)); pal_a <- setNames(PAL_A[seq_along(grp_a)], grp_a)
pal_a["General & Mental Health"] <- "#F28E2B"
pal_a["Sensory & Pain"] <- "#4C78A8"
vol_q <- make_volcano_polished(
  da, pal_a, "A", cfg = PANEL_A,
  core_labels = A_CORE_LABELS, short_label = A_SHORT_LABEL,
  point_nudge = A_POINT_NUDGE, y_dense = 10, y_focus = 20, y_break = 40,
  y_weights = c(0.46, 0.18, 0.24, 0.12),
  jitter_cfg = list(
    ns = list(x_width = 0.012, y_width = 0.45, x_bin = 0.050, y_bin = 2, seed = 13),
    fdr = list(x_width = 0.014, y_width = 0.80, x_bin = 0.045, y_bin = 4, seed = 11,
               stack_x_bin = 0.030, stack_y_bin = 0.032, step = 0.018, max_shift = 0.075),
    bonf = list(x_width = 0.018, y_width = 0.95, x_bin = 0.050, y_bin = 5, seed = 12,
                stack_x_bin = 0.034, stack_y_bin = 0.036, step = 0.024, max_shift = 0.095)
  ),
  point_cfg = list(ns_size = 0.32, fdr_size = 1.20, bonf_size = 1.90,
                   fdr_stroke = 0.16, bonf_stroke = 0.16,
                   fdr_ring = 0.26, bonf_ring = 0.40),
  top_expand = 0.16, repel_force = 7.0, repel_box_pad = 0.66
)
bar_q <- make_sigbar(da, pal_a, "B")

# ---- panels C and D : disease-wide scan --------------------------------------
defs <- fread("input/reference/phecode_definitions_v1_2.csv",
              select = c("phecode","phenotype"))
defs[, variable := paste0("phe", phecode)]
db <- merge(ex[layer == "DWAS"], defs[, .(variable, phenotype)],
            by = "variable", all.x = TRUE)
db[!is.na(phenotype) & phenotype != "", pretty_label := phenotype]
db[, grp := icd3_to_chapter(phecode_icd_primary(variable))]   # accurate WHO ICD-10 chapter
db[is.na(grp), grp := phecode_to_icd_chapter(variable)]        # fig3 range fallback for unmapped
db[is.na(grp) | grp == "Unmapped" | !grp %in% names(PAL_ICD), grp := "Other"]
present <- names(PAL_ICD)[names(PAL_ICD) %in% unique(db$grp)]
db[, grp := factor(grp, levels = present)]
cat("panel b ICD-10 chapters:
"); print(db[, .N, by = grp][order(-N)])
vol_d <- make_volcano_polished(
  db, PAL_ICD, "C", cfg = PANEL_C,
  core_labels = C_CORE_LABELS, short_label = C_SHORT_LABEL,
  point_nudge = C_POINT_NUDGE, y_dense = 5, y_focus = 10, y_break = 20,
  y_weights = c(0.54, 0.22, 0.15, 0.09),
  jitter_cfg = list(
    ns = list(x_width = 0.016, y_width = 0.50, x_bin = 0.050, y_bin = 2, seed = 31),
    fdr = list(x_width = 0.018, y_width = 0.80, x_bin = 0.045, y_bin = 3, seed = 32,
               stack_x_bin = 0.045, stack_y_bin = 0.055, step = 0.040, max_shift = 0.160),
    bonf = list(x_width = 0.022, y_width = 1.05, x_bin = 0.050, y_bin = 4, seed = 33,
                stack_x_bin = 0.050, stack_y_bin = 0.060, step = 0.050, max_shift = 0.200)
  ),
  point_cfg = list(ns_size = 0.42, fdr_size = 1.55, bonf_size = 2.35,
                   fdr_stroke = 0.18, bonf_stroke = 0.18,
                   fdr_ring = 0.30, bonf_ring = 0.44),
  threshold_band = c(3.8, 7.5),
  top_expand = 0.06, repel_force = 4.0, repel_box_pad = 0.50
)
bar_d <- make_sigbar(db, PAL_ICD, "D")

# ---- assemble ----------------------------------------------------------------
fig <- vol_q + bar_q + vol_d + bar_d +
  plot_layout(design = "AB
CD", widths = PANEL_RATIO)

outdir <- "output/figures"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(outdir, "Fig3_exwas_dwas.pdf"), fig,
       width = FIG_W_MM, height = FIG_H_MM, units = "mm", device = cairo_pdf)
ggsave(file.path(outdir, "Fig3_exwas_dwas.png"), fig,
       width = FIG_W_MM, height = FIG_H_MM, units = "mm", dpi = 400, bg = "white")
cat("Saved Fig3_exwas_dwas.{pdf,png} to", outdir, "\n")
