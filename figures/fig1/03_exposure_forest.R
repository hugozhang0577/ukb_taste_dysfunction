#> in   input/reference/exwas_baseline_labels.csv
#> in   input/reference/exwas_followup_labels.csv
#> in   input/reference/phecode_definitions_v1_2.csv
#> in   output/evidence_tiering/exwas_continuous_SD.csv
#> in   output/evidence_tiering/tier_flagged_results.csv
#> out  output/figures/fig1_exposure_forest.pdf
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
source("figures/_common/style.R")
suppressPackageStartupMessages(library(patchwork))

## ---- 1. load and per-SD standardisation ------------------------------------
t  <- fread("output/evidence_tiering/tier_flagged_results.csv")
ex <- t[layer %in% c("ExWAS-A", "ExWAS-B", "ExWAS-C")]
cat(sprintf("[load] ExWAS rows: A=%d B=%d C=%d\n",
            ex[layer == "ExWAS-A", .N], ex[layer == "ExWAS-B", .N], ex[layer == "ExWAS-C", .N]))

sdlk <- fread("output/evidence_tiering/exwas_continuous_SD.csv")
ex <- merge(ex, sdlk, by = "variable", all.x = TRUE, sort = FALSE)
ex[, `:=`(lo_ = log(OR), lolo_ = log(CI_lo), lohi_ = log(CI_hi))]
ex[var_type == "continuous" & !is.na(sd),
   `:=`(lo_ = lo_ * sd, lolo_ = lolo_ * sd, lohi_ = lohi_ * sd)]
ex[, `:=`(ORd = exp(lo_), LOd = exp(lolo_), HId = exp(lohi_))]
cat(sprintf("[per-SD] %d continuous variables standardised; FDR-sig OR range %.2f-%.2f (raw max %.0f)\n",
            ex[var_type == "continuous" & !is.na(sd), .N],
            min(ex[pass_fdr == TRUE, ORd]), max(ex[pass_fdr == TRUE, ORd]),
            max(ex[pass_fdr == TRUE, OR])))

lk <- rbindlist(list(
  fread("input/reference/exwas_baseline_labels.csv",
        select = c("variable", "pretty_label")),
  fread("input/reference/exwas_followup_labels.csv",
        select = c("variable", "pretty_label"))), use.names = TRUE)
defs <- fread("input/reference/phecode_definitions_v1_2.csv",
              select = c("phecode", "phenotype"))
defs[, variable := paste0("phe", phecode)]
lk <- unique(rbind(lk, defs[, .(variable, pretty_label = phenotype)]), by = "variable")
ex <- merge(ex, lk, by = "variable", all.x = TRUE, sort = FALSE)
ex[is.na(pretty_label) | pretty_label == "", pretty_label := variable]

## ---- 2. labels: strip prefixes, apply short names, then truncate -----------
STRIP <- c("^Pain \\(month\\): ", "^Somatic symptoms: ", "^Employment: ", "^Qualification: ",
           "^Doctor-diagnosed ", "^Preference: ", "^Milk: ", "^Never eats: ",
           "^Workplace ", "^Work ", "^Cognitive: ", "^Disorders of ", "^Frequency of ")
strip_prefix <- function(x) { for (p in STRIP) x <- sub(p, "", x); x }
sentence1 <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))

SHORT <- c(
  # --- oral (fixed row set) ---
  oral_painful_gums_baseline = "Painful gums", oral_mouth_ulcers_baseline = "Mouth ulcers",
  oral_bleeding_gums_baseline = "Bleeding gums", periodontal_indicator_baseline = "Periodontal composite",
  oral_toothache_baseline = "Toothache", any_oral_problem_baseline = "Any oral problem",
  oral_loose_teeth_baseline = "Loose teeth", oral_dentures_baseline = "Denture wearing",
  taste_drug_user_baseline_numeric = "Taste-affecting medication",
  surg_taste_affecting_baseline = "Taste-affecting surgery",
  surg_non_taste_affecting_baseline = "Surgery, other sites",
  surg_tooth_extraction_baseline = "Tooth extraction",
  "phe521.1" = "Dental caries, recorded", phe522 = "Pulp/periapical, recorded",
  phe527 = "Salivary gland, recorded", phe529 = "Tongue disorder, recorded",
  # --- pain and nerve ---
  pain_month_pain_all_over_the_body = "Pain all over the body",
  pain_month_stomach_or_abdominal_pain = "Stomach/abdominal pain",
  pain_month_neck_or_shoulder_pain = "Neck or shoulder pain",
  pain_month_count = "Monthly pain-site count", pain_site_count = "Pain-site count",
  pain_chronic = "Chronic pain", pain_widespread = "Widespread pain",
  pain_interference_mean = "Pain interference", pain_dn4_positive = "DN4-positive neuropathic pain",
  pain_dn4_score = "DN4 neuropathic pain score", cond_neuropathy = "Peripheral neuropathy",
  cond_fibromyalgia = "Fibromyalgia", cond_migraine = "Migraine (self-reported)",
  hearing_diff = "Hearing difficulty",
  # --- somatic and digestive symptoms ---
  phq15_somatic_total = "Somatic symptoms (4-item)", phq15_fatigue = "Fatigue",
  phq15_dizziness = "Dizziness", phq15_nausea = "Nausea", phq15_headache = "Headache",
  ibs_ever = "Irritable bowel syndrome", abdom_distension = "Abdominal distension",
  sensitive_stomach = "Sensitive stomach", bowel_satisfaction = "Bowel satisfaction",
  abdom_discomfort_freq = "Abdominal discomfort",
  # --- mood and cognition ---
  phq9_binary = "Screen-positive depression", gad7_binary = "Screen-positive anxiety",
  phq9_total = "PHQ-9 depression score", gad7_total = "GAD-7 anxiety score",
  phq2_score = "PHQ-2 depression score", neuroticism = "Neuroticism",
  loneliness = "Loneliness", mood_swings = "Mood swings", depressed_mood = "Depressed mood",
  unenthusiasm = "Unenthusiasm", mania_irritable = "Mania: irritability",
  cog_symbol_digit = "Symbol-digit substitution", cog_tmt_b = "Trail-making B",
  cog_fluid_iq = "Fluid intelligence",
  # --- general health and medication ---
  health_rating = "Self-rated health", pain_health_today = "Self-rated health today",
  longstanding_ill = "Longstanding illness", n_medications = "Number of medications",
  # --- diet and supplements ---
  diet_change_ill = "Diet change due to illness", diet_change = "Major dietary change (5 y)",
  diet_variation = "Dietary variation", salt_added = "Salt added to food",
  supp_vitamin_b12 = "Vitamin B12 supplement", supp_calcium = "Calcium supplement",
  never_eat_wheat_products = "Never eats wheat products",
  milk_type_never_rarely_have_milk = "Never/rarely has milk",
  fpq_taste_salty = "Prefers salty foods", fpq_taste_fatty = "Prefers fatty foods",
  nutr_sodium_mg = "Dietary sodium", bread_type_white = "White bread",
  # --- work and sleep ---
  work_chemical_fumes = "Chemical fumes at work", work_dusty = "Dusty workplace",
  work_paints_thinners = "Paints or thinners at work", work_pesticides = "Pesticides at work",
  work_diesel = "Diesel exhaust at work", work_noisy = "Noisy workplace",
  work_chemical_score = "Workplace chemical score",
  insomnia = "Insomnia", getting_up = "Ease of getting up",
  daytime_dozing = "Daytime dozing", snoring = "Snoring", nap_daytime = "Daytime napping",
  # --- body composition and biochemistry ---
  body_fat_pct = "Body fat percentage", fat_free_mass = "Fat-free mass",
  waist_circ = "Waist circumference", whole_body_fat = "Whole-body fat mass",
  grip_right = "Grip strength", grip_max = "Grip strength", dbp = "Diastolic blood pressure",
  hdl = "HDL cholesterol", triglycerides = "Triglycerides", urate = "Urate",
  cystatin_c = "Cystatin C", hba1c = "HbA1c", crp = "C-reactive protein",
  apoa = "Apolipoprotein A1", hls_retic = "High light scatter retics",
  vitamin_d = "Vitamin D", igf1 = "IGF-1", shbg = "SHBG",
  # --- self-reported history ---
  vascular_angina = "Angina", vascular_heart_attack = "Heart attack",
  vascular_stroke = "Stroke", vascular_high_blood_pressure = "High blood pressure",
  diabetes_sr = "Diabetes", vascular_count = "Vascular/heart conditions",
  fluid_intel = "Fluid intelligence",   # "(baseline)" dropped: 29 characters would truncate
  employment_unable_to_work_because_of_sickness_or_disability = "Unable to work: sickness",
  income = "Household income",
  # --- diagnoses ---
  phe332 = "Parkinson disease", phe471 = "Nasal polyps", phe475 = "Chronic sinusitis",
  phe495 = "Asthma", phe465 = "Upper respiratory infection", phe465.2 = "Acute pharyngitis",
  phe496.3 = "Bronchiectasis", phe496 = "Chronic airway obstruction",
  phe352 = "Cranial nerve disorders", phe340 = "Migraine (diagnosed)",
  phe369 = "Eye infection", phe369.5 = "Infectious conjunctivitis",
  phe380.1 = "Otitis externa", phe379 = "Other eye disorders",
  phe530.11 = "Gastroesophageal reflux", phe564 = "Functional GI disorder",
  phe564.1 = "Irritable bowel syndrome", phe411.3 = "Angina pectoris",
  phe427.2 = "Atrial fibrillation", phe433 = "Cerebrovascular disease",
  phe296 = "Mood disorders", phe296.2 = "Depression",
  phe276 = "Fluid/electrolyte disorders", phe458 = "Hypotension",
  # short names for labels that would truncate at 27 characters
  general_pain_3mo = "Pain all over body, 3 mo+",
  knee_pain_3mo  = "Knee pain, 3 months or more",
  neck_pain_3mo  = "Neck/shoulder pain, 3 mo+",
  abdom_pain_6mo = "Abdominal pain, past 6 mo",
  stressor_serious_illness_injury_or_assault_to_yourself = "Stressor: serious illness",
  # --- short names for the rows the Results names, which would truncate ---
  never_eat_dairy_products = "Never eats dairy",
  bread_type_wholemeal_or_wholegrain = "Bread: wholemeal",
  employment_in_paid_employment_or_self_employed = "In paid employment",
  qualification_college_or_university_degree = "University degree")

TRUNC <- 27
pretty_of <- function(v, lab) {
  out <- ifelse(v %in% names(SHORT), SHORT[v], sentence1(strip_prefix(lab)))
  ifelse(nchar(out) > TRUNC, paste0(substr(out, 1, TRUNC - 2), ".."), out)
}

## ---- 3. redundancy groups --------------------------------------------------------
PER_GRP_ALLOW <- c("Adiposity" = 2L, "Pain site" = 2L, "Pain burden" = 3L,
                   "Pain condition" = 2L, "Somatic symptom" = 3L,
                   "Occupational exposure" = 3L, "Sleep" = 2L, "Salt" = 1L,
                   "Supplement" = 2L, "Nutrient intake" = 1L)
redundancy_group <- function(v, lab) {
  fcase(
    grepl("^phe", v), paste0("phe_", sub("\\..*$", "", sub("^phe", "", v))),
    v %in% c("health_rating", "pain_health_today"), "Self-rated health",
    # The four salt variables (liking salty food, adding salt at follow-up and
    # at baseline, dietary sodium) are one construct and share a group.
    grepl("salt|^nutr_sodium", v),                  "Salt",
    grepl("^cond_", v),                             "Pain condition",
    grepl("^pain_month_[a-z]", v),                  "Pain site",
    v %in% c("pain_month_count", "pain_site_count"), "Pain count",
    grepl("^pain_(chronic|widespread|interference|avg|dn4)", v), "Pain burden",
    grepl("^phq15_", v),                            "Somatic symptom",
    v %in% c("phq9_binary", "gad7_binary"),         "Screen-positive threshold",
    grepl("^phq9|^phq2|^depressed_mood", v),        "Depression score",
    grepl("^gad7", v),                              "Anxiety score",
    grepl("^mania_", v),                            "Mania",
    grepl("^cog_|fluid_intel", v),                  "Cognition",
    grepl("^work_", v),                             "Occupational exposure",
    grepl("fat|weight|waist|hip_circ|^whr$|^bmi$", v), "Adiposity",
    grepl("grip", v),                               "Grip strength",
    grepl("retic", v),                              "Reticulocytes",
    grepl("^insomnia$|^getting_up$|dozing|^nap_|^snoring$|sleep_hours|chronotype", v), "Sleep",
    grepl("^diet_", v),                             "Dietary change",
    grepl("^nutr_", v),                             "Nutrient intake",
    grepl("^supp_", v),                             "Supplement",
    grepl("^qualification_|^employment_in_paid", v), "Education/employment",
    grepl("^vascular_(count|high_blood)", v),       "Vascular burden",
    default = v)
}
ex[, redgrp := redundancy_group(variable, pretty_label)]

## ---- 4. the eleven sub-panels ----------------------------------------------------
ORAL_FIXED <- names(SHORT)[1:16]   # 8 oral measures + 4 treatments or surgery + 4 hospital dental diagnoses

PANELS <- list(
  # ---- left column: dental and symptom panels ----
  list(title = "Oral health, dental treatment and recorded dental diagnoses",
       fixed = ORAL_FIXED, col = 1),
  list(title = "Pain and neuropathic features",                 dom = c("A7", "B5"),   k = 11,
       drop_grp = c("Self-rated health", "general_pain_3mo"),
       must = c("pain_month_pain_all_over_the_body", "pain_month_facial_pain",
                "pain_widespread", "cond_fibromyalgia", "hearing_diff"), col = 1),
  list(title = "Somatic and digestive symptoms",                dom = c("B3"),         k = 8,
       must = c("ibs_ever"), col = 1),
  list(title = "Mood, mental health and cognition",             dom = c("B4", "B7"),   k = 5,
       drop_grp = "audit_total", must = c("phq9_binary"), col = 1),
  # A1 has its own panel in the right column, so only A6 is left here.
  list(title = "General health, medication and circumstances",  dom = c("A6"),         k = 9,
       must = c("health_rating", "longstanding_ill", "loneliness"), col = 1),
  # ---- right column: exposure and biology panels ----
  list(title = "Socioeconomic position, work and sleep",  dom = c("A1", "B8", "A2"), k = 10,
       must = c("employment_unable_to_work_because_of_sickness_or_disability",
                "income", "employment_in_paid_employment_or_self_employed",
                "qualification_college_or_university_degree",
                "insomnia", "daytime_dozing"), col = 2),
  list(title = "Diet, taste preference and supplements",        dom = c("A3", "B1", "B2"), k = 8,
       must = c("diet_change_ill", "supp_vitamin_b12", "never_eat_dairy_products",
                "bread_type_wholemeal_or_wholegrain"), col = 2),
  list(title = "Body composition and blood biochemistry",       dom = c("A4", "A9"),   k = 12,
       must = c("body_fat_pct", "waist_circ", "hdl", "apoa", "urate", "hba1c", "crp",
                "height", "ldl"), col = 2),
  list(title = "Self-reported medical history",                 dom = c("A8"),         k = 5,
       drop_grp = c("taste_drug_user_baseline_numeric", "fluid_intel"),
       must = c("vascular_angina", "diabetes_sr"), col = 2),
  list(title = "Prior diagnoses: neurological, sensory, sinonasal",
       dom = c("C-Neuro", "C-Respiratory"), k = 9,
       must = c("phe332", "phe471", "phe475", "phe352"), col = 2),   # Parkinson / nasal polyps / chronic sinusitis / cranial nerve
  list(title = "Prior diagnoses: other body systems",
       dom = c("C-Digestive", "C-Circulatory", "C-Endocrine", "C-Mental"), k = 5, col = 2))

build_rows <- function(spec) {
  if (!is.null(spec$fixed)) {
    d <- ex[variable %in% spec$fixed][order(match(variable, spec$fixed))]
    n_sig <- d[pass_fdr == TRUE, .N]; n_tot <- length(spec$fixed)
  } else {
    pool <- ex[domain %in% spec$dom]
    cand <- pool[pass_fdr == TRUE][order(pval)]
    must <- if (is.null(spec$must)) character(0) else spec$must
    stopifnot(all(must %in% pool$variable))   # a misspelt name must fail here rather than silently drop a row
    ns_must <- setdiff(must, cand$variable)
    if (length(ns_must)) cand <- rbind(cand, pool[variable %in% ns_must])[order(pval)]
    cand[, rk := seq_len(.N), by = redgrp]
    if (!is.null(spec$drop_grp)) cand <- cand[!(redgrp %in% spec$drop_grp | variable %in% spec$drop_grp)]
    cand[, allow := fifelse(redgrp %in% names(PER_GRP_ALLOW), PER_GRP_ALLOW[redgrp], 1L)]
    pick <- cand[rk <= allow | variable %in% must]
    pick[, prio := as.integer(!(variable %in% must))]     # must-rows first
    d <- pick[order(prio, pval)][seq_len(min(spec$k, .N))][order(pval)]
    n_sig <- pool[pass_fdr == TRUE, .N]; n_tot <- nrow(pool)
  }
  d[, `:=`(panel = spec$title, colblk = spec$col, sig_frac = n_sig / n_tot,
           lab = pretty_of(variable, pretty_label),
           sig = pass_fdr == TRUE, hdr = sprintf("%d of %d", n_sig, n_tot))]
  d[]
}
D <- rbindlist(lapply(PANELS, build_rows), fill = TRUE)
cat(sprintf("\n[rows] %d rows / %d sub-panels | col1 %d rows, col2 %d rows\n",
            nrow(D), uniqueN(D$panel), D[colblk == 1, .N], D[colblk == 2, .N]))
print(D[, .(rows = .N, sig = sum(sig), domain = hdr[1]), by = panel])
cat(sprintf("[label] longest label %d chars; %d truncated\n",
            max(nchar(D$lab)), sum(grepl("\\.\\.$", D$lab))))

# Row dump: check against this file that every variable the Results names is
# actually on the figure, rather than reading the PDF by eye.
fwrite(D[, .(panel, colblk, variable, lab, OR = ORd, lo = LOd, hi = HId, pval, sig)],
       file.path(FIG_OUT, "_fig1_rows.csv"))

## ---- 5. draw ---------------------------------------------------------------
XLIM <- c(0.60, 2.40); X_CLIP <- 2.28
XBRK <- c(0.6, 0.8, 1.0, 1.5, 2.0)
COL_UP <- "#B23A48"; COL_DN <- "#4E79A7"; COL_NS <- "grey62"

D[, colr := fcase(!sig, COL_NS, ORd > 1, COL_UP, default = COL_DN)]
D[, out_hi := ORd > X_CLIP]
D[, `:=`(ORc = pmin(pmax(ORd, XLIM[1]), X_CLIP),
         LOc = pmin(pmax(LOd, XLIM[1]), X_CLIP),
         HIc = pmin(pmax(HId, XLIM[1]), X_CLIP))]
D[, oor_txt := fifelse(out_hi, sprintf("%.2f", ORd), NA_character_)]
cat(sprintf("[clip] %d rows exceed the axis and are drawn as triangles with the value printed: %s\n",
            D[out_hi == TRUE, .N], paste(D[out_hi == TRUE, lab], collapse = ", ")))

SZA <- SZ$anno / .pt

grDevices::pdfFonts(Arial        = grDevices::pdfFonts()$Helvetica)
grDevices::postscriptFonts(Arial = grDevices::postscriptFonts()$Helvetica)

PIE_D  <- 2.4   # diameter, mm
PIE_UP <- 2.4   # lift above the panel edge, mm
PIE_IN   <- 3.0   # keep at least this much clear of the right edge, mm
pie_grob <- function(frac, fill = "#4E79A7", bg = "grey88") {
  frac <- max(0, min(1, frac))
  th <- seq(pi / 2, pi / 2 - 2 * pi * frac, length.out = max(3L, round(72 * frac) + 2L))
  grid::grobTree(grid::grobTree(
    grid::circleGrob(r = unit(0.5, "npc"), gp = grid::gpar(fill = bg, col = NA)),
    if (frac > 0) grid::polygonGrob(
      x = unit(c(0.5, 0.5 + 0.5 * cos(th)), "npc"),
      y = unit(c(0.5, 0.5 + 0.5 * sin(th)), "npc"),
      gp = grid::gpar(fill = fill, col = NA)) else grid::nullGrob(),
    grid::circleGrob(r = unit(0.5, "npc"),
                     gp = grid::gpar(fill = NA, col = "grey45", lwd = 0.35)),
    vp = grid::viewport(x = unit(1, "npc") - unit(PIE_IN, "mm"),
                        y = unit(1, "npc") + unit(PIE_UP, "mm"),
                        width = unit(PIE_D, "mm"), height = unit(PIE_D, "mm"),
                        just = c("right", "centre"))))
}

sub_panel <- function(ttl, show_axis) {
  d <- D[panel == ttl][, y := rev(seq_len(.N))]
  ggplot(d, aes(x = ORc, y = y)) +
    # Only the last panel in a column draws an x axis; the others are read
    # against these faint rules, which stay aligned down the column.
    geom_vline(xintercept = XBRK, colour = "grey93", linewidth = 0.16) +
    geom_vline(xintercept = 1, linetype = "22", linewidth = 0.26, colour = "grey35") +
    geom_errorbar(aes(xmin = LOc, xmax = HIc), orientation = "y",
                  width = 0, linewidth = 0.40, colour = d$colr) +
    geom_point(aes(shape = out_hi), size = 1.25, colour = d$colr, fill = d$colr, stroke = 0.3) +
    geom_text(data = d[out_hi == TRUE][, vj := fifelse(y == max(d$y), 1.6, -0.9)],
              aes(x = X_CLIP, y = y, label = oor_txt, vjust = vj),
              hjust = 1.15, size = SZA * 0.88,
              family = FIG_FONT, colour = COL_UP) +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17), guide = "none") +
    scale_x_continuous(transform = "log", limits = XLIM, breaks = XBRK,
                       labels = function(x) format(x, drop0trailing = TRUE), expand = c(0, 0)) +
    scale_y_continuous(breaks = d$y, labels = d$lab,
                       limits = c(0.4, nrow(d) + 0.6), expand = c(0, 0)) +
    labs(x = if (show_axis) "Odds ratio per SD or per level" else NULL) +
    annotation_custom(pie_grob(d$sig_frac[1])) +
    coord_cartesian(clip = "off") +
    ggtitle(ttl) +
    theme_fig(grid = "none") +
    theme(axis.title.y = element_blank(), axis.line.y = element_blank(),
          axis.ticks.y = element_blank(), axis.text.y = element_text(hjust = 1),
          axis.title.x = if (show_axis) element_text(size = SZ$axtit) else element_blank(),
          axis.text.x  = if (show_axis) element_text(size = SZ$axtxt) else element_blank(),
          axis.ticks.x = if (show_axis) element_line(linewidth = LW$tick) else element_blank(),
          plot.title.position = "plot",
          plot.title = element_text(size = SZ$anno, face = "bold",
                                    family = FIG_FONT, hjust = 0, margin = margin(b = 0.1)),
          plot.subtitle = element_blank(),   # replaced by the pie chart at the right of the title
          plot.margin = margin(0.3, 1.5, 0.3, 1, "mm"))
}

## ---- canvas constants, aligned with Figure 3 -------------------------------------
CANVAS_MM  <- 205.0  # canvas height, mm; cap 207
# Weighting by row count already balances the two columns to 0.33%, so the
# spacer stays at zero.
CAL_SPACER <- 0.0

titles <- vapply(PANELS, function(x) x$title, character(1))
cols   <- vapply(PANELS, function(x) x$col,  numeric(1))
# Only the last panel in a column draws an x axis, saving nine repeated axes.
last_in_col <- vapply(seq_along(titles), function(i) i == max(which(cols == cols[i])), logical(1))
plots <- Map(sub_panel, titles, last_in_col)
# row height + title placeholder; panels with an axis add the axis height
nrows <- vapply(titles, function(x) D[panel == x, .N], numeric(1))
hgt   <- nrows

i1 <- which(cols == 1); i2 <- which(cols == 2)
h1 <- sum(hgt[i1]) + CAL_SPACER; h2 <- sum(hgt[i2]); H_MM <- CANVAS_MM
stopifnot(H_MM <= 207)                     # print frame
# The shorter column takes a spacer so both columns are the same height; the
# row pitch is equal by construction.
stack_col <- function(idx, extra) {
  ps <- plots[idx]; hs <- unname(hgt[idx])
  if (extra > 0.01) { ps <- c(ps, list(patchwork::plot_spacer())); hs <- c(hs, extra) }
  patchwork::wrap_plots(ps, ncol = 1, heights = hs)
}
fig <- stack_col(i1, CAL_SPACER) | stack_col(i2, 0)

cat(sprintf("[layout] col1 %d rows (+%.2f spacer) | col2 %d rows | canvas %.1f mm\n",
            D[colblk == 1, .N], CAL_SPACER, D[colblk == 2, .N], H_MM))
fig_save(fig, "fig1_exposure_forest", 168, H_MM)

cat("\n[legend] the caption must state:\n")
cat(sprintf("  1. Grey = not significant at 5%% FDR (%d rows): 4 hospital-recorded dental\n", D[sig == FALSE, .N]))
cat("     diagnoses, 4 self-reported denture wear, surgery and extraction, and height\n")
cat("     and LDL, the negative controls named in the Results.\n")
cat("  2. Triangle = the effect runs past the axis; the value is printed beside it.\n")
cat("  3. The pie beside each sub-panel title is the share of that domain's tests\n")
cat("     reaching FDR significance. It gives the proportion only, so the counts\n")
cat("     belong in the caption; the panel denominators do not sum to 234+107+454,\n")
cat("     because one domain of 206 tests is not shown.\n")
cat("  4. The oral panel is the fixed set of all 16 items, 8 of them not significant.\n")
cat("     Other panels select among their domain's FDR-significant results by\n")
cat("     clinical relevance and de-redundancy, not by the k smallest P values.\n")
cat("     Every variable the Results names is forced on by must.\n")
cat("  5. Continuous exposures are per SD; binary, ordinal and diagnosis rows are\n")
cat("     per level. The full scans are in the appendix.\n")
cat(sprintf("\n[check] row dump written to %s\n",
            file.path(FIG_OUT, "_fig1_rows.csv")))
