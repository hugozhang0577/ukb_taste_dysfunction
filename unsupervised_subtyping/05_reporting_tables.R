#!/usr/bin/env Rscript
# =============================================================================
# 05_reporting_tables.R  (subtype reporting tables, k=4)
#
# Produces:
#   (1) signature tables (per sex): per cluster, top-5 features (q<0.05) by
#       abs_effect; union defines rows; values = per-cluster median.
#   (2) label-rationale tables (per sex): per cluster, top-20 feature breakdown
#       by source_analysis.
#   (3) baseline Table 1: one wide table across subtypes x sexes, rows grouped by
#       biological domain; cells mean+/-SD (continuous) or n(%) (binary); per-row
#       omnibus p (Kruskal-Wallis / chi-square within sex), BH-adjusted.
#
# Subtype column ordering is fixed (A,B,C,D) via the sex-specific cluster mapping
# in _subtype_map.R (CODE_DIR, default = current dir; captured before setwd).
#
# Input:  subtyping/clusters/cluster_assignments_g1_{m,f}_k4.rds +
#         subtyping/reports/cluster_profiles_{m,f}_k4.csv +
#         subtyping/reports/fig_tables/table_cluster_topfeatures_{m,f}_k4.csv +
#         data/ml_ready/group1_full.rds + manifest/master_feature_manifest_final.csv
# Output: subtyping/reports/fig_tables/{pptx_signature_k4,pptx_label_rationale_k4}_{m,f}.csv +
#         baseline_table1_k4.csv
# =============================================================================

CODE_DIR <- Sys.getenv("CODE_DIR", unset = getwd())
rm_list_keep <- c("CODE_DIR")
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
suppressPackageStartupMessages(library(data.table))

P7_DIR  <- "output/subtyping"
RPT_DIR <- file.path(P7_DIR, "reports")
TBL_DIR <- file.path(RPT_DIR, "fig_tables")
ML_DIR  <- "output/ml_ready"
MF_DIR  <- "output/feature_manifest"
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)
TOP_N <- 5; Q_THR <- 0.05

source(file.path(CODE_DIR, "_subtype_map.R"))   # SUBTYPE_MAP / SUBTYPE_MAP_INV
SUBTYPE_MAP <- SUBTYPE_MAP_INV                       # letter -> cluster_id
SUBTYPE_LABEL <- c(A = "Aging frailty", B = "Psychosomatic", C = "Cardiometabolic", D = "Young idiopathic")

while (sink.number() > 0) sink()
sink(file.path(RPT_DIR, "07_reporting_tables_k4_log.txt"), split = TRUE)
on.exit(sink(), add = TRUE)
cat("=== 07: reporting tables (k=4) ===\n")

manifest  <- fread(file.path(MF_DIR, "master_feature_manifest_final.csv"))
cases_all <- as.data.table(readRDS(file.path(ML_DIR, "group1_full.rds")))[taste_2w_strict == 1]

# ---- [1] signature + rationale per sex --------------------------------------
build_pptx_and_rationale <- function(sex) {
  suf <- tolower(sex)
  cl_obj <- readRDS(file.path(P7_DIR, "clusters", sprintf("cluster_assignments_g1_%s_k4.rds", suf)))
  prof   <- fread(file.path(RPT_DIR, sprintf("cluster_profiles_%s_k4.csv", suf)))
  topfn  <- fread(file.path(TBL_DIR, sprintf("table_cluster_topfeatures_%s_k4.csv", suf)))
  K <- cl_obj$final_k
  eids <- as.integer(rownames(cl_obj$Z_active)); lab <- cl_obj$final_labels
  cases <- cases_all[match(eids, cases_all$eid)]; cases[, cluster := factor(lab, levels = 1:K)]

  picks <- prof[q < Q_THR][order(cluster, -abs_effect), head(.SD, TOP_N), by = cluster]
  feat_union <- unique(picks$feature_id)
  med_long <- rbindlist(lapply(feat_union, function(f) {
    if (is.null(cases[[f]])) return(NULL)
    cases[, .(median_val = round(median(get(f), na.rm = TRUE), 3)), by = cluster][, feature_id := f][]
  }))
  wide <- dcast(med_long, feature_id ~ cluster, value.var = "median_val")
  setnames(wide, as.character(1:K), paste0("c", 1:K))
  vals <- as.matrix(wide[, paste0("c", 1:K), with = FALSE])
  dev_mat <- abs(vals - rowMeans(vals, na.rm = TRUE))
  wide[, highlight_cluster := paste0("c", apply(dev_mat, 1, function(r) if (all(is.na(r))) NA_integer_ else which.max(r)))]
  wide <- merge(wide, manifest[, .(feature_id, source_analysis, description)], by = "feature_id", all.x = TRUE)
  max_eff <- prof[q < Q_THR, .(max_eff = max(abs_effect), min_q = min(q)), by = feature_id]
  wide <- merge(wide, max_eff, by = "feature_id", all.x = TRUE)[order(highlight_cluster, -max_eff)]
  setcolorder(wide, c("feature_id", "source_analysis", paste0("c", 1:K), "highlight_cluster", "max_eff", "min_q", "description"))
  fwrite(wide, file.path(TBL_DIR, sprintf("pptx_signature_k4_%s.csv", suf)))

  rationale <- merge(topfn[, .(cluster, feature_id)], manifest[, .(feature_id, source_analysis)], by = "feature_id", all.x = TRUE)
  rationale <- rationale[, .(n = .N), by = .(cluster, source_analysis)][order(cluster, -n)]
  fwrite(rationale, file.path(TBL_DIR, sprintf("pptx_label_rationale_k4_%s.csv", suf)))
  invisible(list(cases = cases, K = K))
}
pm <- build_pptx_and_rationale("M")
pf <- build_pptx_and_rationale("F")

# ---- [2] baseline Table 1 ---------------------------------------------------
FEAT_SPEC <- rbind(
  data.table(domain = "Demographics",        feature_id = "age",                 display = "Age at taste assessment", unit = "years",    type = "cont"),
  data.table(domain = "Demographics",        feature_id = "age_baseline",        display = "Age at UKB baseline",     unit = "years",    type = "cont"),
  data.table(domain = "Body composition",    feature_id = "BMI",                 display = "BMI",                     unit = "kg/m^2",   type = "cont"),
  data.table(domain = "Body composition",    feature_id = "waist_circ",          display = "Waist circumference",     unit = "cm",       type = "cont"),
  data.table(domain = "Body composition",    feature_id = "whole_body_fat",      display = "Whole-body fat mass",     unit = "kg",       type = "cont"),
  data.table(domain = "Body composition",    feature_id = "body_fat_pct",        display = "Body fat",                unit = "%",        type = "cont"),
  data.table(domain = "Body composition",    feature_id = "grip_max",            display = "Handgrip, max",           unit = "kg",       type = "cont"),
  data.table(domain = "Metabolic / hepatic", feature_id = "urate",               display = "Urate",                   unit = "umol/L",   type = "cont"),
  data.table(domain = "Metabolic / hepatic", feature_id = "shbg",                display = "SHBG",                    unit = "nmol/L",   type = "cont"),
  data.table(domain = "Metabolic / hepatic", feature_id = "ggt",                 display = "GGT",                     unit = "U/L",      type = "cont"),
  data.table(domain = "Metabolic / hepatic", feature_id = "hba1c",               display = "HbA1c",                   unit = "mmol/mol", type = "cont"),
  data.table(domain = "Metabolic / hepatic", feature_id = "vitamin_d",           display = "25-OH vitamin D",         unit = "nmol/L",   type = "cont"),
  data.table(domain = "Psychiatric",         feature_id = "phq9_total",          display = "PHQ-9 (depression)",      unit = "0-27",     type = "cont"),
  data.table(domain = "Psychiatric",         feature_id = "gad7_total",          display = "GAD-7 (anxiety)",         unit = "0-21",     type = "cont"),
  data.table(domain = "Psychiatric",         feature_id = "neuroticism",         display = "Neuroticism (EPQ-R-S)",   unit = "0-12",     type = "cont"),
  data.table(domain = "Psychiatric",         feature_id = "phq15_somatic_total", display = "PHQ-15 somatic symptoms", unit = "0-30",     type = "cont"),
  data.table(domain = "Comorbidity / meds",  feature_id = "dx_count_total",      display = "Total PheCode diagnoses", unit = "count",    type = "cont"),
  data.table(domain = "Comorbidity / meds",  feature_id = "n_medications",       display = "Number of medications",   unit = "count",    type = "cont"),
  fill = TRUE
)
summarise_cell <- function(x, type) {
  x <- x[!is.na(x)]; if (length(x) < 5) return("--")
  if (type == "cont") sprintf("%.1f +/- %.1f", mean(x), sd(x))
  else { n <- sum(x == 1); sprintf("%d (%.1f%%)", n, 100 * n / length(x)) }
}
omnibus_p <- function(values, labels, type) {
  ok <- !is.na(values) & !is.na(labels); if (sum(ok) < 20) return(NA_real_)
  v <- values[ok]; g <- as.factor(labels[ok]); if (nlevels(g) < 2) return(NA_real_)
  if (type == "cont") tryCatch(kruskal.test(v ~ g)$p.value, error = function(e) NA_real_)
  else { tab <- table(v, g); if (any(dim(tab) < 2)) return(NA_real_)
    tryCatch(suppressWarnings(chisq.test(tab)$p.value), error = function(e) NA_real_) }
}
build_row <- function(spec_row) {
  fid <- spec_row$feature_id; type <- spec_row$type; cells <- list(); pvals <- c(M = NA_real_, F = NA_real_)
  for (sex in c("M", "F")) {
    cases <- (if (sex == "M") pm else pf)$cases; x <- cases[[fid]]
    if (is.null(x)) { for (a in c("A", "B", "C", "D")) cells[[paste(sex, a, sep = "_")]] <- "--"; next }
    pvals[[sex]] <- omnibus_p(x, cases$cluster, type)
    amap <- SUBTYPE_MAP[[tolower(sex)]]
    for (a in c("A", "B", "C", "D")) cells[[paste(sex, a, sep = "_")]] <- summarise_cell(x[cases$cluster == amap[[a]]], type)
  }
  data.table(domain = spec_row$domain,
    display = sprintf("%s%s", spec_row$display, ifelse(nzchar(spec_row$unit), sprintf(" (%s)", spec_row$unit), "")),
    feature_id = fid, M_A = cells$M_A, M_B = cells$M_B, M_C = cells$M_C, M_D = cells$M_D,
    F_A = cells$F_A, F_B = cells$F_B, F_C = cells$F_C, F_D = cells$F_D, p_M = pvals[["M"]], p_F = pvals[["F"]])
}
rows <- rbindlist(lapply(seq_len(nrow(FEAT_SPEC)), function(i) build_row(FEAT_SPEC[i])))
rows[, p_M_BH := p.adjust(p_M, method = "BH")][, p_F_BH := p.adjust(p_F, method = "BH")]
fmtp <- function(p) ifelse(is.na(p), "--", ifelse(p < 1e-3, "<0.001", formatC(p, format = "f", digits = 3)))
rows[, p_M_display := fmtp(p_M_BH)][, p_F_display := fmtp(p_F_BH)]

n_by <- list()
for (sex in c("M", "F")) {
  tb <- table((if (sex == "M") pm else pf)$cases$cluster); amap <- SUBTYPE_MAP[[tolower(sex)]]
  for (a in c("A", "B", "C", "D")) n_by[[paste(sex, a, sep = "_")]] <- tb[[as.character(amap[[a]])]]
}
mk_hdr <- function(display, fid, vals) data.table(domain = "_header", display = display, feature_id = fid,
  M_A = vals[1], M_B = vals[2], M_C = vals[3], M_D = vals[4], F_A = vals[1], F_B = vals[2], F_C = vals[3], F_D = vals[4],
  p_M = NA_real_, p_F = NA_real_, p_M_BH = NA_real_, p_F_BH = NA_real_, p_M_display = "", p_F_display = "")
hdr_label <- mk_hdr("Subtype label", "_label", SUBTYPE_LABEL[c("A", "B", "C", "D")])
hdr_cluster <- data.table(domain = "_header", display = "Cluster ID (sex-specific)", feature_id = "_cid",
  M_A = paste0("M-c", SUBTYPE_MAP$m["A"]), M_B = paste0("M-c", SUBTYPE_MAP$m["B"]), M_C = paste0("M-c", SUBTYPE_MAP$m["C"]),
  M_D = paste0("M-c", SUBTYPE_MAP$m["D"]), F_A = paste0("F-c", SUBTYPE_MAP$f["A"]), F_B = paste0("F-c", SUBTYPE_MAP$f["B"]),
  F_C = paste0("F-c", SUBTYPE_MAP$f["C"]), F_D = paste0("F-c", SUBTYPE_MAP$f["D"]), p_M = NA_real_, p_F = NA_real_,
  p_M_BH = NA_real_, p_F_BH = NA_real_, p_M_display = "", p_F_display = "")
hdr_n <- data.table(domain = "_header", display = "N per subtype", feature_id = "_n",
  M_A = sprintf("n=%d", n_by$M_A), M_B = sprintf("n=%d", n_by$M_B), M_C = sprintf("n=%d", n_by$M_C), M_D = sprintf("n=%d", n_by$M_D),
  F_A = sprintf("n=%d", n_by$F_A), F_B = sprintf("n=%d", n_by$F_B), F_C = sprintf("n=%d", n_by$F_C), F_D = sprintf("n=%d", n_by$F_D),
  p_M = NA_real_, p_F = NA_real_, p_M_BH = NA_real_, p_F_BH = NA_real_, p_M_display = "", p_F_display = "")
out <- rbind(hdr_label, hdr_cluster, hdr_n, rows, fill = TRUE)
setcolorder(out, c("domain", "display", "feature_id", "M_A", "M_B", "M_C", "M_D", "F_A", "F_B", "F_C", "F_D",
                   "p_M_display", "p_F_display", "p_M", "p_F", "p_M_BH", "p_F_BH"))
fwrite(out, file.path(TBL_DIR, "baseline_table1_k4.csv"), bom = TRUE)
cat(sprintf("  baseline_table1_k4.csv written (%d rows)\n", nrow(out)))
cat(sprintf("subtype cluster mapping: M A=c%s B=c%s C=c%s D=c%s | F A=c%s B=c%s C=c%s D=c%s\n",
            SUBTYPE_MAP$m["A"], SUBTYPE_MAP$m["B"], SUBTYPE_MAP$m["C"], SUBTYPE_MAP$m["D"],
            SUBTYPE_MAP$f["A"], SUBTYPE_MAP$f["B"], SUBTYPE_MAP$f["C"], SUBTYPE_MAP$f["D"]))
cat("=== 07 reporting tables complete ===\n")
