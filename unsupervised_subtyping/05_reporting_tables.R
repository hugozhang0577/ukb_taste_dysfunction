# Baseline characteristics of the four subtypes by sex (mean +/- SD); Kruskal-Wallis test across subtypes
# within sex, BH-adjusted over the table rows.
library(data.table)
source("unsupervised_subtyping/_subtype_map.R")
setwd(Sys.getenv("PROJECT_DIR", "."))

cases_all <- as.data.table(readRDS("output/ml_ready/group1_full.rds"))[taste_2w_strict == 1]
cases_of <- function(sx) {
  cl <- readRDS(sprintf("output/subtyping/clusters/cluster_assignments_g1_%s_k4.rds", sx))
  cases_all[match(as.integer(rownames(cl$Z_active)), cases_all$eid)][, cluster := factor(cl$final_labels, levels = 1:cl$final_k)][]
}
CASES <- list(m = cases_of("m"), f = cases_of("f"))

FEAT_SPEC <- data.table(
  domain = rep(c("Demographics", "Body composition", "Metabolic / hepatic", "Psychiatric", "Comorbidity / meds"),
               c(2, 5, 5, 4, 2)),
  feature_id = c("age", "age_baseline", "BMI", "waist_circ", "whole_body_fat", "body_fat_pct", "grip_max",
                 "urate", "shbg", "ggt", "hba1c", "vitamin_d", "phq9_total", "gad7_total", "neuroticism",
                 "phq12_somatic_total", "dx_count_total", "n_medications"),
  display = c("Age at taste assessment", "Age at UKB baseline", "BMI", "Waist circumference", "Whole-body fat mass",
              "Body fat", "Handgrip, max", "Urate", "SHBG", "GGT", "HbA1c", "25-OH vitamin D", "PHQ-9 (depression)",
              "GAD-7 (anxiety)", "Neuroticism score", "Somatic symptoms (4-item)", "Total recorded diagnoses",
              "Number of medications"),
  unit = c("years", "years", "kg/m^2", "cm", "kg", "%", "kg", "umol/L", "nmol/L", "U/L", "mmol/mol", "nmol/L",
           "0-27", "0-21", "0-12", "0-8", "count", "count"))

cell <- function(x) { x <- x[!is.na(x)]; if (length(x) < 5) "--" else sprintf("%.1f +/- %.1f", mean(x), sd(x)) }
omnibus_p <- function(v, g) {
  ok <- !is.na(v) & !is.na(g); if (sum(ok) < 20) return(NA_real_)
  v <- v[ok]; g <- as.factor(g[ok]); if (nlevels(g) < 2) return(NA_real_)
  tryCatch(kruskal.test(v ~ g)$p.value, error = function(e) NA_real_)
}
COLS <- c(paste0("M_", c("A", "B", "C", "D")), paste0("F_", c("A", "B", "C", "D")))
rows <- rbindlist(lapply(seq_len(nrow(FEAT_SPEC)), function(i) {
  s <- FEAT_SPEC[i]
  r <- data.table(domain = s$domain, display = sprintf("%s (%s)", s$display, s$unit), feature_id = s$feature_id)
  for (sx in c("m", "f")) {
    cs <- CASES[[sx]]; x <- cs[[s$feature_id]]; SX <- toupper(sx)
    for (a in c("A", "B", "C", "D"))
      r[, (paste0(SX, "_", a)) := if (is.null(x)) "--" else cell(x[cs$cluster == SUBTYPE_MAP_INV[[sx]][[a]]])]
    r[, (paste0("p_", SX)) := if (is.null(x)) NA_real_ else omnibus_p(x, cs$cluster)]
  }
  r
}), use.names = TRUE)
fmtp <- function(p) ifelse(is.na(p), "--", ifelse(p < 1e-3, "<0.001", formatC(p, format = "f", digits = 3)))
rows[, `:=`(p_M_BH = p.adjust(p_M, method = "BH"), p_F_BH = p.adjust(p_F, method = "BH"))]
rows[, `:=`(p_M_display = fmtp(p_M_BH), p_F_display = fmtp(p_F_BH))]

# header rows: subtype name, sex-specific cluster id, n per subtype
hdr <- function(display, fid, vals)
  data.table(domain = "_header", display = display, feature_id = fid, p_M_display = "", p_F_display = "")[
    , (COLS) := as.list(vals)][]
n_of <- function(sx) as.integer(table(CASES[[sx]]$cluster)[as.character(SUBTYPE_MAP_INV[[sx]])])
out <- rbind(hdr("Subtype label", "_label", rep(SUBTYPE_NAME[c("A", "B", "C", "D")], 2)),
             hdr("Cluster ID (sex-specific)", "_cid", c(paste0("M-c", SUBTYPE_MAP_INV$m), paste0("F-c", SUBTYPE_MAP_INV$f))),
             hdr("N per subtype", "_n", sprintf("n=%d", c(n_of("m"), n_of("f")))),
             rows, fill = TRUE)
setcolorder(out, c("domain", "display", "feature_id", COLS, "p_M_display", "p_F_display", "p_M", "p_F", "p_M_BH", "p_F_BH"))
fwrite(out, "output/subtyping/reports/fig_tables/baseline_table1_k4.csv", bom = TRUE)
