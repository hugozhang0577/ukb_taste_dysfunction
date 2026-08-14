#!/usr/bin/env Rscript
# =============================================================================
# 08_fairness_audit.R  (subgroup / fairness audit)
#
# Pool calibrated predictions from G1 OOF (06h) + G2/G3 (06i), join demographic
# covariates by eid, and stratify by cohort (ancestry proxy), sex, age band, and
# Townsend SES quintile. Per subgroup: AUC (95% DeLong CI), calibration
# intercept/slope, Brier, sens@spec=0.9, PPV@top10%, case rate, N. Per-model
# fairness-gap summary (max-min across subgroups) + a forest plot for M1_TierA.
#
# Input:  models/xgboost/{mid}_oof_calibrated.rds + predictions/{mid}_{g2,g3}_preds.rds
#         data/ml_ready/group{1,2,3}_full.rds (demographics)
# Output: model_reports/eval/fairness_subgroup_metrics.csv + fairness_gaps_summary.csv
#         model_reports/eval/plots/M1_TierA_fairness_forest.png + fairness_audit_log.txt
# =============================================================================

rm(list = ls()); gc()
while (sink.number() > 0) sink()
suppressMessages({ library(data.table); library(pROC) })

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)

ML_DIR    <- "output/ml_ready"
MODEL_DIR <- "output/models/xgboost"
PRED_DIR  <- file.path(MODEL_DIR, "predictions")
P6_REP    <- "output/model_reports"
EVAL_DIR  <- file.path(P6_REP, "eval")
PLOT_DIR  <- file.path(EVAL_DIR, "plots")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)
sink(file.path(P6_REP, "fairness_audit_log.txt"), split = TRUE)
on.exit({ while (sink.number() > 0) sink() }, add = TRUE)

cat("=== fairness / subgroup audit ===\n")
EPS <- 1e-15; MIN_CASES <- 20L; PLOT_MODEL <- "M1_TierA"
models <- fread(file.path(P6_REP, "xgb_best_params.csv"))$model_id

# ---- demographics -----------------------------------------------------------
load_demo <- function(group_num) {
  dt <- as.data.table(readRDS(file.path(ML_DIR, sprintf("group%d_full.rds", group_num))))
  eid_col <- intersect(c("eid", "eid_code", "ID"), names(dt))[1]
  keep <- c(eid_col, intersect(c("sex", "age", "townsend", "townsend_deprivation_index", "smell_any"), names(dt)))
  dt <- dt[, keep, with = FALSE]; setnames(dt, eid_col, "eid")
  if ("townsend_deprivation_index" %in% names(dt)) setnames(dt, "townsend_deprivation_index", "townsend")
  dt[, cohort := paste0("g", group_num)]; dt
}
demo <- rbindlist(list(load_demo(1), load_demo(2), load_demo(3)), fill = TRUE, use.names = TRUE)
demo[, age_band := cut(age, breaks = c(-Inf, 55, 65, Inf), labels = c("<55", "55-64", "65+"), right = FALSE)]
demo[, townsend_q := cut(townsend, breaks = quantile(townsend, seq(0, 1, 0.2), na.rm = TRUE),
                         include.lowest = TRUE, labels = paste0("Q", 1:5))]
demo[, sex_label := fifelse(sex == 1, "male", fifelse(sex == 0, "female", NA_character_))]
if ("smell_any" %in% names(demo))
  demo[, smell_stratum := fifelse(smell_any == 1, "has_smell_issue", fifelse(smell_any == 0, "no_smell_issue", NA_character_))]
cat(sprintf("demographics: N=%d (G1=%d G2=%d G3=%d)\n", nrow(demo),
            sum(demo$cohort == "g1"), sum(demo$cohort == "g2"), sum(demo$cohort == "g3")))

# ---- helpers ----------------------------------------------------------------
clip <- function(p) pmax(pmin(p, 1 - EPS), EPS); safe_logit <- function(p) qlogis(clip(p))
subgroup_metrics <- function(y, p) {
  na_row <- data.table(n = length(y), n_cases = sum(y), auc = NA_real_, auc_lo = NA_real_,
    auc_hi = NA_real_, brier = NA_real_, cal_intercept = NA_real_, cal_slope = NA_real_,
    sens_at_spec90 = NA_real_, ppv_top10 = NA_real_, case_rate = NA_real_)
  if (length(y) < 50 || sum(y) < MIN_CASES || sum(y == 0) < MIN_CASES) return(na_row)
  roc_obj <- pROC::roc(y, p, quiet = TRUE, direction = "<")
  auc_ci <- tryCatch(as.numeric(pROC::ci.auc(roc_obj, method = "delong")), error = function(e) rep(NA_real_, 3))
  lp <- safe_logit(p)
  cal_int <- tryCatch(unname(coef(glm(y ~ offset(lp), family = binomial()))[1]), error = function(e) NA_real_)
  cal_slp <- tryCatch(unname(coef(glm(y ~ lp, family = binomial()))[2]), error = function(e) NA_real_)
  ss <- tryCatch(pROC::coords(roc_obj, x = 0.9, input = "specificity", ret = "sensitivity", transpose = FALSE),
                 error = function(e) data.frame(sensitivity = NA_real_))
  sel <- p >= quantile(p, 0.9, na.rm = TRUE)
  data.table(n = length(y), n_cases = sum(y), auc = auc_ci[2], auc_lo = auc_ci[1], auc_hi = auc_ci[3],
    brier = mean((p - y)^2), cal_intercept = cal_int, cal_slope = cal_slp,
    sens_at_spec90 = as.numeric(ss$sensitivity[1]), ppv_top10 = sum(y == 1 & sel) / max(sum(sel), 1),
    case_rate = mean(y))
}
load_preds <- function(mid) {
  g1 <- as.data.table(readRDS(file.path(MODEL_DIR, sprintf("%s_oof_calibrated.rds", mid))))
  g1 <- g1[, .(eid, y, pred = pred_platt)]; g1[, cohort := "g1"]; out <- list(g1)
  for (ch in c("g2", "g3")) {
    f <- file.path(PRED_DIR, sprintf("%s_%s_preds.rds", mid, ch)); if (!file.exists(f)) next
    x <- as.data.table(readRDS(f))[, .(eid, y, pred = pred_platt)]; x[, cohort := ch]
    out[[length(out) + 1]] <- x
  }
  rbindlist(out, fill = TRUE)
}

# ---- main loop --------------------------------------------------------------
all_metrics <- list(); gap_rows <- list()
for (mid in models) {
  cat(strrep("-", 60), "\nModel:", mid, "\n", sep = "")
  dat <- merge(load_preds(mid), demo, by = c("eid", "cohort"), all.x = TRUE)
  strata <- list(cohort = "cohort", sex = "sex_label", age_band = "age_band", townsend_q = "townsend_q")
  if ("smell_stratum" %in% names(dat)) strata[["smell_stratum"]] <- "smell_stratum"
  mm <- list()
  for (s_name in names(strata)) {
    s_col <- strata[[s_name]]; levels_s <- unique(dat[[s_col]]); levels_s <- levels_s[!is.na(levels_s)]
    for (lv in levels_s) {
      sub <- dat[get(s_col) == lv & !is.na(pred) & !is.na(y)]
      m <- subgroup_metrics(sub$y, sub$pred)
      m[, `:=`(model_id = mid, subgroup_var = s_name, subgroup_level = as.character(lv))]
      mm[[paste(s_name, lv, sep = "::")]] <- m
    }
  }
  mdt <- rbindlist(mm); setcolorder(mdt, c("model_id", "subgroup_var", "subgroup_level", "n", "n_cases"))
  all_metrics[[mid]] <- mdt
  for (s_name in names(strata)) {
    g <- mdt[subgroup_var == s_name & !is.na(auc)]; if (nrow(g) < 2) next
    gap_rows[[paste(mid, s_name, sep = "_")]] <- data.table(model_id = mid, subgroup_var = s_name,
      n_levels = nrow(g), auc_min = min(g$auc), auc_max = max(g$auc), auc_gap = max(g$auc) - min(g$auc),
      ppv10_min = min(g$ppv_top10, na.rm = TRUE), ppv10_max = max(g$ppv_top10, na.rm = TRUE),
      ppv10_gap = max(g$ppv_top10, na.rm = TRUE) - min(g$ppv_top10, na.rm = TRUE),
      sens90_gap = max(g$sens_at_spec90, na.rm = TRUE) - min(g$sens_at_spec90, na.rm = TRUE))
  }
}
metrics_dt <- rbindlist(all_metrics, fill = TRUE)
gaps_dt <- rbindlist(gap_rows, fill = TRUE)
fwrite(metrics_dt, file.path(EVAL_DIR, "fairness_subgroup_metrics.csv"))
fwrite(gaps_dt, file.path(EVAL_DIR, "fairness_gaps_summary.csv"))

# ---- forest plot (M1_TierA) -------------------------------------------------
fm <- metrics_dt[model_id == PLOT_MODEL & !is.na(auc)]
if (nrow(fm) > 0) {
  png(file.path(PLOT_DIR, sprintf("%s_fairness_forest.png", PLOT_MODEL)), 1200, 800, res = 150)
  par(mar = c(4, 11, 3, 2)); y_pos <- seq_len(nrow(fm))
  plot(fm$auc, y_pos, xlim = c(0.6, 1.0), ylim = c(0.5, nrow(fm) + 0.5), pch = 19, col = "#1f77b4",
       yaxt = "n", xlab = "AUC (95% CI)", ylab = "", main = sprintf("%s - subgroup AUC (fairness audit)", PLOT_MODEL))
  segments(fm$auc_lo, y_pos, fm$auc_hi, y_pos, col = "#1f77b4", lwd = 2)
  axis(2, at = y_pos, labels = paste(fm$subgroup_var, fm$subgroup_level, sep = ": "), las = 1, cex.axis = 0.8)
  abline(v = mean(fm[subgroup_var == "cohort"]$auc, na.rm = TRUE), lty = 2, col = "grey50")
  dev.off()
}
cat("\n=== 06 complete ===\nNext: 07_sensitivity_no_smell.R\n")
