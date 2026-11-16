# Subgroup (fairness) audit of the Platt-calibrated predictions (group1 out-of-fold, group2/3
# held-out): discrimination and calibration by cohort, sex, age band, Townsend quintile and
# reported smell change; gaps between the best and worst subgroup.
library(data.table)
library(pROC)
setwd(Sys.getenv("PROJECT_DIR", "."))

ML_DIR    <- "output/ml_ready"
MODEL_DIR <- "output/models/xgboost"
EVAL_DIR  <- "output/model_reports/eval"
EPS <- 1e-15; MIN_CASES <- 20L
models <- fread("output/model_reports/xgb_best_params.csv")$model_id

demo <- rbindlist(lapply(1:3, function(g)
  readRDS(file.path(ML_DIR, sprintf("group%d_full.rds", g)))[
    , .(eid, sex, age, townsend, smell_any, cohort = paste0("g", g))]))
demo[, age_band := cut(age, breaks = c(-Inf, 55, 65, Inf), labels = c("<55", "55-64", "65+"), right = FALSE)]
demo[, townsend_q := cut(townsend, breaks = quantile(townsend, seq(0, 1, 0.2), na.rm = TRUE),
                         include.lowest = TRUE, labels = paste0("Q", 1:5))]
demo[, sex_label := fifelse(sex == 1, "male", fifelse(sex == 0, "female", NA_character_))]
demo[, smell_stratum := fifelse(smell_any == 1, "has_smell_issue",
                                fifelse(smell_any == 0, "no_smell_issue", NA_character_))]
STRATA <- c(cohort = "cohort", sex = "sex_label", age_band = "age_band",
            townsend_q = "townsend_q", smell_stratum = "smell_stratum")

clip <- function(p) pmax(pmin(p, 1 - EPS), EPS)
# Metrics only where a subgroup has >= 50 participants and >= 20 cases and controls
subgroup_metrics <- function(y, p) {
  out <- data.table(n = length(y), n_cases = sum(y), auc = NA_real_, auc_lo = NA_real_, auc_hi = NA_real_,
                    brier = NA_real_, cal_intercept = NA_real_, cal_slope = NA_real_,
                    sens_at_spec90 = NA_real_, ppv_top10 = NA_real_, case_rate = NA_real_)
  if (length(y) < 50 || sum(y) < MIN_CASES || sum(y == 0) < MIN_CASES) return(out)
  roc_obj <- pROC::roc(y, p, quiet = TRUE, direction = "<")
  auc_ci  <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
  lp  <- qlogis(clip(p))
  sel <- p >= quantile(p, 0.9, na.rm = TRUE)
  out[, `:=`(auc = auc_ci[2], auc_lo = auc_ci[1], auc_hi = auc_ci[3], brier = mean((p - y)^2),
             cal_intercept = unname(coef(glm(y ~ offset(lp), family = binomial()))[1]),
             cal_slope = unname(coef(glm(y ~ lp, family = binomial()))[2]),
             sens_at_spec90 = as.numeric(pROC::coords(roc_obj, x = 0.9, input = "specificity",
                                                      ret = "sensitivity", transpose = FALSE)$sensitivity[1]),
             ppv_top10 = sum(y == 1 & sel) / max(sum(sel), 1), case_rate = mean(y))]
}

load_preds <- function(mid) rbindlist(c(
  list(readRDS(file.path(MODEL_DIR, sprintf("%s_oof_calibrated.rds", mid)))[, .(eid, y, pred = pred_platt, cohort = "g1")]),
  lapply(c("g2", "g3"), function(ch)
    readRDS(file.path(MODEL_DIR, "predictions", sprintf("%s_%s_preds.rds", mid, ch)))[
      , .(eid, y, pred = pred_platt, cohort = ch)])))

all_metrics <- list(); gap_rows <- list()
for (mid in models) {
  dat <- merge(load_preds(mid), demo, by = c("eid", "cohort"), all.x = TRUE)
  mdt <- rbindlist(lapply(names(STRATA), function(s) {
    col <- STRATA[[s]]; lv_all <- unique(dat[[col]])
    rbindlist(lapply(as.character(lv_all[!is.na(lv_all)]), function(lv) {
      sub <- dat[get(col) == lv & !is.na(pred) & !is.na(y)]
      cbind(data.table(model_id = mid, subgroup_var = s, subgroup_level = lv),
            subgroup_metrics(sub$y, sub$pred))
    }))
  }))
  all_metrics[[mid]] <- mdt
  for (s in names(STRATA)) {
    g <- mdt[subgroup_var == s & !is.na(auc)]
    if (nrow(g) < 2) next
    gap_rows[[paste(mid, s)]] <- data.table(model_id = mid, subgroup_var = s, n_levels = nrow(g),
      auc_min = min(g$auc), auc_max = max(g$auc), auc_gap = max(g$auc) - min(g$auc),
      ppv10_min = min(g$ppv_top10, na.rm = TRUE), ppv10_max = max(g$ppv_top10, na.rm = TRUE),
      ppv10_gap = max(g$ppv_top10, na.rm = TRUE) - min(g$ppv_top10, na.rm = TRUE),
      sens90_gap = max(g$sens_at_spec90, na.rm = TRUE) - min(g$sens_at_spec90, na.rm = TRUE))
  }
}
fwrite(rbindlist(all_metrics, fill = TRUE), file.path(EVAL_DIR, "fairness_subgroup_metrics.csv"))
fwrite(rbindlist(gap_rows, fill = TRUE), file.path(EVAL_DIR, "fairness_gaps_summary.csv"))
