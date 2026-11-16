# SARS-CoV-2 positivity before the taste assessment among tested cases: each subtype vs the other three by logistic
# regression adjusted for sex, then also for age (linear, and natural spline with 3 df); raw positivity by age band.
library(data.table)
library(splines)
source("unsupervised_subtyping/_subtype_map.R")
setwd(Sys.getenv("PROJECT_DIR", "."))
OUT_DIR <- "output/subtyping/evidence_covid"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

labels_of <- function(sx) {
  fl <- readRDS(sprintf("output/subtyping/clusters/cluster_assignments_g1_%s_k4.rds", sx))$final_labels
  data.table(eid = as.integer(names(fl)), sex = toupper(sx), subtype = unname(SUBTYPE_MAP[[sx]][as.character(as.integer(fl))]))
}
d <- merge(merge(rbind(labels_of("m"), labels_of("f")), fread("input/analysis_ready/covid_temporal_flags.csv"), by = "eid"),
           as.data.table(readRDS("output/ml_ready/group1_full.rds"))[, .(eid, age)], by = "eid")
tst <- d[!is.na(phe_ever_positive) & !is.na(age) & !is.na(phe_pos_before_taste)]   # cases with a SARS-CoV-2 test record
tst[, sex := factor(sex)]

ovr <- function(letter, rhs) {
  m  <- glm(as.formula(paste("phe_pos_before_taste ~ grp +", rhs)),
            data = copy(tst)[, grp := as.integer(subtype == letter)], family = binomial())
  s  <- summary(m)$coefficients["grp", ]
  ci <- exp(s[1] + c(-1, 1) * 1.96 * s[2])
  data.table(subtype = letter, OR = round(exp(s[1]), 3), lo = round(ci[1], 3), hi = round(ci[2], 3), P = signif(s[4], 3))
}
MODELS <- c(sex_only = "sex", age_sex = "age + sex", ns_age3_sex = "ns(age, 3) + sex")
res <- rbindlist(lapply(names(MODELS), function(mn)
  rbindlist(lapply(c("A", "B", "C", "D"), ovr, rhs = MODELS[[mn]]))[, model := mn]))
res <- merge(res, tst[, .(n_tested = .N, n_pos_before = sum(phe_pos_before_taste),
                          raw_pos_pct = round(100 * mean(phe_pos_before_taste), 1)), by = subtype], by = "subtype")
res[, subtype_name := SUBTYPE_NAME[subtype]]
setcolorder(res, c("model", "subtype", "subtype_name", "OR", "lo", "hi", "P", "n_tested", "n_pos_before", "raw_pos_pct"))
fwrite(res[order(factor(model, levels = names(MODELS)), subtype)], file.path(OUT_DIR, "covid_onevsrest_OR.csv"))

tst[, ageband := cut(age, c(0, 60, 65, 70, 75, 200), labels = c("<60", "60-65", "65-70", "70-75", "75+"), right = FALSE)]
fwrite(dcast(tst[, .(pos_pct = round(100 * mean(phe_pos_before_taste), 0)), by = .(ageband, subtype)],
             ageband ~ subtype, value.var = "pos_pct"),
       file.path(OUT_DIR, "covid_positivity_by_ageband.csv"))
