#!/usr/bin/env Rscript
# =============================================================================
# 05_covid_association.R  (SARS-CoV-2 positivity by subtype)
# =============================================================================
#
# Pre-assessment SARS-CoV-2 positivity, one subtype versus the rest, crude and
# age-adjusted. This is the quantity Figure 5E draws on the left.
#
# The point of the script is the contrast between the two models. The subtypes
# differ in age by roughly thirteen years by construction — the young-idiopathic
# subtype D is the youngest — and testing and infection both track age heavily
# over this period. So the crude comparison is confounded by age, and the
# adjusted one is the interpretable estimate:
#
#   sex_only      positivity ~ subtype + sex                   (crude)
#   age_sex       positivity ~ subtype + age + sex             (reported)
#   ns_age3_sex   positivity ~ subtype + ns(age, 3) + sex      (no linearity assumption)
#
# All three are written. `age_sex` is the model the manuscript reports; the
# spline model is there so that the adjustment cannot be dismissed as an artefact
# of forcing age in linearly.
#
# age here is age at the taste assessment, not age at recruitment: the exposure
# and the outcome are both anchored to the questionnaire.
#
# Engine is plain logistic regression. Every cell is large (roughly 500-1,200
# tested per subtype) and there is no separation, so Firth is not needed.
#
# Reading it: a subtype's positivity is a statement about documented infection in
# a tested subgroup, not about infection in the cohort. Testing was not random.
#
# CODE_DIR (env, default = current dir) must hold _subtype_map.R, captured before
# setwd(PROJECT_DIR).
#
# Input:  output/subtyping/clusters/cluster_assignments_g1_{m,f}_k4.rds
#         output/ml_ready/group1_full.rds                (age at assessment)
#         input/analysis_ready/covid_temporal_flags.csv  (eid-keyed timing flags)
# Output: output/subtyping/evidence_covid/
#           covid_onevsrest_OR.csv          -> Figure 5E, Methods, Supplementary
#           covid_subtype_labels.csv        -> Figure 5E timing join
#           covid_positivity_by_ageband.csv -> model-free age check
# =============================================================================

CODE_DIR <- Sys.getenv("CODE_DIR", unset = getwd())
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
suppressPackageStartupMessages({ library(data.table); library(splines) })

source(file.path(CODE_DIR, "_subtype_map.R"))   # SUBTYPE_MAP

CL      <- "output/subtyping/clusters"
OUT_DIR <- "output/subtyping/evidence_covid"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
SUBTYPE_NAME <- c(A = "Aging frailty", B = "Psychosomatic",
                  C = "Cardiometabolic", D = "Young idiopathic")

# ---- [1] subtype labels, derived here rather than read in -------------------
# Taking the labels from the cluster assignments each run is deliberate: a
# pre-labelled participant file is the one way this analysis can silently go
# stale after the model is refitted.
labels_of <- function(sx) {
  f <- file.path(CL, sprintf("cluster_assignments_g1_%s_k4.rds", sx))
  if (!file.exists(f)) stop("cluster assignments not found: ", f)
  fl <- readRDS(f)$final_labels
  data.table(eid     = as.integer(names(fl)),
             sex     = toupper(sx),
             subtype = unname(SUBTYPE_MAP[[sx]][as.character(as.integer(fl))]))
}
sub <- rbind(labels_of("m"), labels_of("f"))
if (anyNA(sub$subtype)) stop("unmapped cluster id -- check _subtype_map.R")
cat("subtype labels:", nrow(sub), "cases\n")
print(sub[, .N, by = subtype][order(subtype)])

# ---- [2] join COVID timing and age ------------------------------------------
flags <- fread("input/analysis_ready/covid_temporal_flags.csv")
need <- c("eid", "phe_pos_before_taste", "phe_ever_positive")
miss <- setdiff(need, names(flags))
if (length(miss)) stop("covid_temporal_flags.csv is missing: ", paste(miss, collapse = ", "))
if ("subtype" %in% names(flags)) flags[, subtype := NULL]   # relabelled below, never trusted

age <- as.data.table(readRDS("output/ml_ready/group1_full.rds"))[, .(eid, age)]

n0 <- nrow(sub)
d  <- merge(sub, flags, by = "eid")
cat("after COVID-flag join:", n0, "->", nrow(d), "\n")
n1 <- nrow(d)
d  <- merge(d, age, by = "eid")
cat("after age join:       ", n1, "->", nrow(d), "\n")

cat("\n=== age by subtype ===\n")
print(d[, .(N = .N, age_mean = round(mean(age, na.rm = TRUE), 1),
            age_sd = round(sd(age, na.rm = TRUE), 1)), by = subtype][order(subtype)])

# Denominator is the tested subgroup: a participant with no test result cannot
# contribute to a positivity comparison.
t <- d[!is.na(phe_ever_positive) & !is.na(age) & !is.na(phe_pos_before_taste)]
t[, sex := factor(sex)]
cat("\ntested with age:", nrow(t), "of", nrow(d), "\n")
if (nrow(t) == 0) stop("no tested participants after the joins")

cat("\n=== raw pre-assessment positivity by subtype ===\n")
print(t[, .(N = .N, pos_pct = round(100 * mean(phe_pos_before_taste), 1)),
        by = subtype][order(subtype)])

# ---- [3] age is a strong driver of positivity on its own --------------------
m0 <- glm(phe_pos_before_taste ~ age + sex, data = t, family = binomial())
co <- summary(m0)$coefficients["age", ]
cat(sprintf("\nage main effect: OR per +10 yr = %.3f, P = %.2e\n", exp(co[1] * 10), co[4]))

# ---- [4] one-versus-rest, three models --------------------------------------
ovr <- function(letter, rhs) {
  t2 <- copy(t)[, grp := as.integer(subtype == letter)]
  m  <- glm(as.formula(paste("phe_pos_before_taste ~ grp +", rhs)),
            data = t2, family = binomial())
  s  <- summary(m)$coefficients["grp", ]
  ci <- exp(s[1] + c(-1, 1) * 1.96 * s[2])
  data.table(subtype = letter, OR = round(exp(s[1]), 3),
             lo = round(ci[1], 3), hi = round(ci[2], 3), P = signif(s[4], 3))
}
MODELS <- c(sex_only = "sex", age_sex = "age + sex", ns_age3_sex = "ns(age, 3) + sex")
res <- rbindlist(lapply(names(MODELS), function(mn)
  rbindlist(lapply(c("A", "B", "C", "D"), ovr, rhs = MODELS[[mn]]))[, model := mn]))

nsum <- t[, .(n_tested = .N, n_pos_before = sum(phe_pos_before_taste),
              raw_pos_pct = round(100 * mean(phe_pos_before_taste), 1)), by = subtype]
res <- merge(res, nsum, by = "subtype")
res[, subtype_name := SUBTYPE_NAME[subtype]]
setcolorder(res, c("model", "subtype", "subtype_name", "OR", "lo", "hi", "P",
                   "n_tested", "n_pos_before", "raw_pos_pct"))
res <- res[order(factor(model, levels = names(MODELS)), subtype)]
cat("\n=== one-versus-rest odds ratios ===\n"); print(res)

fwrite(res, file.path(OUT_DIR, "covid_onevsrest_OR.csv"))
fwrite(d[, .(eid, sex, subtype, phe_pos_before_taste, phe_ever_positive)],
       file.path(OUT_DIR, "covid_subtype_labels.csv"))
# Export and dx upload to RAP  (Figure 5E and the Methods read the OR table)

# ---- [5] the same contrast without a model ----------------------------------
t[, ageband := cut(age, c(0, 60, 65, 70, 75, 200),
                   labels = c("<60", "60-65", "65-70", "70-75", "75+"), right = FALSE)]
band <- dcast(t[, .(pos_pct = round(100 * mean(phe_pos_before_taste), 0)),
                by = .(ageband, subtype)], ageband ~ subtype, value.var = "pos_pct")
cat("\n=== raw positivity (%) within age band ===\n"); print(band)
cat("N per cell:\n"); print(dcast(t[, .N, by = .(ageband, subtype)],
                                  ageband ~ subtype, value.var = "N"))
fwrite(band, file.path(OUT_DIR, "covid_positivity_by_ageband.csv"))

cat("\ndone ->", OUT_DIR, "\n")
