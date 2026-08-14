#!/usr/bin/env Rscript
# =============================================================================
# 03_severity_orthogonality.R  (severity-orthogonality evidence)
#
# Are the four subtypes orthogonal to chemosensory severity, or just severity
# gradations of one axis? Severity measures (excluded from MOFA clustering, so
# non-circular):
#   smell_time   (ordinal 0-3, smell-change duration band)
#   smell_extent (binary, smell functional impact)
# Breadth proxy (smell_any, taste+smell co-occurrence) is a MOFA input ->
# descriptive only. Taste severity is not gradable within cases (taste_time/
# extent are not retained), so smell severity is the available proxy.
#
# Test across A/B/C/D: omnibus + effect size + age/sex-adjusted; check for a
# monotonic gradient and within-subtype spread. Orthogonal => small effect, no
# A<B<C<D ladder, severity varies within subtype.
#
# CODE_DIR (env, default = current dir) must hold _subtype_map.R.
# =============================================================================

CODE_DIR <- Sys.getenv("CODE_DIR", unset = getwd())
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
suppressPackageStartupMessages(library(data.table))
source(file.path(CODE_DIR, "_subtype_map.R"))

OUT <- "output/subtyping/evidence_severity"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

gl <- function(p, sx) {
  fl <- readRDS(p)$final_labels
  data.table(eid = as.integer(names(fl)), sex = sx,
             subtype = unname(SUBTYPE_MAP[[tolower(sx)]][as.character(as.integer(fl))]))
}
subs <- rbind(gl("output/subtyping/clusters/cluster_assignments_g1_m_k4.rds", "M"),
              gl("output/subtyping/clusters/cluster_assignments_g1_f_k4.rds", "F"))
ml <- as.data.table(readRDS("output/ml_ready/group1_full.rds"))
d <- merge(subs, ml[, .(eid, age, sex_ml = sex, smell_any, smell_time, smell_extent)], by = "eid")
cat("[severity] subtype cases:", nrow(d), "| with smell change:", d[smell_any == 1, .N], "\n")

# ---- breadth: smell co-occurrence (smell_any; MOFA input -> descriptive) -----
cat("=== breadth: smell co-occurrence % by subtype (descriptive) ===\n")
print(d[, .(n = .N, smell_cooccur_pct = round(100 * mean(smell_any == 1, na.rm = TRUE), 1)), by = subtype][order(subtype)])
ct <- table(d$subtype, d$smell_any); chi <- chisq.test(ct)
cramV <- sqrt(chi$statistic / (sum(ct) * (min(dim(ct)) - 1)))
cat(sprintf("  chi2 P=%.2e | Cramer's V=%.3f\n", chi$p.value, cramV))

# ---- clean severity measures (excluded from clustering) ----------------------
ds <- d[smell_any == 1 & !is.na(smell_time)]
cat("=== smell severity subset N:", nrow(ds), "===\n")
cat("\n--- smell_time (duration band 0-3) by subtype ---\n")
print(ds[, .(n = .N, median = as.numeric(median(smell_time)), mean = round(mean(smell_time), 2),
             pct_band3 = round(100 * mean(smell_time == 3), 1)), by = subtype][order(subtype)])
kw <- kruskal.test(smell_time ~ factor(subtype), data = ds)
eps2 <- (kw$statistic - 4 + 1) / (nrow(ds) - 4)
cat(sprintf("  Kruskal-Wallis P=%.2e | epsilon^2=%.4f (<0.01 negligible)\n", kw$p.value, eps2))
fit_t <- lm(smell_time ~ factor(subtype) + age + factor(sex), data = ds)
p_adj_t <- anova(fit_t)["factor(subtype)", "Pr(>F)"]
cat("  age+sex-adjusted subtype effect (lm) anova P =", signif(p_adj_t, 3),
    "| model N =", nobs(fit_t), "\n")

cat("\n--- smell_extent (functional impact 0/1) by subtype ---\n")
de <- d[smell_any == 1 & !is.na(smell_extent)]
print(de[, .(n = .N, impact_pct = round(100 * mean(smell_extent == 1), 1)), by = subtype][order(subtype)])
cte <- table(de$subtype, de$smell_extent); chie <- chisq.test(cte)
cramVe <- sqrt(chie$statistic / (sum(cte) * (min(dim(cte)) - 1)))
cat(sprintf("  chi2 P=%.2e | Cramer's V=%.3f\n", chie$p.value, cramVe))
fit_e <- glm(smell_extent ~ factor(subtype) + age + factor(sex), data = de, family = binomial())
p_adj_e <- anova(fit_e, test = "Chisq")["factor(subtype)", "Pr(>Chi)"]
cat("  age+sex-adjusted subtype LR test P =", signif(p_adj_e, 3),
    "| model N =", nobs(fit_e), "\n")

# ---- save summary -----------------------------------------------------------
summ <- rbind(
  d[, .(measure = "smell_cooccur(any)", n = .N, val = round(100 * mean(smell_any == 1, na.rm = TRUE), 1)), by = subtype],
  ds[, .(measure = "smell_time_mean", n = .N, val = round(mean(smell_time), 2)), by = subtype],
  de[, .(measure = "smell_extent_pct", n = .N, val = round(100 * mean(smell_extent == 1), 1)), by = subtype])
fwrite(dcast(summ, subtype ~ measure, value.var = "val"), file.path(OUT, "severity_by_subtype.csv"))
cat("\n[written]", file.path(OUT, "severity_by_subtype.csv"), "\n")

# ---- save omnibus tests + effect sizes --------------------------------------
# Machine-readable companion so downstream table builders never hard-code
# P values / effect sizes.
omni <- data.table(
  descriptor = c("smell_cooccur(any)", "smell_time_mean", "smell_extent_pct"),
  measure_label = c("Smell co-occurrence (%)",
                    "Smell duration, mean band 0-3",
                    "Smell functional impact (%)"),
  mofa_input = c("yes (descriptive only)", "no (excluded)", "no (excluded)"),
  test = c("Pearson chi-squared (subtype x smell_any, 4x2)",
           "Kruskal-Wallis (smell_time by subtype)",
           "Pearson chi-squared (subtype x smell_extent, 4x2)"),
  statistic = c(unname(chi$statistic), unname(kw$statistic), unname(chie$statistic)),
  df = c(unname(chi$parameter), unname(kw$parameter), unname(chie$parameter)),
  P = c(chi$p.value, kw$p.value, chie$p.value),
  effect_size_name = c("Cramer's V", "epsilon-squared", "Cramer's V"),
  effect_size_value = c(unname(cramV), unname(eps2), unname(cramVe)),
  n = c(sum(ct), nrow(ds), sum(cte)),
  P_adj_age_sex = c(NA_real_, p_adj_t, p_adj_e),
  adj_model = c(NA_character_,
                "lm(smell_time ~ subtype + age + sex), anova F test",
                "glm(smell_extent ~ subtype + age + sex, binomial), LR test"),
  n_adj = c(NA_integer_, nobs(fit_t), nobs(fit_e)))
fwrite(omni, file.path(OUT, "severity_omnibus.csv"))
cat("[written]", file.path(OUT, "severity_omnibus.csv"), "\n")
