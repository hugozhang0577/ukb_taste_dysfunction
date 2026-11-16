# Smell change by subtype (group-1 cases): co-occurrence (a clustering input, descriptive only), and duration and
# functional impact among cases with smell change (not used in clustering), tested across subtypes without
# adjustment (chi-squared / Kruskal-Wallis with effect sizes) and adjusted for age and sex.
library(data.table)
source("unsupervised_subtyping/_subtype_map.R")
setwd(Sys.getenv("PROJECT_DIR", "."))
OUT <- "output/subtyping/evidence_severity"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

labels_of <- function(sx) {
  fl <- readRDS(sprintf("output/subtyping/clusters/cluster_assignments_g1_%s_k4.rds", sx))$final_labels
  data.table(eid = as.integer(names(fl)), sex = toupper(sx), subtype = unname(SUBTYPE_MAP[[sx]][as.character(as.integer(fl))]))
}
d <- merge(rbind(labels_of("m"), labels_of("f")),
           as.data.table(readRDS("output/ml_ready/group1_full.rds"))[, .(eid, age, smell_any, smell_time, smell_extent)],
           by = "eid")
cramers_v <- function(ch, tab) sqrt(ch$statistic / (sum(tab) * (min(dim(tab)) - 1)))

# co-occurrence of smell change
ct <- table(d$subtype, d$smell_any); chi <- chisq.test(ct); cramV <- cramers_v(chi, ct)

# duration band (0-3) among cases with smell change
ds <- d[smell_any == 1 & !is.na(smell_time)]
kw <- kruskal.test(smell_time ~ factor(subtype), data = ds)
eps2 <- (kw$statistic - 4 + 1) / (nrow(ds) - 4)   # epsilon-squared, 4 groups
fit_t <- lm(smell_time ~ factor(subtype) + age + factor(sex), data = ds)
p_adj_t <- anova(fit_t)["factor(subtype)", "Pr(>F)"]

# functional impact (0/1) among cases with smell change
de <- d[smell_any == 1 & !is.na(smell_extent)]
cte <- table(de$subtype, de$smell_extent); chie <- chisq.test(cte); cramVe <- cramers_v(chie, cte)
fit_e <- glm(smell_extent ~ factor(subtype) + age + factor(sex), data = de, family = binomial())
p_adj_e <- anova(fit_e, test = "Chisq")["factor(subtype)", "Pr(>Chi)"]

summ <- rbind(
  d[, .(measure = "smell_cooccur(any)", n = .N, val = round(100 * mean(smell_any == 1, na.rm = TRUE), 1)), by = subtype],
  ds[, .(measure = "smell_time_mean", n = .N, val = round(mean(smell_time), 2)), by = subtype],
  de[, .(measure = "smell_extent_pct", n = .N, val = round(100 * mean(smell_extent == 1), 1)), by = subtype])
fwrite(dcast(summ, subtype ~ measure, value.var = "val"), file.path(OUT, "severity_by_subtype.csv"))

fwrite(data.table(
  descriptor = c("smell_cooccur(any)", "smell_time_mean", "smell_extent_pct"),
  measure_label = c("Smell co-occurrence (%)", "Smell duration, mean band 0-3", "Smell functional impact (%)"),
  mofa_input = c("yes (descriptive only)", "no (excluded)", "no (excluded)"),
  test = c("Pearson chi-squared (subtype x smell_any, 4x2)", "Kruskal-Wallis (smell_time by subtype)",
           "Pearson chi-squared (subtype x smell_extent, 4x2)"),
  statistic = c(unname(chi$statistic), unname(kw$statistic), unname(chie$statistic)),
  df = c(unname(chi$parameter), unname(kw$parameter), unname(chie$parameter)),
  P = c(chi$p.value, kw$p.value, chie$p.value),
  effect_size_name = c("Cramer's V", "epsilon-squared", "Cramer's V"),
  effect_size_value = c(unname(cramV), unname(eps2), unname(cramVe)),
  n = c(sum(ct), nrow(ds), sum(cte)),
  P_adj_age_sex = c(NA_real_, p_adj_t, p_adj_e),
  adj_model = c(NA_character_, "lm(smell_time ~ subtype + age + sex), anova F test",
                "glm(smell_extent ~ subtype + age + sex, binomial), LR test"),
  n_adj = c(NA_integer_, nobs(fit_t), nobs(fit_e))),
  file.path(OUT, "severity_omnibus.csv"))
