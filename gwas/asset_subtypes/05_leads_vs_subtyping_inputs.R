# Are the ASSET lead variants determinants of a subtyping input? In group-1 controls, each MOFA+ input
# (standardized) is regressed on the lead dosage (per ASSET effect allele) adjusted for age at baseline and sex.
# Call per variant: circular if any input has P < 5e-8, possible if any has P < 0.05 / (inputs x variants),
# otherwise clear. APOE e4 dose (itself an input) is the positive control.
library(data.table)
setwd(Sys.getenv("PROJECT_DIR", "."))
OUT <- "output/gwas_subtypes/asset"

LOCI <- data.table(tag = c("chr2_COL4A3", "chr20_MKKS", "chr6_107Mb"), lead = c("rs13028944", "rs7269148", "rs59668897"))
LOCI <- merge(LOCI, fread(file.path(OUT, "asset_suggestive_p1e-6.csv"))[, .(lead = SNP, a_eff, asset_p2, pos_subset, neg_subset)],
              by = "lead")

vp    <- fread("output/subtyping/inputs/view_partition.csv")
COVS  <- c("age_baseline", "sex")
feats <- setdiff(vp$feature_id, c(COVS, "age", "assess_centre_id", "years_baseline_to_taste"))
ml <- as.data.table(readRDS("output/ml_ready/group1_full.rds"))
ml <- ml[taste_2w_strict == 0, c("eid", COVS, setdiff(feats, COVS)), with = FALSE]
for (f in feats) if (is.logical(ml[[f]])) set(ml, j = f, value = as.numeric(ml[[f]]))   # prior-diagnosis inputs

read_dose <- function(lead, a_eff) {
  r  <- fread(sprintf("input/gwas_genotype/lead_dosage_%s.raw", lead))
  vc <- grep(paste0("^", lead, "_"), names(r), value = TRUE)
  d  <- data.table(eid = as.integer(r$IID), g = r[[vc]])
  if (sub(".*_", "", vc) != a_eff) d[, g := 2 - g]   # count the ASSET effect allele
  d
}
test_variant <- function(dt, gname, flist) rbindlist(lapply(flist, function(f) {
  y  <- dt[[f]]
  ok <- is.finite(y) & is.finite(dt$g) & is.finite(dt$age_baseline) & is.finite(dt$sex)
  if (sum(ok) < 500 || var(y[ok]) == 0) return(NULL)
  X <- cbind(1, dt$g[ok], dt$age_baseline[ok], dt$sex[ok])
  if (var(X[, 4]) == 0) X <- X[, 1:3]   # female-only inputs
  fit <- lm.fit(X, (y[ok] - mean(y[ok])) / sd(y[ok]))
  dfr <- sum(ok) - ncol(X)
  se  <- sqrt(sum(fit$residuals^2) / dfr * chol2inv(qr.R(fit$qr))[2, 2])
  b   <- fit$coefficients[2]
  data.table(variant = gname, feature = f, view = vp[feature_id == f, view][1], n = sum(ok),
             beta_sd_per_allele = b, se = se, p = 2 * pt(-abs(b / se), dfr))
}))

res <- lapply(seq_len(nrow(LOCI)), function(i)
  test_variant(merge(ml, read_dose(LOCI$lead[i], LOCI$a_eff[i]), by = "eid"), LOCI$lead[i], feats))
res[[length(res) + 1]] <- test_variant(copy(ml)[, g := as.numeric(APOE_e4_dose)], "APOE_e4_dose", setdiff(feats, "APOE_e4_dose"))
tests <- rbindlist(res)

BONF <- 0.05 / (length(feats) * nrow(LOCI))
summ <- tests[, .(n_features = .N, min_p = min(p), top_feature = feature[which.min(p)],
                top_beta_sd = beta_sd_per_allele[which.min(p)], n_p_gw = sum(p < 5e-8), n_p_bonf = sum(p < BONF)), by = variant]
summ[, call := fifelse(n_p_gw > 0, "circular", fifelse(n_p_bonf > 0, "possible", "clear"))]
summ <- merge(summ, LOCI[, .(variant = lead, tag, asset_p2, pos_subset, neg_subset)], by = "variant", all.x = TRUE)
fwrite(tests, file.path(OUT, "lead_vs_subtyping_inputs_all.csv"))
fwrite(summ, file.path(OUT, "lead_vs_subtyping_inputs_summary.csv"))
