# Held-out evaluation of the polygenic score at each P-value threshold: logistic regression of the outcome on
# the standardized score (OR per SD, 95% CI), Nagelkerke R2 and AUC.
library(data.table)
library(pROC)
setwd(Sys.getenv("PROJECT_DIR", "."))

prs   <- fread("output/prs/prs_combined_all_scores.txt", colClasses = list(character = c("FID", "IID")))
pheno <- fread("output/prs/heldout_prsice_pheno.txt", colClasses = list(character = c("FID", "IID")))
data  <- merge(prs, pheno, by = c("FID", "IID"))[!is.na(pheno)]

rows <- rbindlist(lapply(grep("^Pt_", names(data), value = TRUE), function(col) {
  if (sd(data[[col]]) == 0) return(NULL)
  d <- na.omit(data.frame(y = data$pheno, prs_z = scale(data[[col]])[, 1]))
  fit  <- glm(y ~ prs_z, data = d, family = binomial())
  fit0 <- glm(y ~ 1, data = d, family = binomial())
  s <- summary(fit)$coefficients["prs_z", ]
  n <- nrow(d)
  r2_cox <- 1 - exp(2 * (logLik(fit0)[1] - logLik(fit)[1]) / n)
  data.table(Threshold = col, N = n, Cases = sum(d$y == 1), Controls = sum(d$y == 0),
             Beta = round(s[[1]], 4), SE = round(s[[2]], 4), OR = round(exp(s[[1]]), 3),
             CI_low = round(exp(s[[1]] - 1.96 * s[[2]]), 3), CI_high = round(exp(s[[1]] + 1.96 * s[[2]]), 3),
             Z = round(s[[3]], 3), P = s[[4]],
             R2 = round(r2_cox / (1 - exp(2 * logLik(fit0)[1] / n)), 6),   # Nagelkerke
             AUC = round(as.numeric(auc(roc(d$y, predict(fit, type = "response"), quiet = TRUE))), 3))
}))
fwrite(rows[order(P)], "output/prs/prs_heldout_performance.csv")
