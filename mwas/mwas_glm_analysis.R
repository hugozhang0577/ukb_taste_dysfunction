# Logistic regression of the outcome on each NMR measure plus covariates (glm, Wald 95% CI),
# one model per measure across parallel workers; BH-FDR and Bonferroni over the measures tested.
library(data.table)
library(parallel)

fit_one <- function(feature, data, outcome, covariates) {
  d <- na.omit(data[, c(feature, outcome, covariates), with = FALSE])
  res <- data.table(protein = feature, n_total = nrow(d), n_case = sum(d[[outcome]]),
                    beta = NA_real_, se = NA_real_, or = NA_real_, or_lower = NA_real_,
                    or_upper = NA_real_, z = NA_real_, pval = NA_real_, converged = FALSE,
                    error = "Insufficient samples or cases")
  if (res$n_total < 50 || res$n_case < 5) return(res)
  tryCatch({
    fit <- glm(as.formula(paste0(outcome, " ~ ", feature, " + ", paste(covariates, collapse = " + "))),
               data = d, family = binomial(link = "logit"))
    cf  <- summary(fit)$coefficients
    idx <- which(rownames(cf) == feature)
    if (length(idx) == 0) idx <- 2
    ci  <- confint.default(fit)[feature, ]
    res[, `:=`(beta = cf[idx, "Estimate"], se = cf[idx, "Std. Error"],
               or = exp(cf[idx, "Estimate"]), or_lower = exp(ci[1]), or_upper = exp(ci[2]),
               z = cf[idx, "z value"], pval = cf[idx, "Pr(>|z|)"],
               converged = fit$converged, error = NA_character_)]
  }, error = function(e) res[, error := conditionMessage(e)])
}

mwas_glm <- function(nmr_file, pheno_files, outcome, covariates, factor_vars, output_file, n_jobs) {
  nmr      <- fread(nmr_file)
  features <- setdiff(names(nmr), grep("eid|batch|plate|well|qc|index", names(nmr),
                                       ignore.case = TRUE, value = TRUE))
  pheno <- Reduce(function(x, y) merge(x, y, by = "eid"), lapply(pheno_files, fread))
  dt    <- merge(nmr, pheno[, c("eid", outcome, covariates), with = FALSE], by = "eid")
  set(dt, j = outcome, value = as.numeric(dt[[outcome]]))
  dt <- dt[dt[[outcome]] %in% c(0, 1)]
  # Listed factor variables, and any text column, enter the model as factors
  for (v in covariates) {
    x <- dt[[v]]
    if (v %in% factor_vars || inherits(x, "integer64")) x <- as.character(x)
    if (is.character(x) || is.factor(x)) set(dt, j = v, value = as.factor(x))
  }
  features <- features[features %in% names(dt)]

  cl <- makeCluster(n_jobs)
  clusterEvalQ(cl, library(data.table))
  res <- rbindlist(parLapply(cl, features, fit_one, data = dt, outcome = outcome,
                             covariates = covariates))
  stopCluster(cl)

  res[!is.na(pval), pval_fdr := p.adjust(pval, method = "BH")]
  res[, pval_bonf := pmin(pval * sum(!is.na(pval)), 1)]
  res[, `:=`(sig_nominal = pval < 0.05, sig_suggestive = pval < 0.001,
             sig_fdr = pval_fdr < 0.05, sig_bonf = pval_bonf < 0.05)]
  setorder(res, pval)
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  fwrite(res, output_file)
}
