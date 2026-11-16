# PWAS engine, called by 02_run_pwas.R: Firth logistic regression of a binary outcome on
# each protein (complete cases per protein), BH-FDR and Bonferroni across all tested proteins.
library(optparse)
library(data.table)
library(logistf)
library(parallel)

opt <- parse_args(OptionParser(option_list = list(
  make_option("--olink",       type = "character", help = "QC'd protein matrix (eid + one NPX column per protein)"),
  make_option("--phenotype",   type = "character", help = "comma-separated phenotype files, inner-joined on eid"),
  make_option("--outcome",     type = "character", help = "binary outcome (0/1)"),
  make_option("--covariates",  type = "character", help = "comma-separated"),
  make_option("--factor-vars", type = "character", help = "comma-separated"),
  make_option("--output",      type = "character"),
  make_option("--n-jobs",      type = "integer", default = 1)
)))
covariates  <- trimws(strsplit(opt$covariates, ",")[[1]])
factor_vars <- trimws(strsplit(opt$`factor-vars`, ",")[[1]])

olink    <- fread(opt$olink)
proteins <- grep("eid|batch|plate|well|qc|index", names(olink),
                 ignore.case = TRUE, invert = TRUE, value = TRUE)   # drop ID and technical columns
pheno <- Reduce(function(x, y) merge(x, y, by = "eid"),
                lapply(trimws(strsplit(opt$phenotype, ",")[[1]]), fread))

d <- merge(olink, pheno[, c("eid", opt$outcome, covariates), with = FALSE], by = "eid")
set(d, j = opt$outcome, value = as.numeric(d[[opt$outcome]]))
d <- d[d[[opt$outcome]] %in% c(0, 1)]
for (v in covariates)
  if (v %in% factor_vars || inherits(d[[v]], "integer64") || is.character(d[[v]]) || is.factor(d[[v]]))
    set(d, j = v, value = as.factor(as.character(d[[v]])))
proteins <- intersect(proteins, names(d))

run_firth_single <- function(protein, data, outcome, covariates) {
  ad <- na.omit(data[, c(protein, outcome, covariates), with = FALSE])
  n_total <- nrow(ad)
  n_case  <- sum(ad[[outcome]])
  failed <- function(msg) data.table(
    protein = protein, n_total = n_total, n_case = n_case, beta = NA_real_, se = NA_real_,
    or = NA_real_, or_lower = NA_real_, or_upper = NA_real_, z = NA_real_, pval = NA_real_,
    converged = FALSE, error = msg)
  if (n_total < 50 || n_case < 5) return(failed("Insufficient samples or cases"))
  fml <- as.formula(paste0(outcome, " ~ ", protein, " + ", paste(covariates, collapse = " + ")))
  tryCatch({
    fit  <- logistf(fml, data = ad, control = logistf.control(maxit = 250))
    idx  <- 2L   # the protein term (first after the intercept)
    beta <- fit$coefficients[idx]
    se   <- sqrt(diag(vcov(fit)))[idx]
    data.table(
      protein = protein, n_total = n_total, n_case = n_case, beta = beta, se = se,
      or = exp(beta),
      or_lower = exp(fit$ci.lower[idx]),   # profile-likelihood CI
      or_upper = exp(fit$ci.upper[idx]),
      z = beta / se, pval = fit$prob[idx], converged = TRUE, error = NA_character_)
  }, error = function(e) failed(as.character(e$message)))
}

cl <- makeCluster(opt$`n-jobs`)
clusterExport(cl, c("d", "opt", "covariates", "run_firth_single"))
clusterEvalQ(cl, { library(data.table); library(logistf) })
results <- rbindlist(parLapply(cl, proteins, function(p) run_firth_single(p, d, opt$outcome, covariates)))
stopCluster(cl)

# One multiple-testing family: every protein with a p-value
ok <- !is.na(results$pval)
results[, pval_fdr := NA_real_]
results[ok, pval_fdr := p.adjust(pval, method = "BH")]
results[, pval_bonf := pmin(pval * sum(ok), 1)]
results[, `:=`(sig_nominal = pval < 0.05, sig_suggestive = pval < 0.001,
               sig_fdr = pval_fdr < 0.05, sig_bonf = pval_bonf < 0.05)]
setorder(results, pval)
fwrite(results, opt$output)
