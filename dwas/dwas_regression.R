# DWAS engine, called by 04_run_dwas.R: logistic regression of the outcome on each ICD-10
# code (Firth penalized likelihood if glm warns), BH-FDR and Bonferroni across all tested codes.
library(data.table)
library(optparse)
library(logistf)
library(parallel)

opt <- parse_args(OptionParser(option_list = list(
  make_option("--matrix",      type = "character", help = "ICD-10 exposure matrix (eid + one logical column per code)"),
  make_option("--phenotype",   type = "character"),
  make_option("--dx-count",    type = "character", help = "eid, dx_count_pre"),
  make_option("--dictionary",  type = "character"),
  make_option("--outcome",     type = "character"),
  make_option("--covariates",  type = "character", help = "comma-separated"),
  make_option("--factor-vars", type = "character", help = "comma-separated"),
  make_option("--min-exposed", type = "integer"),
  make_option("--output",      type = "character"),
  make_option("--n-jobs",      type = "integer", default = 1)
)))
covariates  <- strsplit(opt$covariates, ",")[[1]]
factor_vars <- strsplit(opt$`factor-vars`, ",")[[1]]

mat   <- readRDS(opt$matrix)
pheno <- fread(opt$phenotype)
if ("dx_count_pre" %in% names(pheno)) pheno[, dx_count_pre := NULL]
pheno <- merge(pheno, fread(opt$`dx-count`, select = c("eid", "dx_count_pre")), by = "eid", all.x = TRUE)
pheno <- pheno[eid %in% mat$eid & !is.na(get(opt$outcome))]

# Codes tested: exposed count among participants with a non-missing outcome >= min-exposed
all_codes  <- setdiff(names(mat), "eid")
idx_all    <- match(pheno$eid, mat$eid)
n_exp_all  <- vapply(all_codes, function(cd) sum(mat[[cd]][idx_all]), integer(1))
analyzable <- all_codes[n_exp_all >= opt$`min-exposed`]

# Complete cases on all model terms; the exposure has no missing values,
# so the analysis set is the same for every code
d <- pheno[, c("eid", opt$outcome, covariates), with = FALSE]
d <- d[complete.cases(d)]
for (fv in intersect(factor_vars, covariates)) set(d, j = fv, value = as.factor(d[[fv]]))
idx  <- match(d$eid, mat$eid)
rows <- setNames(lapply(analyzable, function(cd) which(mat[[cd]][idx])), analyzable)
d <- as.data.frame(d[, eid := NULL])

formula_base <- as.formula(sprintf("%s ~ x + %s", opt$outcome, paste(covariates, collapse = " + ")))

fit_code <- function(cd) {
  dd <- d; dd$x <- 0L; dd$x[rows[[cd]]] <- 1L
  base <- data.table(icd10 = cd, n_total = nrow(dd), n_case_taste = sum(dd[[opt$outcome]] == 1),
                     n_exposed = length(rows[[cd]]), n_unexposed = nrow(dd) - length(rows[[cd]]))
  failed <- function(msg) cbind(base, data.table(
    beta = NA_real_, se = NA_real_, z = NA_real_, pval = NA_real_,
    or_lower = NA_real_, or_upper = NA_real_, converged = FALSE, method = msg))
  tryCatch({
    co <- summary(glm(formula_base, data = dd, family = binomial()))$coefficients["x", ]
    cbind(base, data.table(
      beta = co[["Estimate"]], se = co[["Std. Error"]], z = co[["z value"]], pval = co[["Pr(>|z|)"]],
      or_lower = exp(co[["Estimate"]] - 1.96 * co[["Std. Error"]]),
      or_upper = exp(co[["Estimate"]] + 1.96 * co[["Std. Error"]]),
      converged = TRUE, method = "glm"))
  }, warning = function(w) {
    # any glm warning (non-convergence, separation) -> Firth logistic regression
    tryCatch({
      ff <- logistf(formula_base, data = dd)
      i  <- which(names(ff$coefficients) == "x")
      cbind(base, data.table(
        beta = ff$coefficients[[i]], se = sqrt(diag(vcov(ff)))[[i]], z = NA_real_,
        pval = ff$prob[[i]], or_lower = exp(ff$ci.lower[[i]]), or_upper = exp(ff$ci.upper[[i]]),
        converged = TRUE, method = "firth"))
    }, error = function(e) failed(paste("failed:", conditionMessage(e))))
  }, error = function(e) failed(paste("error:", conditionMessage(e))))
}

results <- rbindlist(mclapply(analyzable, fit_code, mc.cores = opt$`n-jobs`), fill = TRUE)
results[, `:=`(or = exp(beta), n_exposed_outcome = n_exp_all[icd10])]

# One multiple-testing family: every tested code
results[, pval_fdr  := p.adjust(pval, method = "fdr")]
results[, pval_bonf := p.adjust(pval, method = "bonferroni")]

results <- merge(results, fread(opt$dictionary, select = c("icd10", "description", "chapter", "domain", "oral")),
                 by = "icd10", all.x = TRUE)
setcolorder(results, c("icd10", "description", "chapter", "domain", "oral",
  "beta", "se", "z", "pval", "or", "or_lower", "or_upper", "n_total", "n_case_taste",
  "n_exposed", "n_unexposed", "n_exposed_outcome", "converged", "method",
  "pval_fdr", "pval_bonf"))
setorder(results, pval)
fwrite(results, opt$output)
