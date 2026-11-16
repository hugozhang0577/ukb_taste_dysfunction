# ExWAS regression engine: one logistic regression of the outcome on each exposure plus
# covariates (glm, with Firth fallback), then Benjamini-Hochberg FDR within each exposure family.
library(data.table)
library(logistf)

fit_one <- function(sub, v, covars, var_type, method, cell_firth) {
  fml <- as.formula(paste(".y ~", paste(c(v, covars), collapse = " + ")))
  out <- list(beta = NA_real_, se = NA_real_, OR = NA_real_, OR_lower = NA_real_,
              OR_upper = NA_real_, p_value = NA_real_, method_used = NA_character_,
              converged = FALSE, error_msg = NA_character_)
  force_firth <- method == "firth"
  # Hybrid mode: binary exposures with a 2x2 cell below cell_firth go straight to Firth
  if (method == "hybrid" && var_type %in% c("binary", "derived_binary")) {
    tab <- table(sub$.y, sub[[v]])
    if (nrow(tab) < 2 || ncol(tab) < 2) { out$error_msg <- "degenerate_2x2"; return(out) }
    if (min(tab) < cell_firth) force_firth <- TRUE
  }

  sep_warn <- FALSE; glm_err <- NULL; fit <- NULL
  if (!force_firth) {
    fit <- tryCatch(
      withCallingHandlers(
        glm(fml, data = sub, family = binomial()),
        warning = function(w) {
          if (grepl("fitted probabilities numerically 0 or 1", conditionMessage(w), fixed = TRUE))
            sep_warn <<- TRUE
          invokeRestart("muffleWarning")
        }),
      error = function(e) { glm_err <<- conditionMessage(e); NULL })
  }
  # Firth when forced, when glm fails or does not converge, or on separation
  need_firth <- method != "glm" && (force_firth || is.null(fit) || !isTRUE(fit$converged) || sep_warn)

  if (!need_firth && !is.null(fit)) {
    co <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
    if (is.null(co) || !(v %in% rownames(co))) { out$error_msg <- "glm_coef_missing"; return(out) }
    b <- co[v, "Estimate"]; se <- co[v, "Std. Error"]
    out$beta <- b; out$se <- se; out$p_value <- co[v, "Pr(>|z|)"]
    out$OR <- exp(b); out$OR_lower <- exp(b - 1.96 * se); out$OR_upper <- exp(b + 1.96 * se)
    out$method_used <- "glm"; out$converged <- TRUE
    return(out)
  }
  if (method == "glm") {
    out$error_msg <- if (is.null(fit)) paste0("glm_failed:", glm_err) else "glm_nonconverged"
    out$method_used <- "glm_failed"
    return(out)
  }

  fit2 <- tryCatch(logistf(fml, data = sub, control = logistf.control(maxit = 100)),
                   error = function(e) { out$error_msg <<- paste0("firth_error:", conditionMessage(e)); NULL })
  if (is.null(fit2)) { out$method_used <- "firth_failed"; return(out) }
  idx <- which(names(fit2$coefficients) == v)
  if (length(idx) == 0) { out$error_msg <- "firth_coef_missing"; return(out) }
  b <- unname(fit2$coefficients[idx])
  out$beta <- b; out$se <- unname(sqrt(diag(vcov(fit2)))[idx]); out$p_value <- unname(fit2$prob[idx])
  out$OR <- exp(b)
  out$OR_lower <- exp(unname(fit2$ci.lower[idx]))   # profile-likelihood CI
  out$OR_upper <- exp(unname(fit2$ci.upper[idx]))
  out$method_used <- if (method == "firth") "firth" else "firth_fallback"
  out$converged <- TRUE
  out
}

# omics: exposure matrix CSV (eid + one column per exposure); phenotype: outcome + covariates CSV;
# var_meta: CSV with variable, var_type, fdr_group; covariates / factor_vars: comma-separated names.
run_exwas <- function(omics, phenotype, outcome, covariates, factor_vars, var_meta, output,
                      method = "hybrid", min_cases = 100, cell_drop = 10, cell_firth = 20,
                      eid_col = "eid") {
  covars  <- strsplit(covariates, ",", fixed = TRUE)[[1]]
  factors <- strsplit(factor_vars, ",", fixed = TRUE)[[1]]
  pheno <- fread(phenotype)
  meta  <- fread(var_meta)
  # Exposures are taken from the exposure matrix, not from same-named phenotype columns
  pheno[, setdiff(intersect(names(pheno), meta$variable), eid_col) := NULL]
  dat <- merge(pheno, fread(omics), by = eid_col)
  for (cv in intersect(factors, names(dat))) dat[[cv]] <- factor(dat[[cv]])
  dat[[".y"]] <- as.integer(dat[[outcome]])
  meta <- meta[variable %in% names(dat)]

  res <- lapply(seq_len(nrow(meta)), function(i) {
    v <- meta$variable[i]; vt <- meta$var_type[i]
    sub <- dat[, c(".y", v, covars), with = FALSE]
    sub <- sub[complete.cases(sub)]
    n_cases <- sum(sub$.y == 1)
    if (n_cases < min_cases) return(NULL)            # too few cases with complete data
    vr <- var(as.numeric(sub[[v]]))
    if (is.na(vr) || vr <= 0) return(NULL)           # no variation
    if (vt %in% c("binary", "derived_binary")) {     # sparse 2x2 cell
      tab <- table(sub$.y, sub[[v]])
      if (nrow(tab) < 2 || ncol(tab) < 2 || min(tab) < cell_drop) return(NULL)
    }
    r <- fit_one(sub, v, covars, vt, method, cell_firth)
    data.table(variable = v, var_type = vt, fdr_group = meta$fdr_group[i],
               n_complete = nrow(sub), n_cases_complete = n_cases,
               beta = r$beta, se = r$se, OR = r$OR, OR_lower = r$OR_lower, OR_upper = r$OR_upper,
               p_value = r$p_value, method_used = r$method_used, converged = r$converged,
               error_msg = r$error_msg)
  })
  results <- rbindlist(res, fill = TRUE)

  # BH FDR within each exposure family
  results[, p_adj := NA_real_]
  for (g in unique(results$fdr_group))
    results[fdr_group == g & !is.na(p_value), p_adj := p.adjust(p_value, method = "BH")]
  results[, significant := !is.na(p_adj) & p_adj < 0.05]

  dir.create(dirname(output), showWarnings = FALSE, recursive = TRUE)
  fwrite(results, output)
  invisible(results)
}
