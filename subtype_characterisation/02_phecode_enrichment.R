#!/usr/bin/env Rscript
# =============================================================================
# 02_phecode_enrichment.R  (disease-comorbidity evidence)
#
# Per-subtype full-PheCode-phenome enrichment. Two design red lines:
#   (1) De-circularisation: the discovery universe is the main-analysis-tested
#       PheCodes (cohort case>=100) MINUS the PheCodes that were MOFA inputs
#       (those are enriched by construction): 454 - 33 = 421 candidates.
#   (2) Adjust for age + sex (age was a MOFA input; PheCode burden tracks age).
#
# Two-layer power design (the scan is cluster-vs-rest within the ~5,772 cases):
#   TEST      : focal cell n_subtype_with >= 20 (PheWAS minimum, Denny 2010); Firth
#               (Heinze-Schemper 2002) when any 2x2 cell < 20.
#   INTERPRET : flag well_powered when n_subtype_with >= 50.
#   BH-FDR within subtype over that subtype's tested set.
# Characterisation / hypothesis-generating; interpret well-powered hits.
#
# CODE_DIR (env, default = current dir) must hold _subtype_map.R.
# =============================================================================

CODE_DIR <- Sys.getenv("CODE_DIR", unset = getwd())
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
suppressPackageStartupMessages({ library(data.table); library(logistf) })
source(file.path(CODE_DIR, "_subtype_map.R"))   # SUBTYPE_MAP$m / SUBTYPE_MAP$f

# ---- paths ------------------------------------------------------------------
F_CLUSTER_M  <- "output/subtyping/clusters/cluster_assignments_g1_m_k4.rds"
F_CLUSTER_F  <- "output/subtyping/clusters/cluster_assignments_g1_f_k4.rds"
F_VIEW_PART  <- "output/subtyping/inputs/view_partition.csv"   # MOFA-input PheCodes
F_ANALYZABLE <- "input/assoc_results/dwas_phecode_primary.csv"  # case>=100 tested PheCodes
F_ML_READY   <- "output/ml_ready/group1_full.rds"      # age (at taste assessment)
F_PHECODE    <- "input/analysis_ready/phecode_matrix_group1.rds"

OUT_DIR <- "output/subtyping/evidence_disease"
LOG_DIR <- file.path(OUT_DIR, "logs")
dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)
TEST_FLOOR <- 20L; WELLPOW <- 50L; SEED <- 20260528L

for (f in c(F_CLUSTER_M, F_CLUSTER_F, F_VIEW_PART, F_ANALYZABLE, F_ML_READY, F_PHECODE))
  if (!file.exists(f)) stop("missing input: ", f)
cat("[ok] inputs present\n")

# ---- unify subtype labels ---------------------------------------------------
get_labels <- function(path, sex) {
  fl <- readRDS(path)$final_labels; stopifnot(!is.null(fl), length(fl) > 0)
  map <- SUBTYPE_MAP[[tolower(sex)]]
  data.table(eid = as.integer(names(fl)), sex = sex, subtype = unname(map[as.character(as.integer(fl))]))
}
subs <- rbind(get_labels(F_CLUSTER_M, "M"), get_labels(F_CLUSTER_F, "F"))
stopifnot(!anyNA(subs$subtype))
cat("[subtype] total cases:", nrow(subs), "\n")
print(subs[, .N, by = .(subtype, sex)][order(subtype, sex)])

# ---- covariate (age) + PheCode matrix ---------------------------------------
cov <- as.data.table(readRDS(F_ML_READY))[, .(eid, age)]
phe <- as.data.table(readRDS(F_PHECODE)); pc_all <- setdiff(names(phe), "eid")
cat("[phe] rows:", nrow(phe), "| PheCode cols:", length(pc_all), "\n")

# ---- discovery universe = analyzable - MOFA inputs --------------------------
vp <- fread(F_VIEW_PART)
input_phe <- sub("^phe", "", vp[view == "phecode", feature_id])
stopifnot(length(intersect(input_phe, pc_all)) >= 30)   # de-circularisation must bite
# PheCodes are zero-padded strings ("079"), not numbers: without
# keepLeadingZeros the column parses as numeric and the padded codes no longer
# match the PheCode matrix column names, silently dropping them from the universe.
analyzable <- as.character(fread(F_ANALYZABLE, keepLeadingZeros = TRUE)$phecode)
universe   <- setdiff(analyzable, input_phe)
pc_test    <- intersect(universe, pc_all)
cat(sprintf("[universe] analyzable:%d  MOFA-input excluded:%d  universe:%d  in-matrix:%d\n",
            length(analyzable), length(intersect(input_phe, pc_all)), length(universe), length(pc_test)))
stopifnot(length(pc_test) > 300)

# ---- assemble (fail-visibly) ------------------------------------------------
dat0 <- merge(subs, cov, by = "eid")
cat("[merge age] N:", nrow(subs), "->", nrow(dat0), "\n")
dat <- merge(dat0, phe, by = "eid"); dat[, sex := factor(sex)]
cat("[merge phe] N:", nrow(dat0), "->", nrow(dat), "\n")

# ---- per-subtype x per-PheCode logistic -------------------------------------
set.seed(SEED); out_list <- list()
for (a in c("A", "B", "C", "D")) {
  ia_all <- as.integer(dat$subtype == a)
  fits <- vector("list", length(pc_test)); fail <- character(0); k <- 0L
  for (p in pc_test) {
    yv <- dat[[p]]; ok <- !is.na(yv) & !is.na(dat$age)
    ia <- ia_all[ok]; yy <- as.integer(yv[ok] == TRUE)
    n_aw <- sum(ia == 1 & yy == 1); n_awo <- sum(ia == 1 & yy == 0)
    n_rw <- sum(ia == 0 & yy == 1); n_rwo <- sum(ia == 0 & yy == 0)
    if (n_aw < TEST_FLOOR) next
    use_firth <- min(n_aw, n_awo, n_rw, n_rwo) < TEST_FLOOR
    mf <- data.frame(in_subtype = ia, y = yy, age = dat$age[ok], sex = dat$sex[ok])
    fit <- tryCatch({
      if (use_firth) logistf(in_subtype ~ y + age + sex, data = mf, control = logistf.control(maxit = 250))
      else glm(in_subtype ~ y + age + sex, data = mf, family = binomial())
    }, error = function(e) { fail <<- c(fail, sprintf("%s/%s: %s", a, p, e$message)); NULL })
    if (is.null(fit)) next
    if (use_firth) { b <- coef(fit)[["y"]]; pv <- fit$prob[["y"]]; lo <- fit$ci.lower[["y"]]; hi <- fit$ci.upper[["y"]] }
    else { sm <- summary(fit)$coefficients; b <- sm["y", "Estimate"]; se <- sm["y", "Std. Error"]
           pv <- sm["y", "Pr(>|z|)"]; lo <- b - 1.96 * se; hi <- b + 1.96 * se }
    k <- k + 1L
    fits[[k]] <- data.table(subtype = a, phecode = p, n_subtype_with = n_aw, n_subtype_without = n_awo,
      n_rest_with = n_rw, n_rest_without = n_rwo, engine = if (use_firth) "firth" else "glm",
      OR = exp(b), OR_lo = exp(lo), OR_hi = exp(hi), p = pv, well_powered = (n_aw >= WELLPOW))
  }
  res <- rbindlist(fits[seq_len(k)]); res[, q := p.adjust(p, method = "BH")]; out_list[[a]] <- res
  cat(sprintf("=== %s: tested %d | FDR-sig %d | FDR-sig & well-powered %d | failures %d\n",
              a, nrow(res), sum(res$q < 0.05), sum(res$q < 0.05 & res$well_powered), length(fail)))
  if (length(fail)) writeLines(fail, file.path(LOG_DIR, sprintf("fail_%s.txt", a)))
}
out <- rbindlist(out_list)
fwrite(out, file.path(OUT_DIR, "phecode_enrichment_g1_by_subtype.csv"))

for (a in c("A", "B", "C", "D"))
  fwrite(out[subtype == a & q < 0.05 & OR > 1][order(-well_powered, p)], file.path(OUT_DIR, sprintf("enriched_%s.csv", a)))
sig <- out[q < 0.05 & well_powered & OR > 1][order(subtype, p)]
fwrite(sig, file.path(OUT_DIR, "wellpowered_enriched_signature.csv"))
cat("[written] wellpowered_enriched_signature.csv  rows:", nrow(sig), "\n")
cat("[done]\n")
