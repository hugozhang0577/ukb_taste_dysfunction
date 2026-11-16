# Prior diagnoses not used in clustering: for each subtype vs the other three and each three-character ICD-10 code
# tested in the diagnosis-wide scan (MOFA+ input codes excluded), logistic regression of subtype membership on the
# diagnosis adjusted for age and sex (Firth when any 2 x 2 cell < 20); BH FDR within subtype.
library(data.table)
library(logistf)
source("unsupervised_subtyping/_subtype_map.R")
setwd(Sys.getenv("PROJECT_DIR", "."))
OUT_DIR <- "output/subtyping/evidence_disease"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
TEST_FLOOR <- 20L   # minimum diagnosed participants in the subtype; Firth when any cell is below it
WELLPOW    <- 50L   # well powered: >= 50 diagnosed participants in the subtype
set.seed(20260528)

labels_of <- function(sx) {
  fl <- readRDS(sprintf("output/subtyping/clusters/cluster_assignments_g1_%s_k4.rds", sx))$final_labels
  data.table(eid = as.integer(names(fl)), sex = toupper(sx), subtype = unname(SUBTYPE_MAP[[sx]][as.character(as.integer(fl))]))
}
subs <- rbind(labels_of("m"), labels_of("f"))
cov  <- as.data.table(readRDS("output/ml_ready/group1_full.rds"))[, .(eid, age)]   # age at taste assessment
dx   <- as.data.table(readRDS("input/analysis_ready/icd10_matrix_group1.rds"))

# codes tested in the diagnosis-wide scan, minus the codes that entered MOFA+
input_dx <- sub("^icd_", "", fread("output/subtyping/inputs/view_partition.csv")[view == "prior_dx", feature_id])
dx_test  <- intersect(setdiff(fread("input/assoc_results/dwas_icd10_primary.csv")$icd10, input_dx), setdiff(names(dx), "eid"))

dat <- merge(merge(subs, cov, by = "eid"), dx, by = "eid"); dat[, sex := factor(sex)]
out <- rbindlist(lapply(c("A", "B", "C", "D"), function(a) {
  ia_all <- as.integer(dat$subtype == a)
  res <- rbindlist(lapply(dx_test, function(p) {
    yv <- dat[[p]]; ok <- !is.na(yv) & !is.na(dat$age)
    ia <- ia_all[ok]; yy <- as.integer(yv[ok] == TRUE)
    n_aw <- sum(ia == 1 & yy == 1); n_awo <- sum(ia == 1 & yy == 0)
    n_rw <- sum(ia == 0 & yy == 1); n_rwo <- sum(ia == 0 & yy == 0)
    if (n_aw < TEST_FLOOR) return(NULL)
    use_firth <- min(n_aw, n_awo, n_rw, n_rwo) < TEST_FLOOR
    mf <- data.frame(in_subtype = ia, y = yy, age = dat$age[ok], sex = dat$sex[ok])
    fit <- tryCatch(if (use_firth) logistf(in_subtype ~ y + age + sex, data = mf, control = logistf.control(maxit = 250))
                    else glm(in_subtype ~ y + age + sex, data = mf, family = binomial()), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    if (use_firth) { b <- coef(fit)[["y"]]; pv <- fit$prob[["y"]]; lo <- fit$ci.lower[["y"]]; hi <- fit$ci.upper[["y"]] }
    else { sm <- summary(fit)$coefficients; b <- sm["y", "Estimate"]; se <- sm["y", "Std. Error"]
           pv <- sm["y", "Pr(>|z|)"]; lo <- b - 1.96 * se; hi <- b + 1.96 * se }
    data.table(subtype = a, icd10 = p, n_subtype_with = n_aw, n_subtype_without = n_awo,
               n_rest_with = n_rw, n_rest_without = n_rwo, engine = if (use_firth) "firth" else "glm",
               OR = exp(b), OR_lo = exp(lo), OR_hi = exp(hi), p = pv, well_powered = (n_aw >= WELLPOW))
  }))
  res[, q := p.adjust(p, method = "BH")]
}))
fwrite(out, file.path(OUT_DIR, "diagnosis_enrichment_g1_by_subtype.csv"))
fwrite(out[q < 0.05 & well_powered & OR > 1][order(subtype, p)], file.path(OUT_DIR, "wellpowered_enriched_signature.csv"))
