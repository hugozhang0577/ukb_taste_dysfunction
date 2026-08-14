#!/usr/bin/env Rscript
# =============================================================================
# 01_fullpanel_omics.R  (subtype molecular characterisation, full panels)
# =============================================================================
#
# One-versus-rest effect of each subtype on every measured biomarker, across the
# complete Nightingale NMR panel and the complete Olink panel. This is the
# quantity Figure 5B draws and the Results paragraph reports.
#
#   biomarker ~ is_subtype + age_baseline + sex + assessment_centre
#
# taken as the `is_subtype` coefficient, with BH-FDR within subtype and within
# platform. The two sexes are pooled and sex enters as a covariate, so a subtype's
# effect is not carried by its sex composition.
#
# Why these three covariates and no more:
#   age_baseline  age at the blood draw, not at the taste questionnaire — the
#                 assays were run on the UK Biobank baseline sample, about 15
#                 years earlier. The subtypes differ strongly in age by
#                 construction, so an unadjusted comparison would attribute much
#                 of that age difference to the subtype.
#   sex           the subtypes were fitted separately per sex; pooling without
#                 sex would reintroduce the sex composition.
#   centre        assessment centre absorbs geographic and technical batch.
# BMI, smoking and alcohol are deliberately NOT adjusted: they are subtyping
# inputs, so conditioning on them would remove part of the subtype definition.
#
# Descriptive, not confirmatory. The NMR panel contributed to the factor model
# that defined the subtypes, so the metabolome necessarily separates; the Olink
# panel did not, and its coverage is a case subset.
#
# The two platforms are NOT on a common scale: the metabolite matrix is z-scored,
# so its betas are in SD units, whereas the protein matrix is centred but not
# scaled, so its betas are in native NPX. Figure 5B gives them separate fill
# scales for this reason; do not pool the two sets of betas.
#
# CODE_DIR (env, default = current dir) must hold _subtype_map.R, captured before
# setwd(PROJECT_DIR).
#
# Input:  output/subtyping/clusters/cluster_assignments_g1_{m,f}_k4.rds
#         output/ml_ready/group1_full.rds            (age_baseline, centre)
#         input/analysis_ready/metabolomics_group1.csv
#         input/analysis_ready/proteomics_group1.csv
#         input/reference/olink_panel_mapping.csv
# Output: output/subtyping/evidence_omics/
#           omics_signature_{nmr,olink}_plotdata.csv    -> Figure 5B
#           omics_fullpanel_{nmr,olink}_adjusted.csv    -> complete matrix
# =============================================================================

CODE_DIR <- Sys.getenv("CODE_DIR", unset = getwd())
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
suppressPackageStartupMessages(library(data.table))

source(file.path(CODE_DIR, "_subtype_map.R"))   # SUBTYPE_MAP

CL      <- "output/subtyping/clusters"
OUT_DIR <- "output/subtyping/evidence_omics"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
log <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

MIN_N_TOTAL <- 40L   # rows with a finite value and a centre
MIN_N_IN    <- 10L   # rows inside the subtype being tested

# ---- [1] subtype labels, both sexes pooled ----------------------------------
labels_of <- function(sx) {
  f <- file.path(CL, sprintf("cluster_assignments_g1_%s_k4.rds", sx))
  if (!file.exists(f)) stop("cluster assignments not found: ", f)
  fl <- readRDS(f)$final_labels
  data.table(eid     = as.integer(names(fl)),
             sex     = toupper(sx),
             subtype = unname(SUBTYPE_MAP[[sx]][as.character(as.integer(fl))]))
}
sub <- rbind(labels_of("m"), labels_of("f"))
if (anyNA(sub$subtype)) stop("unmapped cluster id -- check _subtype_map.R")
log("subtype labels: ", nrow(sub), " cases")
print(sub[, .N, by = .(sex, subtype)][order(sex, subtype)])

# ---- [2] covariates ---------------------------------------------------------
cov <- as.data.table(readRDS("output/ml_ready/group1_full.rds"))[
  , .(eid, age_baseline, assess_centre_id)]
n0 <- nrow(sub)
sub <- merge(sub, cov, by = "eid")
log("after covariate join: ", n0, " -> ", nrow(sub), " (dropped ", n0 - nrow(sub), ")")
if (nrow(sub) == 0) stop("no cases survive the covariate join")

# ---- [3] the model ----------------------------------------------------------
# One lm per (subtype, biomarker); the reported effect is the is_subtype term.
effects <- function(d, vars, platform) {
  age <- d$age_baseline
  sx  <- factor(d$sex)
  ce  <- factor(d$assess_centre_id)
  rbindlist(lapply(c("A", "B", "C", "D"), function(s) {
    ia <- as.integer(d$subtype == s)
    rbindlist(lapply(vars, function(v) {
      y  <- d[[v]]
      ok <- is.finite(y) & !is.na(ce)
      if (sum(ok) < MIN_N_TOTAL || sum(ia[ok] == 1) < MIN_N_IN) return(NULL)
      fit <- lm(y[ok] ~ ia[ok] + age[ok] + sx[ok] + ce[ok])
      co  <- summary(fit)$coefficients
      data.table(platform = platform, subtype = s, biomarker = v,
                 beta = co[2, "Estimate"], p = co[2, "Pr(>|t|)"], n = sum(ok))
    }))
  }))
}

# ---- [4] biomarker grouping -------------------------------------------------
# NMR: Nightingale class. Small classes are folded so that no strip is too thin to
# label; the residual class keeps glucose, lactate, pyruvate, citrate, acetate,
# ketones, albumin, creatinine and GlycA.
nmr_class <- function(v) {
  g  <- rep("Glycolysis, ketones & other", length(v))
  sz <- grepl("^(XXL|XL|L|M|S|XS)_", v)
  g[sz & grepl("VLDL", v)]                    <- "VLDL subclasses"
  g[sz & grepl("LDL", v) & !grepl("VLDL", v)] <- "LDL subclasses"
  g[sz & grepl("HDL", v)]                     <- "HDL subclasses"
  g[grepl("^IDL_", v)]                        <- "IDL subclasses"
  g[grepl("_size$", v) | grepl("^Apo", v)]    <- "Apolipoproteins & particle size"
  g[grepl("FA$|^MUFA|^SFA|^PUFA|^LA$|^DHA|Omega|FA_pct|by_MUFA|Omega_6_by", v)] <- "Fatty acids"
  g[v %in% c("Ala","Gln","Gly","His","Ile","Leu","Val","Phe","Tyr","Total_BCAA")] <- "Amino acids"
  o <- g == "Glycolysis, ketones & other"
  g[o & grepl("_C$|_C_", v)] <- "Cholesterol (totals)"
  g[g == "Glycolysis, ketones & other" & grepl("_TG$|_TG_", v)] <- "Triglycerides (totals)"
  g[g == "Glycolysis, ketones & other" &
      grepl("_PL$|_PL_|_CE$|_FC$|_L$|_P$|_CE_|_FC_", v)] <- "Phospholipids & other lipids"
  g
}
NMR_ORDER <- c("VLDL subclasses","IDL subclasses","LDL subclasses","HDL subclasses",
               "Apolipoproteins & particle size","Cholesterol (totals)",
               "Triglycerides (totals)","Phospholipids & other lipids",
               "Fatty acids","Amino acids","Glycolysis, ketones & other")

# Olink: the eight official UK Biobank Pharma Proteomics Project panels. Composite
# multi-gene assay names are matched by normalising the gene set (lower-case, split
# on [_;], sort, rejoin) rather than by string equality.
norm_key <- function(x) vapply(strsplit(tolower(x), "[_;]"),
                               function(t) paste(sort(t), collapse = "_"), character(1))
OLINK_ORDER <- c("Cardiometabolic","Cardiometabolic_II","Inflammation","Inflammation_II",
                 "Neurology","Neurology_II","Oncology","Oncology_II")

# ---- [5] the two exports ----------------------------------------------------
# plotdata           = only the biomarkers reaching FDR in at least one subtype,
#                      i.e. exactly the rows Figure 5B draws.
# fullpanel_adjusted = the same schema for every tested biomarker, FDR-null rows
#                      included, and (Olink) assays with no official panel kept
#                      with group = NA rather than dropped.
row_order <- function(dt, group_levels) {
  dt <- copy(dt)
  dt[, strip := factor(group, levels = group_levels)]
  ord <- dt[, .(mb = mean(beta)), by = .(strip, biomarker)][order(strip, -mb)]
  dt[, biomarker := factor(biomarker, levels = rev(ord$biomarker))]
  dt[, subtype   := factor(subtype, levels = c("A","B","C","D"))]
  setorderv(dt, c("strip","biomarker","subtype"), na.last = TRUE)
  dt
}

write_platform <- function(res, group_levels, stem) {
  fwrite(row_order(res, group_levels),
         file.path(OUT_DIR, sprintf("omics_fullpanel_%s_adjusted.csv", stem)))
  sig <- res[, any(q < 0.05), by = biomarker][V1 == TRUE, biomarker]
  fwrite(row_order(res[biomarker %in% sig & !is.na(group)], group_levels),
         file.path(OUT_DIR, sprintf("omics_signature_%s_plotdata.csv", stem)))
  # Export and dx upload to RAP  (Figure 5B reads the plotdata file)
  log(stem, ": tested ", uniqueN(res$biomarker), " biomarkers; ",
      length(sig), " reach FDR in at least one subtype")
  cat("  FDR-significant per subtype:\n")
  print(res[, .(n_sig = sum(q < 0.05), n_tested = .N), by = subtype][order(subtype)])
}

# ---- [6] NMR ----------------------------------------------------------------
log("NMR: loading")
nmr  <- fread("input/analysis_ready/metabolomics_group1.csv")
mets <- setdiff(names(nmr), "eid")
dN   <- merge(sub, nmr, by = "eid"); rm(nmr); gc(verbose = FALSE)
log("NMR: ", nrow(dN), " cases x ", length(mets), " metabolites")
rN <- effects(dN, mets, "NMR")
rN[, q := p.adjust(p, "BH"), by = subtype]
rN[, group := nmr_class(biomarker)]
write_platform(rN, NMR_ORDER, "nmr")
rm(dN); gc(verbose = FALSE)

# ---- [7] Olink --------------------------------------------------------------
log("Olink: loading")
ol    <- fread("input/analysis_ready/proteomics_group1.csv")
prots <- setdiff(names(ol), "eid")
dO    <- merge(sub, ol, by = "eid"); rm(ol); gc(verbose = FALSE)
log("Olink: ", nrow(dO), " cases x ", length(prots), " assays")
rO <- effects(dO, prots, "Olink")
rO[, q := p.adjust(p, "BH"), by = subtype]

pm <- fread("input/reference/olink_panel_mapping.csv")
pm[, key := norm_key(assay)]; pm <- unique(pm, by = "key")
rO[, key := norm_key(biomarker)]
rO <- merge(rO, pm[, .(key, group = panel)], by = "key", all.x = TRUE)
log("Olink: ", rO[is.na(group), uniqueN(biomarker)],
    " assays have no official panel (kept in the full-panel export only)")
write_platform(rO, OLINK_ORDER, "olink")

log("done -> ", OUT_DIR)
