# Each subtype vs the other three combined on the full NMR metabolomic and Olink proteomic panels (group-1 cases,
# both sexes): linear regression adjusted for age at baseline, sex and assessment center; BH FDR within subtype
# and platform.
library(data.table)
source("unsupervised_subtyping/_subtype_map.R")
setwd(Sys.getenv("PROJECT_DIR", "."))
OUT_DIR <- "output/subtyping/evidence_omics"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
MIN_N_TOTAL <- 40L   # participants with a value and an assessment center
MIN_N_IN    <- 10L   # of whom in the subtype tested

labels_of <- function(sx) {
  fl <- readRDS(sprintf("output/subtyping/clusters/cluster_assignments_g1_%s_k4.rds", sx))$final_labels
  data.table(eid = as.integer(names(fl)), sex = toupper(sx), subtype = unname(SUBTYPE_MAP[[sx]][as.character(as.integer(fl))]))
}
sub <- merge(rbind(labels_of("m"), labels_of("f")),
             as.data.table(readRDS("output/ml_ready/group1_full.rds"))[, .(eid, age_baseline, assess_centre_id)], by = "eid")

effects <- function(d, vars, platform) {
  age <- d$age_baseline; sx <- factor(d$sex); ce <- factor(d$assess_centre_id)
  rbindlist(lapply(c("A", "B", "C", "D"), function(s) {
    ia <- as.integer(d$subtype == s)
    rbindlist(lapply(vars, function(v) {
      y <- d[[v]]; ok <- is.finite(y) & !is.na(ce)
      if (sum(ok) < MIN_N_TOTAL || sum(ia[ok] == 1) < MIN_N_IN) return(NULL)
      co <- summary(lm(y[ok] ~ ia[ok] + age[ok] + sx[ok] + ce[ok]))$coefficients
      data.table(platform = platform, subtype = s, biomarker = v, beta = co[2, "Estimate"], p = co[2, "Pr(>|t|)"], n = sum(ok))
    }))
  }))
}

# biomarker groups (display order of the signature tables)
nmr_class <- function(v) {
  g  <- rep("Glycolysis, ketones & other", length(v))
  sz <- grepl("^(XXL|XL|L|M|S|XS)_", v)
  g[sz & grepl("VLDL", v)]                    <- "VLDL subclasses"
  g[sz & grepl("LDL", v) & !grepl("VLDL", v)] <- "LDL subclasses"
  g[sz & grepl("HDL", v)]                     <- "HDL subclasses"
  g[grepl("^IDL_", v)]                        <- "IDL subclasses"
  g[grepl("_size$", v) | grepl("^Apo", v)]    <- "Apolipoproteins & particle size"
  g[grepl("FA$|^MUFA|^SFA|^PUFA|^LA$|^DHA|Omega|FA_pct|by_MUFA|Omega_6_by", v)] <- "Fatty acids"
  g[v %in% c("Ala", "Gln", "Gly", "His", "Ile", "Leu", "Val", "Phe", "Tyr", "Total_BCAA")] <- "Amino acids"
  o <- g == "Glycolysis, ketones & other"
  g[o & grepl("_C$|_C_", v)] <- "Cholesterol (totals)"
  g[g == "Glycolysis, ketones & other" & grepl("_TG$|_TG_", v)] <- "Triglycerides (totals)"
  g[g == "Glycolysis, ketones & other" & grepl("_PL$|_PL_|_CE$|_FC$|_L$|_P$|_CE_|_FC_", v)] <- "Phospholipids & other lipids"
  g
}
NMR_ORDER <- c("VLDL subclasses", "IDL subclasses", "LDL subclasses", "HDL subclasses", "Apolipoproteins & particle size",
               "Cholesterol (totals)", "Triglycerides (totals)", "Phospholipids & other lipids", "Fatty acids",
               "Amino acids", "Glycolysis, ketones & other")
OLINK_ORDER <- c("Cardiometabolic", "Cardiometabolic_II", "Inflammation", "Inflammation_II",
                 "Neurology", "Neurology_II", "Oncology", "Oncology_II")
norm_key <- function(x) vapply(strsplit(tolower(x), "[_;]"), function(t) paste(sort(t), collapse = "_"), character(1))

row_order <- function(dt, group_levels) {
  dt <- copy(dt)
  dt[, strip := factor(group, levels = group_levels)]
  ord <- dt[, .(mb = mean(beta)), by = .(strip, biomarker)][order(strip, -mb)]
  dt[, biomarker := factor(biomarker, levels = rev(ord$biomarker))]
  dt[, subtype := factor(subtype, levels = c("A", "B", "C", "D"))]
  setorderv(dt, c("strip", "biomarker", "subtype"), na.last = TRUE)
  dt
}
# full results, and biomarkers reaching FDR < 0.05 in at least one subtype (with a panel/class group)
write_platform <- function(res, group_levels, stem) {
  fwrite(row_order(res, group_levels), file.path(OUT_DIR, sprintf("omics_fullpanel_%s_adjusted.csv", stem)))
  sig <- res[, any(q < 0.05), by = biomarker][V1 == TRUE, biomarker]
  fwrite(row_order(res[biomarker %in% sig & !is.na(group)], group_levels),
         file.path(OUT_DIR, sprintf("omics_signature_%s_plotdata.csv", stem)))
}

nmr <- fread("input/analysis_ready/metabolomics_group1.csv")
rN  <- effects(merge(sub, nmr, by = "eid"), setdiff(names(nmr), "eid"), "NMR")
rN[, q := p.adjust(p, "BH"), by = subtype]
rN[, group := nmr_class(biomarker)]
write_platform(rN, NMR_ORDER, "nmr")

ol <- fread("input/analysis_ready/proteomics_group1.csv")
rO <- effects(merge(sub, ol, by = "eid"), setdiff(names(ol), "eid"), "Olink")
rO[, q := p.adjust(p, "BH"), by = subtype]
pm <- fread("input/reference/olink_panel_mapping.csv")
pm[, key := norm_key(assay)]; pm <- unique(pm, by = "key")
rO[, key := norm_key(biomarker)]
rO <- merge(rO, pm[, .(key, group = panel)], by = "key", all.x = TRUE)
write_platform(rO, OLINK_ORDER, "olink")
