# ===========================================================================
# Fig 4 — Subtype Hoadley feature x sample heatmap
#
# Cross-view overview with samples grouped by subtype, after Hoadley 2014
# (Nature 13480) Fig 1:
#   - Fig 1a: feature x sample heatmap with multi-track column annotation
#   - Fig 1c: per-cluster summary bars beneath the main heatmap
#
# Scope is the subtyping structure itself, so the column tracks are:
#        Subtype block | Depression (PHQ-9) | Sex | Age band | BMI band
# Chemosensory severity and COVID exposure are subtype characterisation rather
# than structure, and belong to Figure 5. The six row-block titles already carry
# the view-family level, so there is no separate super-group strip.
# The subtype palette and the user-facing subtype names match Figure 5.
#
# Includes ALL 5 MOFA+ input views:
#   1. Olink (proteins)         - 17 features, ~20-25% sample coverage
#   2. NMR   (metabolites)      - 57 features, ~50-55% sample coverage
#   3. Clinical continuous      - 123 features, full coverage
#   4. Clinical binary          - 84 features,  full coverage
#   5. PheCode (comorbidity)    - 31 features,  full coverage
# ===========================================================================

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ggplot2)
  library(cowplot)
  library(ggplotify)
})

ht_opt$message <- FALSE  # suppress raster info messages

# Subtype colours (A frailty / B psychosomatic / C cardiometabolic /
# D young idiopathic). Must match Fig 5 and the UMAP panels.
SUBTYPE_COL <- c(A = "#4E79A7", B = "#E15759", C = "#F28E2B", D = "#59A14F")
SEX_COL  <- c(M = "#4A6FA5", F = "#D67D8E")

# SUBTYPE_MAP is the single source of truth for the cluster-id <-> subtype-letter
# mapping. Update _subtype_map.R only when a MOFA re-fit changes cluster numbering.
source("unsupervised_subtyping/_subtype_map.R")

# Annotation track palettes (sequential or discrete)
AGE_COL   <- c("<55" = "#FFFFCC", "55-65" = "#FED976",
               "65-75" = "#FD8D3C", "75-85" = "#E31A1C", "85+"  = "#800026")
BMI_COL   <- c("<25" = "#A6DBA0", "25-30" = "#FFFFBF", ">=30" = "#D7191C")
DEP_COL   <- c("None (0-4)"       = "#FEF0D9", "Mild (5-9)"       = "#FDCC8A",
               "Moderate (10-14)" = "#FC8D59", "Severe (15+)"     = "#D7301F")
NA_COL    <- "grey90"

# --- Block schema (6 row blocks total) ---
BLOCK_ORDER <- c(
  "Olink\n(proteins)",
  "NMR\n(metabolites)",
  "Demographics &\nanthropometry",
  "Lifestyle, behaviors &\nquestionnaire",
  "Other follow-up &\nfemale-specific",
  "PheCode\ncomorbidities"
)

# Clinical block: source_analysis -> block name (3-way merge)
CLINICAL_BLOCKS <- list(
  "Demographics &\nanthropometry" = c(
    "covariate", "covariate_derived", "GWAS",
    "ExWAS_A_Demographics & SES", "ExWAS_A_Lifestyle & sleep"),
  "Lifestyle, behaviors &\nquestionnaire" = c(
    "ExWAS_A_Dietary intake & preferences", "ExWAS_A_Anthropometric & physiological", "ExWAS_A_General & mental health",
    "ExWAS_A_Sensory function & pain", "ExWAS_A_Clinical screening & treatment", "ExWAS_A_Blood biochemistry & haematology",
    "ExWAS_B_Digestive health", "ExWAS_B_Mental health", "ExWAS_B_Experience of pain"),
  "Other follow-up &\nfemale-specific" = c(
    "ExWAS_A_Oral health", "ExWAS_A_Reproductive (female-only)", "ExWAS_C_derived",
    "ExWAS_B_Food preferences", "ExWAS_B_Diet (24-hour recall)", "ExWAS_B_Cognitive function", "ExWAS_B_Work environment")
)

# PheCode ICD chapter mapping (for in-block sorting, 4 sub-categories)
phecode_to_chapter <- function(phe_ids) {
  prefix3 <- suppressWarnings(as.integer(substr(sub("^phe", "", phe_ids), 1, 3)))
  out <- ifelse(prefix3 %in% c(276, 411, 427, 433, 458), "1.Cardiometabolic",
        ifelse(prefix3 %in% c(296, 327, 332, 340, 352), "2.Mental_neuro",
        ifelse(prefix3 %in% c(471, 475, 495, 496),      "3.Respiratory",
                                                         "4.Digestive_GU")))
  out
}

OUT_DIR <- "output/figures/fig4"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


# ---------------------------------------------------------------------------
# 1. Load cluster assignments
# ---------------------------------------------------------------------------
cat("[1/7] loading cluster assignments\n")
cl_m <- readRDS("output/subtyping/clusters/cluster_assignments_g1_m_k4.rds")
cl_f <- readRDS("output/subtyping/clusters/cluster_assignments_g1_f_k4.rds")

assign_dt <- rbind(
  data.table(eid = as.integer(names(cl_m$final_labels)),
             sex = "M",
             cluster = as.integer(cl_m$final_labels),
             subtype = SUBTYPE_MAP$m[as.character(cl_m$final_labels)]),
  data.table(eid = as.integer(names(cl_f$final_labels)),
             sex = "F",
             cluster = as.integer(cl_f$final_labels),
             subtype = SUBTYPE_MAP$f[as.character(cl_f$final_labels)])
)
assign_dt[, subtype := factor(subtype, levels = c("A", "B", "C", "D"))]
assign_dt[, sex       := factor(sex,       levels = c("M", "F"))]
setorder(assign_dt, subtype, sex, eid)

cat("  total samples:", nrow(assign_dt), "\n")
print(table(assign_dt$sex, assign_dt$subtype))


# ---------------------------------------------------------------------------
# 2. Load extra covariates per sample (age, BMI, PHQ-9) for the column tracks.
#    Chemosensory and COVID tracks belong to Figure 5, not here.
# ---------------------------------------------------------------------------
cat("[2/7] loading per-sample covariates (age, BMI, PHQ-9)\n")
ml <- readRDS("output/ml_ready/group1_full.rds")
covar_cols <- intersect(c("eid", "age", "BMI", "phq9_total"), colnames(ml))
covar <- ml[, ..covar_cols]
assign_dt <- merge(assign_dt, covar, by = "eid", all.x = TRUE)
cat("  ml_ready covariate join: ", sum(!is.na(assign_dt$age)), "/", nrow(assign_dt),
    " samples with age\n", sep = "")

# Re-sort (merge re-orders); keep the canonical subtype/sex/eid order
setorder(assign_dt, subtype, sex, eid)

# Derived bands for annotation
assign_dt[, age_band := cut(age, breaks = c(-Inf, 55, 65, 75, 85, Inf),
                            labels = c("<55", "55-65", "65-75", "75-85", "85+"),
                            right = FALSE)]
assign_dt[, bmi_band := cut(BMI, breaks = c(-Inf, 25, 30, Inf),
                            labels = c("<25", "25-30", ">=30"), right = FALSE)]
# Depression (PHQ-9): standard PHQ-9 total-score severity bands
assign_dt[, dep_band := cut(phq9_total, breaks = c(-Inf, 5, 10, 15, Inf),
                            labels = c("None (0-4)", "Mild (5-9)",
                                       "Moderate (10-14)", "Severe (15+)"),
                            right = FALSE)]


# ---------------------------------------------------------------------------
# 3. Load MOFA+ inputs (sex-stratified, then combine)
# ---------------------------------------------------------------------------
cat("[3/7] loading MOFA+ inputs (sex-stratified)\n")
input_m <- readRDS("output/subtyping/inputs/clustering_input_g1_m.rds")
input_f <- readRDS("output/subtyping/inputs/clustering_input_g1_f.rds")

sample_ids_m <- names(cl_m$final_labels)
sample_ids_f <- names(cl_f$final_labels)
for (v in names(input_m)) {
  if (ncol(input_m[[v]]) != length(sample_ids_m))
    stop(sprintf("M view '%s': ncol=%d but cluster has %d samples - order mismatch!",
                 v, ncol(input_m[[v]]), length(sample_ids_m)))
  colnames(input_m[[v]]) <- sample_ids_m
}
for (v in names(input_f)) {
  if (ncol(input_f[[v]]) != length(sample_ids_f))
    stop(sprintf("F view '%s': ncol=%d but cluster has %d samples - order mismatch!",
                 v, ncol(input_f[[v]]), length(sample_ids_f)))
  colnames(input_f[[v]]) <- sample_ids_f
}

cat("  views found (M):", paste(names(input_m), collapse = ", "), "\n")


# ---------------------------------------------------------------------------
# 4. Build feature x sample matrix per view, then concatenate
# ---------------------------------------------------------------------------
cat("[4/7] building 6 row blocks\n")
all_samples_chr <- as.character(assign_dt$eid)

extract_subset <- function(mat_m, mat_f, feats, samples_chr) {
  feats <- intersect(feats, union(rownames(mat_m), rownames(mat_f)))
  out <- matrix(NA_real_, nrow = length(feats), ncol = length(samples_chr),
                dimnames = list(feats, samples_chr))
  m_cols <- intersect(colnames(mat_m), samples_chr)
  m_feats <- intersect(rownames(mat_m), feats)
  if (length(m_feats) > 0 && length(m_cols) > 0)
    out[m_feats, m_cols] <- mat_m[m_feats, m_cols, drop = FALSE]
  f_cols <- intersect(colnames(mat_f), samples_chr)
  f_feats <- intersect(rownames(mat_f), feats)
  if (length(f_feats) > 0 && length(f_cols) > 0)
    out[f_feats, f_cols] <- mat_f[f_feats, f_cols, drop = FALSE]
  out
}

manifest <- read.csv("output/feature_manifest/master_feature_manifest_final.csv",
                     encoding = "UTF-8", stringsAsFactors = FALSE)
feat_to_source <- setNames(manifest$source_analysis, manifest$feature_id)

clin_mat_m <- rbind(input_m$clinical_bin, input_m$clinical_cont)
clin_mat_f <- rbind(input_f$clinical_bin, input_f$clinical_cont)

view_mats <- list()
view_mats[["Olink\n(proteins)"]] <-
  extract_subset(input_m$olink, input_f$olink, rownames(input_m$olink), all_samples_chr)
view_mats[["NMR\n(metabolites)"]] <-
  extract_subset(input_m$nmr,   input_f$nmr,   rownames(input_m$nmr),   all_samples_chr)

for (blk_name in names(CLINICAL_BLOCKS)) {
  src_set <- CLINICAL_BLOCKS[[blk_name]]
  feats_in_blk <- names(feat_to_source)[feat_to_source %in% src_set]
  feats_in_blk <- intersect(feats_in_blk, rownames(clin_mat_m))
  if (length(feats_in_blk) == 0) next
  view_mats[[blk_name]] <- extract_subset(clin_mat_m, clin_mat_f,
                                          feats_in_blk, all_samples_chr)
}

phe_feats   <- rownames(input_m$phecode)
phe_chap    <- phecode_to_chapter(phe_feats)
phe_order   <- order(phe_chap, phe_feats)
phe_sorted  <- phe_feats[phe_order]
view_mats[["PheCode\ncomorbidities"]] <-
  extract_subset(input_m$phecode, input_f$phecode, phe_sorted, all_samples_chr)

view_mats <- view_mats[BLOCK_ORDER]
view_mats <- view_mats[!sapply(view_mats, is.null)]

for (blk in names(view_mats)) {
  m <- view_mats[[blk]]
  cat(sprintf("  [%-40s] %3d x %d  (non-NA %.1f%%)\n",
              gsub("\\n", " ", blk), nrow(m), ncol(m),
              100 * sum(!is.na(m)) / length(m)))
}

heat_mat <- do.call(rbind, view_mats)
block_sizes <- sapply(view_mats, nrow)
# Clip to +/- 1.5 so mid-strength signals (|z|~0.7) already reach ~50% saturation
heat_mat[heat_mat >  1.5] <-  1.5
heat_mat[heat_mat < -1.5] <- -1.5
cat("  combined matrix:", nrow(heat_mat), "x", ncol(heat_mat), "\n")

# Per-feature "visual saturation" = mean(|z|) over the displayed matrix.
# Used downstream as a secondary sort within the _zebra subsection so that
# rows transition smoothly from saturated to pale (avoids zebra appearance
# from random shuffling of high-|z| / low-|z| features).
feat_saturation <- data.table(
  feature_id = rownames(heat_mat),
  mean_abs_z = rowMeans(abs(heat_mat), na.rm = TRUE)
)


# ---------------------------------------------------------------------------
# 4b. Pre-compute row order per block and persist to CSV so that the G2/G3
#     supp validation heatmap re-uses the same row order.
#     - Lifestyle block (138 features, 44% of rows): SUPERVISED ordering by
#       best-driving subtype from the Wilcoxon-vs-rest discriminators. Features split into
#       4 sub-sections (A/B/C/D), each sorted by |delta| descending. This
#       block-diagonalises the Lifestyle band visually, addressing the
#       cool-tone dominance from the Ward.D2 mix. Disclosed in caption.
#     - All other row blocks: Ward.D2 hierarchical (unsupervised) as before.
# ---------------------------------------------------------------------------
cat("[4b/7] computing per-block row order (Lifestyle supervised, rest Ward.D2)\n")

LIFESTYLE_BLK <- "Lifestyle, behaviors &\nquestionnaire"
OTHER_FU_BLK  <- "Other follow-up &\nfemale-specific"
# Blocks that use SUPERVISED ordering (subtype sub-sections + _zebra tail);
# all other blocks fall back to Ward.D2 hierarchical clustering.
SUPERVISED_BLKS <- c(LIFESTYLE_BLK, OTHER_FU_BLK)

# Load the Wilcoxon-vs-rest cluster profiles (both sexes)
prof_m <- fread("output/subtyping/reports/cluster_profiles_m_k4.csv")
prof_f <- fread("output/subtyping/reports/cluster_profiles_f_k4.csv")
prof   <- rbind(prof_m[, sex := "M"], prof_f[, sex := "F"])

# NEW metric (replaces exclusivity-based approach): MOFA loading magnitude.
# Rationale (verified empirically): features that look like uniform "zebra"
# rows have IDENTICAL medians across all 4 subtypes (med_range = 0) but
# non-trivial within-subtype variance - abs_effect from the discriminators does not catch
# these. MOFA loadings DO: features with sum-of-squares loadings near zero
# contribute negligibly to any of the 14 latent factors that define the
# clustering, so by construction cannot discriminate subtypes.
# Spearman correlation between MOFA-SS and direct med_range = 0.66 (diagnostic).
mload_m <- readRDS("output/subtyping/mofa/mofa_loadings_g1_m.rds")
mload_f <- readRDS("output/subtyping/mofa/mofa_loadings_g1_f.rds")
collect_ss <- function(load_list) {
  # NOTE: do.call(c, named_list) prefixes names with list element ("view.feat")
  # which breaks downstream feature_id lookups. Use Reduce/c on unnamed list.
  result <- numeric()
  for (W in load_list) result <- c(result, rowSums(W^2))
  result
}
ss_m <- collect_ss(mload_m)
ss_f <- collect_ss(mload_f)
feat_ids_all <- union(names(ss_m), names(ss_f))
mofa_ss <- data.table(
  feature_id = feat_ids_all,
  mofa_ss_m  = ss_m[feat_ids_all],
  mofa_ss_f  = ss_f[feat_ids_all]
)
mofa_ss[, mofa_ss := rowMeans(cbind(mofa_ss_m, mofa_ss_f), na.rm = TRUE)]
mofa_ss[is.nan(mofa_ss), mofa_ss := 0]

# Which subtype each feature most-drives (from the discriminator abs_effect, summarising
# the sign of contribution rather than its magnitude).
prof[, subtype := ifelse(sex == "M",
                            SUBTYPE_MAP$m[as.character(cluster)],
                            SUBTYPE_MAP$f[as.character(cluster)])]
feat_subtype_eff <- prof[, .(eff = mean(abs_effect, na.rm = TRUE)),
                       by = .(feature_id, subtype)]
feat_subtype_eff[is.nan(eff), eff := NA_real_]
feat_drive <- dcast(feat_subtype_eff, feature_id ~ subtype,
                    value.var = "eff", fill = 0)
subtype_cols <- c("A", "B", "C", "D")
feat_mat <- as.matrix(feat_drive[, ..subtype_cols])
feat_drive[, best_subtype := subtype_cols[apply(feat_mat, 1, which.max)]]
feat_drive[, best_eff  := apply(feat_mat, 1, max, na.rm = TRUE)]
feat_drive <- merge(feat_drive, mofa_ss[, .(feature_id, mofa_ss)],
                    by = "feature_id", all.x = TRUE)
feat_drive[is.na(mofa_ss), mofa_ss := 0]

# Threshold: MOFA SS below this means the feature contributes ~nothing to any
# of the 14 latent factors that define subtypes. Truly zebra.
# Diagnostic: ~10 features SS<0.05 (vascular_stroke, milk_type_never_rarely,
# coeliac_gluten, never_eat_*, hot_drink_temp, stressor_death_close); rising
# to ~40 at SS<0.20. Default 0.10 captures the truly broken ones.
MOFA_SS_THRESHOLD <- 0.10

supervised_order <- function(blk_feats) {
  dt <- data.table(feature_id = blk_feats)
  dt <- merge(dt, feat_drive[, .(feature_id, best_subtype, best_eff, mofa_ss)],
              by = "feature_id", all.x = TRUE)
  dt <- merge(dt, feat_saturation, by = "feature_id", all.x = TRUE)
  dt[is.na(mofa_ss),     mofa_ss     := 0]
  dt[is.na(best_eff),    best_eff    := 0]
  dt[is.na(mean_abs_z),  mean_abs_z  := 0]
  dt[is.na(best_subtype),                best_subtype := "_unranked"]
  dt[mofa_ss < MOFA_SS_THRESHOLD,     best_subtype := "_zebra"]
  dt[, best_subtype := factor(best_subtype,
                            levels = c("A", "B", "C", "D",
                                       "_zebra", "_unranked"))]
  # Within-section sort:
  # - pure (A/B/C/D): by MOFA SS desc (strongest discriminator at top)
  # - _zebra: by mean|z| desc (most-saturated rows at top transitioning from
  #     pure D, palest rows at very bottom -> smooth visual gradient,
  #     eliminates intra-section zebra from random saturation shuffling)
  dt[, sort_key := fcase(
        best_subtype %in% c("A", "B", "C", "D"), -mofa_ss,
        best_subtype == "_zebra",                -mean_abs_z,
        default = 0)]
  setorder(dt, best_subtype, sort_key)
  dt$feature_id
}

compute_block_row_order <- function(mat, blocks) {
  out_order <- integer(0)
  out_meta  <- list()
  start <- 1
  for (blk in names(blocks)) {
    n   <- blocks[blk]
    end <- start + n - 1
    sub <- mat[start:end, , drop = FALSE]
    blk_feats <- rownames(sub)
    if (blk %in% SUPERVISED_BLKS) {
      ordered_feats <- supervised_order(blk_feats)
      ordered_idx   <- match(ordered_feats, blk_feats)
      out_order     <- c(out_order, (start:end)[ordered_idx])
      out_meta[[blk]] <- "supervised_by_subtype"
    } else {
      sub_imp <- sub; sub_imp[is.na(sub_imp)] <- 0
      hc <- hclust(dist(sub_imp), method = "ward.D2")
      out_order <- c(out_order, (start:end)[hc$order])
      out_meta[[blk]] <- "ward.D2"
    }
    start <- end + 1
  }
  attr(out_order, "method_per_block") <- out_meta
  out_order
}
row_order_vec <- compute_block_row_order(heat_mat, block_sizes)

# Persist with subtype-section column (NA for non-Lifestyle blocks).
# Apply the same exclusivity-threshold reclassification used in supervised_order.
ordered_feat <- rownames(heat_mat)[row_order_vec]
ordered_blk  <- rep(names(block_sizes), block_sizes)[row_order_vec]
mi <- match(ordered_feat, feat_drive$feature_id)
subtype_section    <- feat_drive$best_subtype[mi]
mofa_ss_lookup  <- feat_drive$mofa_ss[mi]
subtype_section[is.na(subtype_section)] <- "_unranked"
subtype_section[!is.na(mofa_ss_lookup) & mofa_ss_lookup < MOFA_SS_THRESHOLD] <- "_zebra"
subtype_section[!ordered_blk %in% SUPERVISED_BLKS] <- NA  # supervised blocks only
row_order_dt <- data.table(
  feature_id        = ordered_feat,
  block             = ordered_blk,
  position_in_block = unlist(lapply(block_sizes, seq_len), use.names = FALSE),
  subtype_section = subtype_section,
  ordering_method   = vapply(ordered_blk,
                              function(b) attr(row_order_vec, "method_per_block")[[b]],
                              character(1))
)
ROW_ORDER_CSV <- "output/subtyping/reports/fig4_row_order.csv"
fwrite(row_order_dt, ROW_ORDER_CSV)
cat("  saved row order:", ROW_ORDER_CSV, "(", nrow(row_order_dt), "rows)\n")
cat("  Lifestyle sub-section counts:\n")
print(row_order_dt[block %in% SUPERVISED_BLKS, .N, by = .(block, subtype_section)])


# ---------------------------------------------------------------------------
# 5. Splits + annotations (TOP = 6 tracks, LEFT = supergroup)
# ---------------------------------------------------------------------------
cat("[5/7] splits + 6-track top annotation\n")

row_split <- factor(rep(names(block_sizes), block_sizes),
                    levels = names(block_sizes))

subtype_full <- c(A = "A - Aging frailty",
               B = "B - Psychosomatic",
               C = "C - Cardiometabolic",
               D = "D - Young idiopathic")
col_split <- factor(subtype_full[as.character(assign_dt$subtype)],
                    levels = subtype_full[c("A","B","C","D")])

# Subtype track = anno_block: a coloured block per subtype slice with the subtype
# name printed inside it, so the colour key and the slice label are one element.
top_anno <- HeatmapAnnotation(
  Subtype = anno_block(
    gp        = gpar(fill = SUBTYPE_COL[c("A","B","C","D")], col = "white", lwd = 2),
    # 2-line labels so the narrow B/C blocks don't clip the descriptive name
    labels    = c("A\nAging\nfrailty", "B\nPsycho-\nsomatic",
                  "C\nCardio-\nmetabolic", "D\nYoung\nidiopathic"),
    labels_gp = gpar(col = "white", fontface = "bold", fontsize = 8.5)
  ),
  "Depression (PHQ-9)"  = assign_dt$dep_band,
  Sex                   = assign_dt$sex,
  "Age band"            = assign_dt$age_band,
  "BMI band"            = assign_dt$bmi_band,
  col = list(
    "Depression (PHQ-9)"  = DEP_COL,
    Sex                   = SEX_COL,
    "Age band"            = AGE_COL,
    "BMI band"            = BMI_COL
  ),
  na_col = NA_COL,
  annotation_height = unit.c(
    unit(15,  "mm"),  # Subtype block (prominent, holds 3-line name text)
    unit(3.2, "mm"),  # Depression (PHQ-9)
    unit(3.2, "mm"),  # Sex
    unit(3.2, "mm"),  # Age band
    unit(3.2, "mm")   # BMI band
  ),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 11, fontface = "bold"),
  show_legend = TRUE,
  gap = unit(0.8, "mm"),
  border = TRUE,
  annotation_legend_param = list(
    "Depression (PHQ-9)" = list(title_gp = gpar(fontsize = 11, fontface = "bold"),
                                labels_gp = gpar(fontsize = 11)),
    Sex       = list(title_gp = gpar(fontsize = 11, fontface = "bold"),
                     labels_gp = gpar(fontsize = 11)),
    "Age band" = list(title_gp = gpar(fontsize = 11, fontface = "bold"),
                      labels_gp = gpar(fontsize = 11), ncol = 1),
    "BMI band" = list(title_gp = gpar(fontsize = 11, fontface = "bold"),
                      labels_gp = gpar(fontsize = 11), ncol = 1)
  )
)

# Row blocks are split by BLOCK_ORDER below; the six block titles carry the
# view-family level, so no separate left-hand strip is drawn.

# 5-stop diverging ramp with vivid endpoints. Intermediate stops are 50% tints
# of the endpoints toward white.
# - Blue-purple: #3C54A5 (vivid) -> #9EAAD2 (50% tint) -> #FFFFFF
# - Red-orange : #EC2027 (vivid) -> #F69093 (50% tint) -> #FFFFFF
col_fun <- colorRamp2(c(-1.5, -0.6, 0, 0.6, 1.5),
                      c("#3C54A5", "#9EAAD2", "#FFFFFF", "#F69093", "#EC2027"))


# ---------------------------------------------------------------------------
# 6. Build heatmap
# ---------------------------------------------------------------------------
cat("[6/7] drawing heatmap\n")

ht <- Heatmap(
  heat_mat,
  name   = "z-score",
  col    = col_fun,
  na_col = "#E8E0D0",

  row_split    = row_split,
  column_split = col_split,
  row_gap    = unit(1.5, "mm"),   # uniform gap between all 6 row blocks
  column_gap = unit(2.5, "mm"),

  row_order              = row_order_vec,   # pre-computed (block-wise)
  cluster_rows           = FALSE,
  cluster_columns        = FALSE,
  cluster_row_slices     = FALSE,
  show_row_dend          = FALSE,
  use_raster             = TRUE,
  raster_quality         = 4,

  show_row_names      = FALSE,
  show_column_names   = FALSE,
  row_title_gp        = gpar(fontsize = 11, fontface = "bold", col = "#2C3E50"),
  row_title_rot       = 0,
  row_title_side      = "left",
  column_title        = NULL,   # subtype names live inside the anno_block now

  top_annotation  = top_anno,

  width  = unit(17, "cm"),
  height = unit(16, "cm"),

  heatmap_legend_param = list(
    title       = "z-score\n(saturates |z|>=1.5)",
    title_gp    = gpar(fontsize = 12, fontface = "bold"),
    labels_gp   = gpar(fontsize = 11),
    at          = c(-1.5, -0.75, 0, 0.75, 1.5),
    labels      = c("<=-1.5", "-0.75", "0", "0.75", ">=1.5"),
    legend_height = unit(3.5, "cm")
  ),
  border = TRUE
)


# ---------------------------------------------------------------------------
# 7. Save outputs — Fig 4 Panel a only (the feature x sample heatmap).
#    Panels b (MOFA UMAP) and c (cross-sex cosine) are produced by the
#    standalone script fig4_panels_bc.R as separate square figures;
#    final multi-panel assembly is done by the user in external software.
# ---------------------------------------------------------------------------
cat("[7/7] saving outputs (Fig 4 Panel a)\n")

save_heatmap <- function(open_fn) {
  open_fn()
  draw(ht, merge_legend = TRUE,
       heatmap_legend_side    = "right",
       annotation_legend_side = "right",
       padding = unit(c(4, 4, 4, 4), "mm"))
  dev.off()
}

save_heatmap(function() pdf(file.path(OUT_DIR, "fig4_feature_sample_heatmap.pdf"),
                            width = 12, height = 12))
save_heatmap(function() png(file.path(OUT_DIR, "fig4_feature_sample_heatmap.png"),
                            width = 12, height = 12, units = "in", res = 220))

cat("\n=== Done ===\n")
cat("Saved to:", OUT_DIR, "\n")
cat("  fig4_feature_sample_heatmap.{pdf,png}  (Fig 4 Panel a)\n")
cat("Matrix: ", nrow(heat_mat), " features x ", ncol(heat_mat), " samples\n", sep="")
cat("Top column tracks: Subtype (block+text) | Depression (PHQ-9) | ",
    "Sex | Age band | BMI band\n", sep = "")
cat("Subtype sizes (should be A1837 / B754 / C1625 / D1556):\n")
print(table(assign_dt$subtype))
