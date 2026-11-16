# Assemble the modeling matrices for group1/2/3 from the master feature manifest: extract
# values, derive comorbidity counts and the smell composite, split by omics availability,
# define the six nested models, and exclude features by missingness in group1.
library(data.table)
setwd(Sys.getenv("PROJECT_DIR", "."))

MANIFEST_DIR <- "output/feature_manifest"
ML_DIR       <- "output/ml_ready"
REPORT_DIR   <- "output/feature_reports"
for (d in c(ML_DIR, REPORT_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

RAW_DIAGNOSIS <- "input/analysis_ready/icd10_matrix_{group}.rds"
GROUPS  <- c("group1", "group2", "group3")
OUTCOME <- "taste_2w_strict"

# Diagnosis domains = ICD-10 chapter letters (digestive K includes oral K00-K14)
DX_DOMAINS <- list(neuro = "G", psych = "F", digestive = "K")
DX_COUNT_COLS <- c("dx_count_total", paste0("dx_count_", names(DX_DOMAINS)))

# Missingness exclusion (percent missing in group1, worst model sample)
DROP_PCT <- 90; REVIEW_PCT <- 70
DROP_KEEP   <- c("smell_time", "smell_extent")
REVIEW_KEEP <- c("fpq_pca_pc1", "tinnitus", "age_menopause")

expand_group <- function(template, grp) gsub("\\{group\\}", grp, template)
read_source  <- function(path) if (grepl("\\.rds$", path)) as.data.table(readRDS(path)) else fread(path)

manifest <- fread(file.path(MANIFEST_DIR, "master_feature_manifest.csv"))

# Inversely associated diagnoses are collapsed into one count (n_protective_dx)
protective_ids <- manifest[source_analysis == "ExWAS_C" & direction == "-", feature_id]

derived_manifest <- rbindlist(c(
  list(data.table(
    feature_id = "years_baseline_to_taste", layer = "L0", source_analysis = "covariate_derived",
    original_name = "age_at_taste - age_baseline", raw_file = NA_character_, raw_column = NA_character_,
    selection_reason = "derived_temporal", pval = NA_real_, direction = NA_character_,
    data_requirement = "A", female_only = FALSE,
    description = "Years between baseline assessment and taste questionnaire")),
  lapply(list(
    c("dx_count_total",     "Total recorded diagnoses (distinct ICD-10 3-character codes before the taste assessment; comorbidity burden)"),
    c("dx_count_neuro",     "Nervous-system diagnosis count (ICD-10 G00-G99)"),
    c("dx_count_psych",     "Mental and behavioural diagnosis count (ICD-10 F00-F99)"),
    c("dx_count_digestive", "Digestive diagnosis count (ICD-10 K00-K93, includes oral K00-K14)")),
    function(d) data.table(
      feature_id = d[1], layer = "L1", source_analysis = "ExWAS_C_derived", original_name = d[1],
      raw_file = RAW_DIAGNOSIS, raw_column = NA_character_, selection_reason = "derived_comorbidity",
      pval = NA_real_, direction = "+", data_requirement = "B", female_only = FALSE, description = d[2])),
  list(data.table(
    feature_id = "n_protective_dx", layer = "L1", source_analysis = "ExWAS_C_derived",
    original_name = paste0("sum of ", length(protective_ids),
                           " inversely associated ICD-10 diagnoses (OR<1, FDR<0.05)"),
    raw_file = NA_character_, raw_column = NA_character_,
    selection_reason = "derived_utilization_proxy", pval = NA_real_, direction = "-",
    data_requirement = "B", female_only = FALSE,
    description = paste0("Count of ", length(protective_ids),
                         " FDR-significant inversely associated ICD-10 diagnoses; healthcare-utilisation proxy.")),
  data.table(
    feature_id = "smell_2w_strict_OR", layer = "L0", source_analysis = "covariate_derived",
    original_name = "smell_any==1 AND (smell_time>=1 OR smell_extent==1)",
    raw_file = NA_character_, raw_column = NA_character_,
    selection_reason = "derived_smell_severity_OR", pval = NA_real_, direction = NA_character_,
    data_requirement = "A", female_only = FALSE,
    description = "Composite smell-severity binary mirroring the taste_2w_strict OR-logic"))),
  fill = TRUE)

# Distinct three-character ICD-10 codes recorded before the taste assessment, overall and by chapter
dx_counts <- function(path) {
  mat   <- as.data.table(readRDS(path))
  codes <- setdiff(names(mat), "eid")
  dx    <- as.matrix(mat[, ..codes]) == TRUE
  out   <- data.table(eid = mat$eid, dx_count_total = as.integer(rowSums(dx)))
  for (d in names(DX_DOMAINS))
    out[, (paste0("dx_count_", d)) :=
          as.integer(rowSums(dx[, substr(codes, 1, 1) == DX_DOMAINS[[d]], drop = FALSE]))]
  out
}

# Smell composite with the same OR logic as the taste outcome
smell_2w_strict_or <- function(d) {
  is_case <- !is.na(d$smell_any) & d$smell_any == 1L &
    ((!is.na(d$smell_time) & d$smell_time >= 1L) | (!is.na(d$smell_extent) & d$smell_extent == 1L))
  out <- rep(NA_integer_, nrow(d))
  out[is_case] <- 1L; out[!is.na(d$smell_any) & d$smell_any == 0L] <- 0L
  out
}

# Omics availability = non-missing value of the most significant protein / metabolite
SENTINEL_OLINK <- manifest[layer == "L3"][which.min(pval), feature_id]
SENTINEL_NMR   <- manifest[layer == "L4"][which.min(pval), feature_id]

sources <- split(manifest, manifest$raw_file)
for (grp in GROUPS) {
  eids_keep <- fread(file.path("input/eids", paste0("eids_", grp, ".txt")))[[1]]
  dt <- Reduce(function(a, b) merge(a, b, by = "eid", all = TRUE), lapply(names(sources), function(tmpl) {
    m   <- sources[[tmpl]]
    raw <- read_source(expand_group(tmpl, grp))[eid %in% eids_keep]
    avail <- intersect(m$raw_column, names(raw))
    out <- raw[, c("eid", avail), with = FALSE]
    setnames(out, avail, m$feature_id[match(avail, m$raw_column)])
    out
  }))

  dt[, years_baseline_to_taste := age - age_baseline]
  dt <- merge(dt, dx_counts(expand_group(RAW_DIAGNOSIS, grp)), by = "eid", all.x = TRUE)
  for (dc in DX_COUNT_COLS) dt[is.na(get(dc)), (dc) := 0L]
  prot_cols <- intersect(protective_ids, names(dt))  # matrices hold only codes seen in that group
  prot <- as.matrix(dt[, ..prot_cols]) == 1
  prot[is.na(prot)] <- FALSE
  dt[, n_protective_dx := as.integer(rowSums(prot))]
  dt[, smell_2w_strict_OR := smell_2w_strict_or(dt)]
  dt[, (prot_cols) := NULL]

  has_olink <- !is.na(dt[[SENTINEL_OLINK]]); has_nmr <- !is.na(dt[[SENTINEL_NMR]])
  saveRDS(dt,                        file.path(ML_DIR, sprintf("%s_full.rds", grp)))
  saveRDS(dt[has_olink],             file.path(ML_DIR, sprintf("%s_olink.rds", grp)))
  saveRDS(dt[has_nmr],               file.path(ML_DIR, sprintf("%s_nmr.rds", grp)))
  saveRDS(dt[has_olink & has_nmr],   file.path(ML_DIR, sprintf("%s_olink_nmr.rds", grp)))
}

manifest_final <- rbindlist(list(manifest[!feature_id %in% protective_ids], derived_manifest), fill = TRUE)
fwrite(manifest_final, file.path(MANIFEST_DIR, "master_feature_manifest_final.csv"))

# Six nested models; the outcome and genetic PCs enter none of them
excluded_ids <- c(manifest_final[grepl("outcome", var_role, ignore.case = TRUE), feature_id],
                  grep("^PC[0-9]+$", manifest_final$feature_id, value = TRUE))
use   <- manifest_final[!feature_id %in% excluded_ids]
L0    <- use[layer == "L0", feature_id]
A     <- use[data_requirement == "A", feature_id]
A_f   <- use[data_requirement == "A_female", feature_id]
B     <- use[data_requirement == "B", feature_id]
C     <- use[data_requirement == "C", feature_id]
D     <- use[data_requirement == "D", feature_id]
D_prot <- use[data_requirement == "D" & layer == "L3", feature_id]
D_met  <- use[data_requirement == "D" & layer == "L4", feature_id]

model_features <- list(
  unique(c(L0, A, A_f)),
  unique(c(L0, A, A_f, B)),
  unique(c(L0, A, A_f, B, C)),
  unique(c(L0, A, A_f, B, C, D_prot)),
  unique(c(L0, A, A_f, B, C, D_met)),
  unique(c(L0, A, A_f, B, C, D)))
models <- data.table(
  model_id = paste0("model", 1:6),
  description = c("Questionnaire and routine biochemistry",
                  "Model 1 + prior diagnoses",
                  "Model 2 + APOE genotype",
                  "Model 3 + plasma proteins",
                  "Model 3 + metabolic biomarkers",
                  "Model 3 + plasma proteins and metabolic biomarkers"),
  sample_subset = c("full", "full", "full", "olink", "nmr", "olink_nmr"),
  data_requirement = c("A", "A+B", "A+B+C", "A+B+C+D(Olink)", "A+B+C+D(NMR)", "A+B+C+D(Full)"),
  n_features = lengths(model_features),
  feature_ids = vapply(model_features, paste, character(1), collapse = ";"))
fwrite(models, file.path(MANIFEST_DIR, "model_definitions.csv"))

# Exclude features by their worst missingness across the group1 model samples
miss <- rbindlist(lapply(seq_len(nrow(models)), function(i) {
  dt <- readRDS(file.path(ML_DIR, sprintf("group1_%s.rds", models$sample_subset[i])))
  f  <- intersect(model_features[[i]], names(dt))
  data.table(feature = f,
             pct = round(vapply(f, function(x) sum(is.na(dt[[x]])), integer(1)) / nrow(dt) * 100, 2))
}))[, .(pct_missing_ref = max(pct, na.rm = TRUE)), by = feature]

excluded <- merge(miss, manifest_final[, .(feature = feature_id, layer, source_analysis, data_requirement)],
                  by = "feature")[
  (pct_missing_ref >= DROP_PCT & !feature %in% DROP_KEEP) |
  (pct_missing_ref >= REVIEW_PCT & pct_missing_ref < DROP_PCT & !feature %in% REVIEW_KEEP)]
excluded[, removal_reason := fifelse(
  pct_missing_ref >= DROP_PCT,
  sprintf(">=%d%% missing (conditional branch / sparsity)", DROP_PCT),
  sprintf("%d-%d%% missing (skip pattern / redundant alternative)", REVIEW_PCT, DROP_PCT))]
fwrite(excluded[order(-pct_missing_ref)], file.path(REPORT_DIR, "excluded_features.csv"))
