#!/usr/bin/env Rscript
# =============================================================================
# 01_select_features.R
# =============================================================================
#
# Builds the master feature manifest: the single table that declares every
# variable entering the machine-learning feature matrix, where its values come
# from, and which deployability tier it belongs to.
#
# Seven selection sources feed the manifest:
#
#   PWAS        FDR-significant Olink proteins from the primary model
#   MWAS        de-redundant Nightingale NMR cluster representatives
#   ExWAS-A     FDR-significant baseline questionnaire exposures (within domain)
#   ExWAS-B     FDR-significant follow-up questionnaire exposures (within source)
#   DWAS        FDR-significant PheCodes from the disease phenome-wide scan
#   L0          demographic / lifestyle covariates, genetic PCs, outcome, smell
#   APOE        the APOE e4 dose genetic anchor
#
# Sources are combined and de-duplicated across sources by keeping the
# minimum-P occurrence of each feature.
#
# INPUT  (under $PROJECT_DIR/input/, see README for the expected file layout)
#   assoc_results/   the primary result CSV of each association scan
#   dictionaries/    optional variable dictionary (descriptions only)
#   analysis_ready/  the per-cohort matrices the values are later read from;
#                    referenced here as path templates only, not opened
#                    (exception: the PheCode matrix, opened to read its column
#                    names so result codes can be matched to matrix columns)
#
# OUTPUT (under $PROJECT_DIR/output/)
#   feature_manifest/master_feature_manifest.csv
#
# Convention: data.table; relative paths under $PROJECT_DIR; fail-visibly counts.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR))
  stop("PROJECT_DIR does not exist: ", PROJECT_DIR,
       "\n  set it to the project root, e.g. Sys.setenv(PROJECT_DIR = '/mnt/project')")
setwd(PROJECT_DIR)

IN_RESULTS <- "input/assoc_results"
IN_DICT    <- "input/dictionaries"
IN_READY   <- "input/analysis_ready"
OUT_DIR    <- "output/feature_manifest"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_FILE <- file.path(OUT_DIR, "master_feature_manifest.csv")

# Association-scan results (one primary model per scan)
PWAS_RESULT    <- file.path(IN_RESULTS, "pwas_primary.csv")
MWAS_REPS      <- file.path(IN_RESULTS, "mwas_cluster_representatives.csv")
MWAS_ANNOT     <- file.path(IN_RESULTS, "mwas_metabolite_annotation.csv")
EXWAS_A_MAIN   <- file.path(IN_RESULTS, "exwas_baseline_primary.csv")
EXWAS_A_FEMALE <- file.path(IN_RESULTS, "exwas_baseline_female_primary.csv")
EXWAS_B_MAIN   <- file.path(IN_RESULTS, "exwas_followup_primary.csv")
DWAS_RESULT    <- file.path(IN_RESULTS, "dwas_phecode_primary.csv")
EXWAS_A_DICT   <- file.path(IN_DICT,    "exwas_baseline_dictionary.csv")

# Per-cohort value sources. {group} is expanded to group1/2/3 at extraction.
RAW_PROTEOMICS   <- file.path(IN_READY, "proteomics_{group}.csv")
RAW_METABOLOMICS <- file.path(IN_READY, "metabolomics_{group}.csv")
RAW_EXWAS_A      <- file.path(IN_READY, "exwas_baseline_{group}.csv")
RAW_EXWAS_B      <- file.path(IN_READY, "exwas_followup_{group}.csv")
RAW_PHECODE      <- file.path(IN_READY, "phecode_matrix_{group}.rds")
RAW_PHENOTYPE    <- file.path(IN_READY, "phenotype_{group}.csv")
RAW_APOE         <- file.path(IN_READY, "apoe_{group}.csv")

# The PheCode matrix of the discovery cohort supplies the real column names
# (all cohorts share the same columns).
PHECODE_REF <- gsub("\\{group\\}", "group1", RAW_PHECODE)

FDR_THRESHOLD <- 0.05

# =============================================================================
# 1. PWAS — FDR-significant proteins
# =============================================================================
# Layer L3, deployability tier D (research-grade; requires the Olink platform).

build_manifest_pwas <- function() {
  cat("\n[PWAS] ", PWAS_RESULT, "\n", sep = "")
  if (!file.exists(PWAS_RESULT)) stop("PWAS result file not found: ", PWAS_RESULT)
  pwas <- fread(PWAS_RESULT)

  required_cols <- c("protein", "or", "pval", "pval_fdr", "beta", "sig_fdr", "converged")
  missing_cols <- setdiff(required_cols, names(pwas))
  if (length(missing_cols) > 0)
    stop("PWAS results missing columns: ", paste(missing_cols, collapse = ", "),
         "\n  actual columns: ", paste(names(pwas), collapse = ", "))

  cat("  proteins tested :", nrow(pwas), "\n")
  cat("  converged       :", sum(pwas$converged, na.rm = TRUE), "\n")
  cat("  FDR-significant :", sum(pwas$sig_fdr, na.rm = TRUE), "\n")

  out <- pwas[
    !is.na(pval_fdr) & pval_fdr < FDR_THRESHOLD & converged == TRUE,
    .(
      feature_id         = protein,
      layer              = "L3",
      source_analysis    = "PWAS",
      original_name      = protein,
      raw_file           = RAW_PROTEOMICS,
      raw_column         = protein,
      selection_reason   = "FDR",
      pval               = pval,
      direction          = fifelse(or > 1, "+", "-"),
      deployability_tier = "D"
    )
  ]
  cat("  selected        :", nrow(out), "\n")
  out
}

# =============================================================================
# 2. MWAS — de-redundant metabolite cluster representatives
# =============================================================================
# Layer L4, tier D. The representative table is produced upstream by the MWAS
# de-redundancy step: within each correlation cluster the representative is the
# CE-IVD-certified member whose effect direction agrees with the cluster lead
# (minimum-P member); where no certified member agrees, the lead is kept.
# FDR is applied under the full 327-metabolite panel, not within clusters.

build_manifest_mwas <- function() {
  cat("\n[MWAS] ", MWAS_REPS, "\n", sep = "")
  if (!file.exists(MWAS_REPS)) stop("representative table not found: ", MWAS_REPS)
  m4 <- fread(MWAS_REPS)
  cat("  representatives :", nrow(m4), " | clusters:", uniqueN(m4$cluster_id), "\n")

  # Significance / FDR columns carry either name depending on how the
  # representative table was written; accept both.
  if ("sig_fdr" %in% names(m4)) {
    m4[, sig_for_ml := sig_fdr]
  } else if ("is_sig_global" %in% names(m4)) {
    m4[, sig_for_ml := is_sig_global]
  } else {
    stop("no sig_fdr or is_sig_global column in the representative table")
  }

  if ("pval_fdr_global" %in% names(m4)) {
    m4[, pval_fdr_for_ml := pval_fdr_global]
  } else if ("pval_fdr" %in% names(m4)) {
    m4[, pval_fdr_for_ml := pval_fdr]
  } else {
    warning("no FDR p-value column; manifest pval_fdr will be empty")
    m4[, pval_fdr_for_ml := NA_real_]
  }

  # Every representative in this panel is FDR-significant by construction.
  n_nonsig <- sum(!m4$sig_for_ml, na.rm = TRUE)
  if (n_nonsig > 0) warning(n_nonsig, " non-significant representatives found")

  # Since every representative is significant, the only distinction the manifest
  # needs is whether it is a CE-IVD certified biomarker; that orders the panel below.
  m4[, selection_reason := fifelse(is_ce, "CE_sig", "nonCE_sig")]

  # Pathway annotation is descriptive only; absence does not change selection.
  if (file.exists(MWAS_ANNOT)) {
    ann <- fread(MWAS_ANNOT)
    ann_keep <- unique(ann[, .(metabolite, Display_Category, Nightingale_Group)],
                       by = "metabolite")
    m4 <- merge(m4, ann_keep, by.x = "protein", by.y = "metabolite", all.x = TRUE)
  } else {
    cat("  [WARN] annotation table missing; pathway columns left empty\n")
    m4[, Display_Category := NA_character_]
    m4[, Nightingale_Group := NA_character_]
  }

  out <- m4[, .(
    feature_id         = protein,
    layer              = "L4",
    source_analysis    = "MWAS",
    original_name      = protein,
    raw_file           = RAW_METABOLOMICS,
    raw_column         = protein,
    selection_reason   = selection_reason,
    pval               = pval,
    pval_fdr           = pval_fdr_for_ml,
    direction          = fifelse(log2or > 0, "+", "-"),
    deployability_tier = "D",
    cluster_id         = cluster_id,
    is_ce              = is_ce,
    sig_fdr            = sig_for_ml,
    display_category   = Display_Category,
    nightingale_group  = Nightingale_Group,
    or                 = or,
    or_lower           = or_lower,
    or_upper           = or_upper
  )]

  dup_check <- out[, .N, by = feature_id][N > 1]
  if (nrow(dup_check) > 0) {
    warning("duplicate feature_id: ", paste(dup_check$feature_id, collapse = ", "))
    out <- unique(out, by = "feature_id")
  }

  priority_order <- c("CE_sig", "nonCE_sig")
  out[, .priority := match(selection_reason, priority_order)]
  setorder(out, .priority, pval)
  out[, .priority := NULL]

  cat("  selected        :", nrow(out),
      " | CE-IVD:", sum(out$is_ce),
      " | non-CE:", sum(!out$is_ce), "\n")
  cat("  reported panel (Methods): 57 representatives, 13 CE-IVD certified\n")
  out
}

# =============================================================================
# 3. Baseline questionnaire exposures (ExWAS A)
# =============================================================================
# Layer L1, tier A (questionnaire-only, universally available). Female-only
# reproductive variables form their own FDR family and their own tier
# (A_female); three reproductive anchors are retained regardless of
# significance so the female-specific models are not left without them.

FEMALE_ANCHORS <- c("age_menarche", "age_menopause", "menopause")

build_manifest_exwas_a <- function() {
  cat("\n[ExWAS-A] ", EXWAS_A_MAIN, "\n", sep = "")
  res_main <- fread(EXWAS_A_MAIN)
  res_main[, female_only := FALSE]
  cat("  all-sample variables :", nrow(res_main),
      " | significant:", sum(res_main$significant, na.rm = TRUE), "\n")

  if (file.exists(EXWAS_A_FEMALE)) {
    res_female <- fread(EXWAS_A_FEMALE)
    res_female[, female_only := TRUE]
    cat("  female-only variables:", nrow(res_female),
        " | significant:", sum(res_female$significant, na.rm = TRUE), "\n")
    res <- rbind(res_main, res_female, fill = TRUE)
  } else {
    warning("female-only result file not found: ", EXWAS_A_FEMALE)
    res <- res_main
  }

  if (file.exists(EXWAS_A_DICT)) {
    dict <- fread(EXWAS_A_DICT)
    if (all(c("var_name", "description") %in% names(dict)))
      res <- merge(res, dict[, .(variable = var_name, dict_description = description)],
                   by = "variable", all.x = TRUE)
  } else {
    res[, dict_description := NA_character_]
  }

  selected <- res[
    significant == TRUE & converged == TRUE,
    .(
      feature_id         = variable,
      layer              = "L1",
      source_analysis    = paste0("ExWAS_A_", fdr_group),
      original_name      = variable,
      raw_file           = RAW_EXWAS_A,
      raw_column         = variable,
      selection_reason   = "FDR_within_domain",
      pval               = p_value,
      direction          = fifelse(beta > 0, "+", "-"),
      deployability_tier = fifelse(female_only, "A_female", "A"),
      domain             = fdr_group,
      var_type           = var_type,
      female_only        = female_only
    )
  ]
  cat("  selected (FDR)       :", nrow(selected),
      " | main:", sum(selected$female_only == FALSE),
      " | female-only:", sum(selected$female_only == TRUE), "\n")

  # Force-keep the reproductive anchors that FDR did not select.
  to_force <- setdiff(FEMALE_ANCHORS, selected[female_only == TRUE, feature_id])
  if (length(to_force) > 0) {
    forced <- res[variable %in% to_force,
      .(feature_id = variable, layer = "L1", source_analysis = "ExWAS_A_Reproductive (female-only)",
        original_name = variable, raw_file = RAW_EXWAS_A, raw_column = variable,
        selection_reason = "forced_female_anchor", pval = p_value,
        direction = fifelse(beta > 0, "+", "-"), deployability_tier = "A_female",
        domain = "Reproductive (female-only)", var_type = var_type, female_only = TRUE)]
    cat("  forced female anchors:", nrow(forced), "\n")
  } else {
    forced <- selected[0]
    cat("  forced female anchors: 0 (all already FDR-selected)\n")
  }

  unique(rbind(selected, forced, fill = TRUE), by = "feature_id")
}

# =============================================================================
# 4. Follow-up questionnaire exposures (ExWAS B)
# =============================================================================
# Layer L1, tier A. Same schema as the baseline scan; FDR within source.

build_manifest_exwas_b <- function() {
  cat("\n[ExWAS-B] ", EXWAS_B_MAIN, "\n", sep = "")
  res <- fread(EXWAS_B_MAIN)
  cat("  variables   :", nrow(res),
      " | converged:", sum(res$converged, na.rm = TRUE),
      " | significant:", sum(res$significant, na.rm = TRUE), "\n")

  out <- res[
    significant == TRUE & converged == TRUE,
    .(
      feature_id         = variable,
      layer              = "L1",
      source_analysis    = paste0("ExWAS_B_", fdr_group),
      original_name      = variable,
      raw_file           = RAW_EXWAS_B,
      raw_column         = variable,
      selection_reason   = "FDR_within_domain",
      pval               = p_value,
      direction          = fifelse(beta > 0, "+", "-"),
      deployability_tier = "A",
      domain             = fdr_group,
      var_type           = var_type,
      female_only        = FALSE
    )
  ]
  cat("  selected    :", nrow(out), "\n")
  out
}

# =============================================================================
# 5. DWAS — disease phenome-wide scan (PheCodes)
# =============================================================================
# Layer L1, tier B (requires linked clinical records). Result files carry
# unpadded PheCodes ("79"); matrix columns are zero-padded ("079"), so result
# codes are matched to real column names after stripping leading zeros.

build_manifest_dwas <- function() {
  cat("\n[DWAS] ", DWAS_RESULT, "\n", sep = "")
  res <- fread(DWAS_RESULT)
  cat("  PheCodes tested :", nrow(res),
      " | FDR <", FDR_THRESHOLD, ":", sum(res$pval_fdr < FDR_THRESHOLD, na.rm = TRUE), "\n")
  res[, phecode_chr := sub("\\.0$", "", as.character(phecode))]

  if (!file.exists(PHECODE_REF)) stop("PheCode matrix not found: ", PHECODE_REF)
  mat_names <- setdiff(names(readRDS(PHECODE_REF)), "eid")
  cat("  matrix columns  :", length(mat_names), "\n")

  lookup <- data.table(raw_col = mat_names,
                       match_key = sub("^0+([0-9])", "\\1", mat_names))
  res <- merge(res, lookup, by.x = "phecode_chr", by.y = "match_key", all.x = TRUE)
  cat("  matched         :", sum(!is.na(res$raw_col)),
      " | unmatched:", sum(is.na(res$raw_col)), "\n")

  out <- res[
    !is.na(pval_fdr) & pval_fdr < FDR_THRESHOLD & converged == TRUE & !is.na(raw_col),
    .(
      feature_id         = paste0("phe", phecode_chr),
      layer              = "L1",
      source_analysis    = "ExWAS_C",
      original_name      = phecode_chr,
      raw_file           = RAW_PHECODE,
      raw_column         = raw_col,
      selection_reason   = "FDR",
      pval               = pval,
      direction          = fifelse(or > 1, "+", "-"),
      deployability_tier = "B",
      phecode            = phecode,
      or_value           = or,
      n_case_phecode     = n_case_phecode,
      female_only        = FALSE
    )
  ]

  n_fdr_no_match <- sum(res$pval_fdr < FDR_THRESHOLD & is.na(res$raw_col), na.rm = TRUE)
  if (n_fdr_no_match > 0) {
    cat("  [WARN] FDR-significant without a matrix column (excluded):", n_fdr_no_match, "\n")
    print(res[pval_fdr < FDR_THRESHOLD & is.na(raw_col), .(phecode_chr, pval, pval_fdr)])
  }

  cat("  selected        :", nrow(out),
      " | risk (OR>1):", sum(out$direction == "+"),
      " | protective (OR<1):", sum(out$direction == "-"), "\n")
  out
}

# =============================================================================
# 6. L0 — covariates, genetic PCs, outcome, smell features
# =============================================================================
# Declared rather than selected: no association result is read. Lifestyle
# variables that the association models adjust away are restored here as ML
# features. The outcome is carried in the manifest so extraction picks it up;
# it is used as y, never as a predictor.

build_manifest_L0 <- function() {
  cat("\n[L0] declared covariates / outcome / smell features\n")
  L0 <- data.table(
    feature_id = c(
      # demographics
      "age_baseline", "age", "sex", "townsend", "BMI",
      # lifestyle covariates (adjusted away in the scans; restored for ML)
      "smoking", "drink", "assess_centre_id",
      # genetic principal components
      paste0("PC", 1:10),
      # outcome
      "taste_2w_strict",
      # smell features
      "smell_any", "smell_time", "smell_extent"
    ),
    deployability_tier = c(
      rep("A", 5),   # demographics
      rep("A", 3),   # lifestyle
      rep("C", 10),  # PCs require genotyping
      "outcome",
      "A", "A", "A"  # smell
    ),
    var_role = c(
      rep("covariate_demographic", 5),
      rep("covariate_lifestyle", 3),
      rep("covariate_genetic", 10),
      "outcome_primary",
      "feature_sensory", "feature_sensory", "feature_sensory"
    )
  )
  L0[, `:=`(
    layer            = "L0",
    source_analysis  = "covariate",
    original_name    = feature_id,
    raw_file         = RAW_PHENOTYPE,
    raw_column       = feature_id,
    selection_reason = "forced",
    pval             = NA_real_,
    direction        = NA_character_,
    female_only      = FALSE
  )]
  cat("  declared        :", nrow(L0), "\n")
  L0
}

# =============================================================================
# 7. APOE — single genetic anchor
# =============================================================================
# Tier C (targeted genotyping). Effect size quoted from the reported APOE
# diplotype association: OR = 1.15 per e4 allele, P = 2.71e-08.

build_manifest_apoe <- function() {
  cat("\n[APOE] declared genetic anchor\n")
  out <- data.table(
    feature_id         = "APOE_e4_dose",
    layer              = "L1",
    source_analysis    = "GWAS",
    original_name      = "APOE_e4_dose",
    raw_file           = RAW_APOE,
    raw_column         = "APOE_e4_dose",
    selection_reason   = "GWS_anchor",
    pval               = 2.71e-08,
    direction          = "+",
    deployability_tier = "C",
    female_only        = FALSE
  )
  cat("  declared        : 1 (APOE_e4_dose)\n")
  out
}

# =============================================================================
# 8. Combine, de-duplicate across sources, write
# =============================================================================

cat("=== building the master feature manifest ===\n")

sources <- list(
  L0      = build_manifest_L0(),
  APOE    = build_manifest_apoe(),
  ExWAS_A = build_manifest_exwas_a(),
  ExWAS_B = build_manifest_exwas_b(),
  DWAS    = build_manifest_dwas(),
  MWAS    = build_manifest_mwas(),
  PWAS    = build_manifest_pwas()
)

all_manifests <- rbindlist(sources, fill = TRUE)
cat("\n[combine] rows from all sources:", nrow(all_manifests), "\n")
print(all_manifests[, .N, by = .(layer, source_analysis)][order(layer, source_analysis)])
cat("\n  by deployability tier:\n")
print(all_manifests[, .N, by = deployability_tier][order(deployability_tier)])

# ---- cross-source de-duplication -------------------------------------------
# A feature can be selected by more than one scan (e.g. a metabolite that is
# also an ExWAS exposure). Keep the occurrence with the smallest P; declared
# features (P = NA) rank last, so a tested occurrence always wins.
cat("\n[dedup] cross-source duplicates\n")
all_manifests[, dup_count := .N, by = feature_id]
duplicates <- all_manifests[dup_count > 1]
if (nrow(duplicates) > 0) {
  cat("  duplicated features:", uniqueN(duplicates$feature_id), "\n")
  for (fid in unique(duplicates$feature_id)) {
    rows <- duplicates[feature_id == fid]
    cat(sprintf("    %-30s appears %d times:\n", fid, nrow(rows)))
    for (j in seq_len(nrow(rows)))
      cat(sprintf("      [%s]  P=%.2e  reason=%s\n",
                  rows$source_analysis[j], rows$pval[j], rows$selection_reason[j]))
  }
  all_manifests <- all_manifests[order(feature_id, pval, na.last = TRUE)]
  master <- all_manifests[, .SD[1], by = feature_id]
  cat("  after de-dup:", nrow(master), "unique features\n")
} else {
  cat("  none\n")
  master <- copy(all_manifests)
}
master[, dup_count := NULL]

# ---- checks -----------------------------------------------------------------
cat("\n[check] anchor features and path form\n")
anchors <- c("PON3", "LEP", "GlycA", "APOE_e4_dose",
             "age_baseline", "sex", "BMI", "smoking", "drink",
             "taste_2w_strict", "smell_any", "phe332")
missing_anchors <- setdiff(anchors, master$feature_id)
cat("  anchors present:", length(anchors) - length(missing_anchors), "/", length(anchors), "\n")
if (length(missing_anchors) > 0)
  cat("  [WARN] missing:", paste(missing_anchors, collapse = ", "), "\n")

abs_paths <- master[grepl("^[A-Za-z]:", raw_file) | grepl("^/", raw_file)]
if (nrow(abs_paths) > 0) {
  warning("absolute paths in raw_file; all sources must be relative to PROJECT_DIR:\n",
          paste(unique(abs_paths$raw_file), collapse = "\n"))
} else {
  cat("  all raw_file paths are relative\n")
}
cat("  single-file sources:", master[!grepl("\\{group\\}", raw_file), uniqueN(raw_file)], "\n")
cat("  per-cohort templates:", master[grepl("\\{group\\}", raw_file), uniqueN(raw_file)], "\n")

# ---- PheCode zero-padding ---------------------------------------------------
# A handful of PheCodes below 100 survive the leading-zero match with the
# unpadded form; align them to the padded matrix column name.
phecode_pad <- data.table(unpadded = c("79", "79.1", "53", "78", "41"),
                          padded   = c("079", "079.1", "053", "078", "041"))
n_padded <- 0
for (i in seq_len(nrow(phecode_pad))) {
  idx <- master$raw_column == phecode_pad$unpadded[i]
  if (any(idx)) {
    master[idx, raw_column := phecode_pad$padded[i]]
    n_padded <- n_padded + sum(idx)
  }
}
cat("  PheCode raw_column zero-padded:", n_padded, "\n")

# ---- write ------------------------------------------------------------------
fwrite(master, OUT_FILE)
# Export and dx upload to RAP  (master_feature_manifest.csv drives every
# downstream extraction step)

cat("\n=== master feature manifest ===\n")
print(master[, .(
  N_features = .N,
  N_FDR      = sum(grepl("FDR", selection_reason)),
  N_forced   = sum(grepl("forced", selection_reason))
), by = .(layer, deployability_tier)][order(layer, deployability_tier)])

cat("\n  file : ", OUT_FILE, "\n", sep = "")
cat("  total:", nrow(master), "features\n")
cat("    L0 (covariate) :", sum(master$layer == "L0"), "\n")
cat("    L1 (clinical)  :", sum(master$layer == "L1"), "\n")
cat("    L3 (protein)   :", sum(master$layer == "L3"), "\n")
cat("    L4 (metabolite):", sum(master$layer == "L4"), "\n")
cat("\nNext: 02_build_ml_ready.R\n")
