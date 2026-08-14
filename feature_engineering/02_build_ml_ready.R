#!/usr/bin/env Rscript
# =============================================================================
# 02_build_ml_ready.R
# =============================================================================
#
# Turns the master feature manifest into the frozen analysis-ready matrices the
# machine-learning benchmark trains on, and into the exclusion lists that tell
# each model family which features it may use.
#
# Stages (run in order; each depends on the one before):
#
#   A  extract      pull every manifest feature out of its source file, per cohort
#   B  derive       add the derived features (interval, ratio, dx counts, the
#                   composite smell-severity indicator)
#   C  split        write the four ML-ready subsets per cohort
#   D  prune        collapse protective PheCodes; freeze the matrices and the
#                   final manifest
#   E  tier map     define the six nested tiered models from the frozen manifest
#   F  missingness  per-feature missingness within each model's own context
#   G  collinearity flag |Spearman rho| > 0.85 pairs (diagnostic)
#   H  exclusions   the structural exclusion list the models read
#
# The prune runs before the tier map, so each model's feature list names
# n_protective_dx directly and nothing has to be reconciled afterwards.
#
# INPUT  (under $PROJECT_DIR/)
#   output/feature_manifest/master_feature_manifest.csv   (from 01_*)
#   input/eids/eids_group{1,2,3}.txt                      cohort membership
#   input/analysis_ready/                                 per-cohort sources
#
# OUTPUT (under $PROJECT_DIR/output/)
#   ml_ready/{group}_{subset}.rds        frozen matrices (immutable after stage D)
#   feature_manifest/master_feature_manifest_final.csv   what is in them
#   feature_manifest/tier_model_definitions_final.csv    which model uses what
#   feature_reports/                     logs, missingness, exclusion lists
#
# Those four are the whole output: every downstream script reads one of them.
#
# Convention: data.table; relative paths under $PROJECT_DIR; fail-visibly counts.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR))
  stop("PROJECT_DIR does not exist: ", PROJECT_DIR,
       "\n  set it to the project root, e.g. Sys.setenv(PROJECT_DIR = '/mnt/project')")
setwd(PROJECT_DIR)

IN_READY     <- "input/analysis_ready"
EIDS_DIR     <- "input/eids"
MANIFEST_DIR <- "output/feature_manifest"
VALUES_DIR   <- "output/feature_values"
ML_DIR       <- "output/ml_ready"
REPORT_DIR   <- "output/feature_reports"
for (d in c(VALUES_DIR, ML_DIR, REPORT_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

MANIFEST_FILE <- file.path(MANIFEST_DIR, "master_feature_manifest.csv")
FINAL_FILE    <- file.path(MANIFEST_DIR, "master_feature_manifest_final.csv")
TIER_FINAL    <- file.path(MANIFEST_DIR, "tier_model_definitions_final.csv")

# Per-cohort PheCode matrix; the derived comorbidity counts are computed from it
# directly rather than from the selected PheCode columns.
RAW_PHECODE <- file.path(IN_READY, "phecode_matrix_{group}.rds")

GROUPS  <- c("group1", "group2", "group3")
SUBSETS <- c("full", "olink", "nmr", "olink_nmr")
OUTCOME <- "taste_2w_strict"

# PheCode domain ranges (PheCode Map 1.2)
PHECODE_DOMAINS <- list(neuro = c(320, 349.99), psych = c(290, 319.99),
                        digestive = c(520, 579.99))

# Missingness classification thresholds (percent missing in the discovery cohort)
P1_THRESHOLD <- 90; P2_THRESHOLD <- 70; P3_THRESHOLD <- 50
P4_THRESHOLD <- 20; P5_THRESHOLD <- 5
REF_GROUP    <- "group1"
RHO_THRESHOLD <- 0.85

expand_group <- function(template, grp) gsub("\\{group\\}", grp, template)

manifest <- fread(MANIFEST_FILE)
cat("=== assembling ML-ready data from", nrow(manifest), "manifest features ===\n")

# =============================================================================
# Stage A — extract values
# =============================================================================
# Each manifest row names a source file and a column inside it. Sources whose
# path carries {group} are per-cohort and read once per cohort; sources without
# it are shared and read once, then filtered by cohort membership. Every source
# column is renamed to its feature_id so that identically-named columns in
# different files cannot collide. Sources are full-outer-joined on eid.

read_any <- function(path) {
  if (!file.exists(path)) stop("file not found: ", path)
  cat("      reading:", path, "")
  t0 <- Sys.time()
  dt <- if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    as.data.table(readRDS(path))
  } else if (grepl("\\.tsv$|\\.txt$", path, ignore.case = TRUE)) {
    fread(path, sep = "\t")
  } else {
    fread(path)
  }
  cat(sprintf("(%d x %d, %.1fs)\n", nrow(dt), ncol(dt),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  dt
}

harmonise_eid <- function(dt) {
  if (!"eid" %in% names(dt)) {
    for (alt in c("IID", "id", "FID")) {
      if (alt %in% names(dt)) { setnames(dt, alt, "eid"); return(dt) }
    }
    warning("no eid/IID/id/FID column; using first column '", names(dt)[1], "'")
    setnames(dt, names(dt)[1], "eid")
  }
  dt
}

cat("\n", strrep("=", 70), "\n[A] value extraction\n", strrep("=", 70), "\n", sep = "")

file_groups <- split(manifest, manifest$raw_file)
single_file_templates <- names(file_groups)[!grepl("\\{group\\}", names(file_groups))]
cat("  distinct sources:", length(file_groups),
    " | shared (read once):", length(single_file_templates), "\n")

single_file_cache <- list()
for (tmpl in single_file_templates) {
  cat(sprintf("  [CACHE] %s\n", tmpl))
  single_file_cache[[tmpl]] <- harmonise_eid(read_any(tmpl))
}

for (grp in GROUPS) {
  cat("\n  --- extracting ", grp, " ---\n", sep = "")
  eid_file <- file.path(EIDS_DIR, paste0("eids_", grp, ".txt"))
  if (!file.exists(eid_file)) stop("cohort eid list not found: ", eid_file)
  eids_keep <- fread(eid_file)[[1]]
  cat("    cohort size:", length(eids_keep), "\n")

  extracted_list <- list(); extract_log <- data.table()

  for (tmpl in names(file_groups)) {
    grp_manifest <- file_groups[[tmpl]]
    needed_cols  <- grp_manifest$raw_column
    is_template  <- grepl("\\{group\\}", tmpl)
    cat(sprintf("    [SRC] %s (%d cols)\n", substr(tmpl, 1, 60), length(needed_cols)))

    if (is_template) {
      raw <- harmonise_eid(read_any(expand_group(tmpl, grp)))
      n_before <- nrow(raw)
      raw <- raw[eid %in% eids_keep]
      cat(sprintf("      eid filter: %d -> %d\n", n_before, nrow(raw)))
    } else {
      cached <- single_file_cache[[tmpl]]
      available_in_cache <- intersect(needed_cols, names(cached))
      raw <- cached[eid %in% eids_keep, c("eid", available_in_cache), with = FALSE]
      cat(sprintf("      eid filter: %d -> %d\n", nrow(cached), nrow(raw)))
    }

    available    <- intersect(needed_cols, names(raw))
    missing_cols <- setdiff(needed_cols, names(raw))
    if (length(missing_cols) > 0)
      warning(length(missing_cols), " missing columns in ", tmpl, ": ",
              paste(head(missing_cols, 10), collapse = ", "))
    if (length(available) == 0) { warning("  no extractable columns; skipping"); next }

    out <- if (is_template) raw[, c("eid", available), with = FALSE] else raw
    setnames(out, available, grp_manifest$feature_id[match(available, grp_manifest$raw_column)])
    cat(sprintf("      extracted: %d samples x %d features\n", nrow(out), ncol(out) - 1))

    extracted_list[[tmpl]] <- out
    extract_log <- rbind(extract_log, data.table(
      source = tmpl, n_requested = length(needed_cols),
      n_extracted = length(available), n_missing = length(missing_cols),
      n_samples = nrow(out)))
    if (is_template) rm(raw)
  }

  if (length(extracted_list) == 0) stop("no data extracted for ", grp)
  cat(sprintf("    [MERGE] joining %d sources\n", length(extracted_list)))
  master_values <- Reduce(function(a, b) merge(a, b, by = "eid", all = TRUE), extracted_list)
  cat(sprintf("      merged: %d samples x %d columns (incl eid)\n",
              nrow(master_values), ncol(master_values)))

  actual_features <- ncol(master_values) - 1
  if (actual_features != nrow(manifest)) {
    cat(sprintf("      [WARN] expected %d features, got %d\n", nrow(manifest), actual_features))
    for (mf in setdiff(manifest$feature_id, setdiff(names(master_values), "eid")))
      cat(sprintf("        - %-30s (from %s)\n", mf, manifest[feature_id == mf, source_analysis]))
  } else {
    cat(sprintf("      [OK] all %d features present\n", actual_features))
  }

  cat("      missingness by layer:\n")
  for (lyr in c("L0", "L1", "L3", "L4")) {
    lyr_present <- intersect(manifest[layer == lyr, feature_id], names(master_values))
    if (length(lyr_present) > 0)
      cat(sprintf("        %s (%3d features): %.1f%% NA\n", lyr, length(lyr_present),
                  round(100 * mean(is.na(as.matrix(master_values[, ..lyr_present]))), 1)))
  }

  saveRDS(master_values, file.path(VALUES_DIR, paste0("master_values_", grp, ".rds")))
  fwrite(extract_log, file.path(REPORT_DIR, paste0("extract_log_", grp, ".csv")))
  rm(extracted_list, master_values); gc(verbose = FALSE)
}
rm(single_file_cache); gc(verbose = FALSE)

# =============================================================================
# Stage B — derived variables
# =============================================================================
#   years_baseline_to_taste  interval between baseline assessment and the taste
#                            questionnaire (age at questionnaire - age at baseline)
#   ApoB_ApoA1_ratio         atherogenic index; derived only if the ratio is not
#                            already a selected NMR feature
#   dx_count_{total,neuro,psych,digestive}
#                            comorbidity burden, summed over the whole PheCode
#                            matrix (not only the selected PheCodes)
#   n_protective_dx          count of FDR-significant protective PheCodes; a
#                            healthcare-utilisation proxy that replaces the
#                            individual protective columns at stage D
#   smell_2w_strict_OR       composite smell-severity binary, built with the same
#                            OR-logic as the taste outcome, so the smell side is
#                            represented at the same severity threshold as the
#                            outcome rather than by presence alone

cat("\n", strrep("=", 70), "\n[B] derived variables\n", strrep("=", 70), "\n", sep = "")

existing_features   <- manifest$feature_id
protective_phecodes <- manifest[source_analysis == "ExWAS_C" & direction == "-"]
risk_phecodes       <- manifest[source_analysis == "ExWAS_C" & direction == "+"]
cat("  PheCodes: total", nrow(manifest[source_analysis == "ExWAS_C"]),
    "| risk", nrow(risk_phecodes), "| protective", nrow(protective_phecodes), "\n")
protective_feature_ids <- protective_phecodes$feature_id
fwrite(protective_phecodes[, .(feature_id, phecode, or_value, pval, n_case_phecode, selection_reason)],
       file.path(REPORT_DIR, "protective_phecode_list.csv"))

need_apob_apoa1 <- !any(grepl("ApoB_by_ApoA1|ApoB_ApoA1", existing_features))
can_derive_apob_apoa1 <- need_apob_apoa1 &&
  ("ApoB" %in% existing_features) && ("ApoA1" %in% existing_features)
cat("  ApoB_ApoA1_ratio derivable:", can_derive_apob_apoa1, "\n")

TASTE_AGE_COL <- if ("age_at_taste" %in% existing_features) "age_at_taste" else
  if ("age" %in% existing_features) "age" else NA_character_
can_derive_years <- !is.na(TASTE_AGE_COL) && ("age_baseline" %in% existing_features)
cat("  years_baseline_to_taste derivable:", can_derive_years,
    "(taste-age column =", TASTE_AGE_COL, ")\n")

compute_dx_counts <- function(phecode_rds_path) {
  if (!file.exists(phecode_rds_path)) {
    cat("      [WARN] PheCode matrix not found:", phecode_rds_path, "\n"); return(NULL)
  }
  mat <- as.data.table(readRDS(phecode_rds_path))
  cat("      PheCode matrix:", nrow(mat), "x", ncol(mat), "\n")
  phe_cols    <- setdiff(names(mat), "eid")
  phe_numeric <- suppressWarnings(as.numeric(phe_cols)); names(phe_numeric) <- phe_cols
  valid_mask  <- !is.na(phe_numeric)
  phe_cols_valid    <- phe_cols[valid_mask]
  phe_numeric_valid <- phe_numeric[valid_mask]
  if (sum(!valid_mask) > 0) cat("      [WARN]", sum(!valid_mask), "non-numeric columns\n")

  phe_binary <- as.matrix(mat[, ..phe_cols_valid]) == 1
  phe_binary[is.na(phe_binary)] <- FALSE
  result <- data.table(eid = mat$eid, dx_count_total = as.integer(rowSums(phe_binary)))
  for (domain_name in names(PHECODE_DOMAINS)) {
    lo <- PHECODE_DOMAINS[[domain_name]][1]; hi <- PHECODE_DOMAINS[[domain_name]][2]
    domain_cols <- phe_cols_valid[phe_numeric_valid >= lo & phe_numeric_valid <= hi]
    col_name <- paste0("dx_count_", domain_name)
    if (length(domain_cols) == 0) {
      result[, (col_name) := 0L]
    } else {
      idx <- match(domain_cols, phe_cols_valid)
      result[, (col_name) := as.integer(rowSums(phe_binary[, idx, drop = FALSE]))]
      cat(sprintf("      %-22s %3d PheCodes  median=%d  max=%d  >0: %d\n",
                  col_name, length(domain_cols), median(result[[col_name]]),
                  max(result[[col_name]]), sum(result[[col_name]] > 0)))
    }
  }
  result
}

derived_rows <- list()
if (can_derive_years) {
  derived_rows[["years_baseline_to_taste"]] <- data.table(
    feature_id = "years_baseline_to_taste", layer = "L0",
    source_analysis = "covariate_derived", original_name = "age_at_taste - age_baseline",
    raw_file = NA_character_, raw_column = NA_character_,
    selection_reason = "derived_temporal", pval = NA_real_, direction = NA_character_,
    deployability_tier = "A", female_only = FALSE,
    description = "Years between baseline assessment and taste questionnaire")
}
if (can_derive_apob_apoa1) {
  derived_rows[["ApoB_ApoA1_ratio"]] <- data.table(
    feature_id = "ApoB_ApoA1_ratio", layer = "L4",
    source_analysis = "MWAS_derived", original_name = "ApoB / ApoA1",
    raw_file = NA_character_, raw_column = NA_character_,
    selection_reason = "derived_ratio", pval = NA_real_, direction = NA_character_,
    deployability_tier = "D", female_only = FALSE,
    description = "Atherogenic index (Walldius 2001 Lancet); conditional on not already in the NMR panel")
}
for (dd in list(
  list(id = "dx_count_total",     desc = "Total PheCode diagnoses (comorbidity burden)"),
  list(id = "dx_count_neuro",     desc = "Neurological PheCode count (320-349)"),
  list(id = "dx_count_psych",     desc = "Mental-disorder PheCode count (290-319)"),
  list(id = "dx_count_digestive", desc = "Digestive PheCode count (520-579)"))) {
  derived_rows[[dd$id]] <- data.table(
    feature_id = dd$id, layer = "L1", source_analysis = "ExWAS_C_derived",
    original_name = dd$id, raw_file = RAW_PHECODE, raw_column = NA_character_,
    selection_reason = "derived_comorbidity", pval = NA_real_, direction = "+",
    deployability_tier = "B", female_only = FALSE, description = dd$desc)
}
derived_rows[["n_protective_dx"]] <- data.table(
  feature_id = "n_protective_dx", layer = "L1", source_analysis = "ExWAS_C_derived",
  original_name = paste0("sum of ", nrow(protective_phecodes), " protective PheCodes (OR<1, FDR<0.05)"),
  raw_file = NA_character_, raw_column = NA_character_,
  selection_reason = "derived_utilization_proxy", pval = NA_real_, direction = "-",
  deployability_tier = "B", female_only = FALSE,
  description = paste0("Count of ", nrow(protective_phecodes),
                       " FDR-significant protective PheCodes; healthcare-utilisation proxy."))

SMELL_ITEMS <- c("smell_any", "smell_time", "smell_extent")
can_derive_smell_or <- all(SMELL_ITEMS %in% existing_features)
cat("  smell_2w_strict_OR derivable:", can_derive_smell_or, "\n")
if (can_derive_smell_or) {
  derived_rows[["smell_2w_strict_OR"]] <- data.table(
    feature_id = "smell_2w_strict_OR", layer = "L0",
    source_analysis = "covariate_derived",
    original_name = "smell_any==1 AND (smell_time>=1 OR smell_extent==1)",
    raw_file = NA_character_, raw_column = NA_character_,
    selection_reason = "derived_smell_severity_OR", pval = NA_real_,
    direction = NA_character_, deployability_tier = "A", female_only = FALSE,
    description = paste("Composite smell-severity binary mirroring the",
                        "taste_2w_strict OR-logic"))
}

# smell_any==1 AND (>=2 weeks OR daily-life impact); controls are smell_any==0.
# A participant reporting smell change with both sub-items missing stays NA
# rather than being counted as a mild case.
derive_smell_2w_strict_OR <- function(d) {
  sa <- d$smell_any; st <- d$smell_time; se <- d$smell_extent
  is_case <- !is.na(sa) & sa == 1L & ((!is.na(st) & st >= 1L) | (!is.na(se) & se == 1L))
  is_ctrl <- !is.na(sa) & sa == 0L
  out <- rep(NA_integer_, length(sa))
  out[is_case] <- 1L; out[is_ctrl] <- 0L
  out
}

derived_manifest <- rbindlist(derived_rows, fill = TRUE)
cat("  derived features:", nrow(derived_manifest), "\n")

dx_count_report <- list(); derive_summary <- list()
for (grp in GROUPS) {
  cat("\n  --- deriving ", grp, " ---\n", sep = "")
  rds_path <- file.path(VALUES_DIR, sprintf("master_values_%s.rds", grp))
  if (!file.exists(rds_path)) { cat("    [WARN] not found; skipping\n"); next }
  dt <- readRDS(rds_path)
  n_rows <- nrow(dt); n_cols_before <- ncol(dt)
  cat("    values:", n_rows, "x", n_cols_before, "\n")

  if (can_derive_years && TASTE_AGE_COL %in% names(dt) && "age_baseline" %in% names(dt)) {
    dt[, years_baseline_to_taste := get(TASTE_AGE_COL) - age_baseline]
    n_valid <- sum(!is.na(dt$years_baseline_to_taste))
    cat("    years_baseline_to_taste valid:", n_valid, "/", n_rows, "\n")
    if (n_valid > 0) {
      med <- median(dt$years_baseline_to_taste, na.rm = TRUE)
      if (med < 10 || med > 20) cat("      [WARN] median outside the expected 12-17 year range\n")
    }
  }

  if (can_derive_apob_apoa1 && "ApoB" %in% names(dt) && "ApoA1" %in% names(dt)) {
    dt[, ApoB_ApoA1_ratio := fifelse(!is.na(ApoA1) & ApoA1 != 0, ApoB / ApoA1, NA_real_)]
    cat("    ApoB_ApoA1_ratio valid:", sum(!is.na(dt$ApoB_ApoA1_ratio)), "/", n_rows, "\n")
  }

  cat("    dx_count variables\n")
  dx_dt <- compute_dx_counts(expand_group(RAW_PHECODE, grp))
  dx_cols <- c("dx_count_total", "dx_count_neuro", "dx_count_psych", "dx_count_digestive")
  if (!is.null(dx_dt)) {
    n_before <- nrow(dt)
    dt <- merge(dt, dx_dt, by = "eid", all.x = TRUE)
    if (nrow(dt) != n_before) cat("      [WARN] row count changed:", n_before, "->", nrow(dt), "\n")
    for (dc in dx_cols) if (any(is.na(dt[[dc]]))) dt[is.na(get(dc)), (dc) := 0L]
    for (dc in dx_cols) {
      vals <- dt[[dc]]
      cat(sprintf("      %-22s mean=%.1f  median=%d  P95=%d  max=%d  >0: %d\n",
                  dc, mean(vals), median(vals), quantile(vals, 0.95), max(vals), sum(vals > 0)))
    }
    dx_count_report[[grp]] <- data.table(
      group = grp, n_samples = nrow(dt),
      dx_total_mean = round(mean(dt$dx_count_total), 2),
      dx_total_median = median(dt$dx_count_total), dx_total_max = max(dt$dx_count_total),
      dx_neuro_gt0 = sum(dt$dx_count_neuro > 0), dx_psych_gt0 = sum(dt$dx_count_psych > 0),
      dx_digest_gt0 = sum(dt$dx_count_digestive > 0))
  } else {
    for (dc in dx_cols) dt[, (dc) := NA_integer_]
  }

  prot_cols_present <- intersect(protective_feature_ids, names(dt))
  cat("    n_protective_dx: columns present", length(prot_cols_present), "/",
      length(protective_feature_ids), "\n")
  if (length(prot_cols_present) > 0) {
    prot_binary <- as.matrix(dt[, ..prot_cols_present]) == 1
    prot_binary[is.na(prot_binary)] <- FALSE
    dt[, n_protective_dx := as.integer(rowSums(prot_binary))]
    if (OUTCOME %in% names(dt))
      cat(sprintf("      case mean=%.2f  control mean=%.2f\n",
                  mean(dt[get(OUTCOME) == 1]$n_protective_dx, na.rm = TRUE),
                  mean(dt[get(OUTCOME) == 0]$n_protective_dx, na.rm = TRUE)))
  } else {
    dt[, n_protective_dx := 0L]
  }

  if (can_derive_smell_or && all(SMELL_ITEMS %in% names(dt))) {
    dt[, smell_2w_strict_OR := derive_smell_2w_strict_OR(dt)]
    cat(sprintf("    smell_2w_strict_OR: 1=%d  0=%d  NA=%d\n",
                sum(dt$smell_2w_strict_OR == 1L, na.rm = TRUE),
                sum(dt$smell_2w_strict_OR == 0L, na.rm = TRUE),
                sum(is.na(dt$smell_2w_strict_OR))))
  }

  n_cols_after <- ncol(dt)
  saveRDS(dt, file.path(VALUES_DIR, sprintf("master_values_%s_derived.rds", grp)))
  cat("    columns:", n_cols_before, "->", n_cols_after, "\n")
  derive_summary[[grp]] <- data.table(
    group = grp, n_samples = n_rows, n_features_original = n_cols_before - 1,
    n_features_derived = n_cols_after - 1, n_new_derived = n_cols_after - n_cols_before)
  rm(dt, dx_dt); gc(verbose = FALSE)
}
print(rbindlist(derive_summary))
if (length(dx_count_report) > 0) {
  dx_summary <- rbindlist(dx_count_report); print(dx_summary)
  fwrite(dx_summary, file.path(REPORT_DIR, "dx_count_summary.csv"))
}
cat("\n  feature accounting: manifest", nrow(manifest), "+ derived", nrow(derived_manifest),
    "- pruned at stage D", nrow(protective_phecodes),
    "= final", nrow(manifest) + nrow(derived_manifest) - nrow(protective_phecodes), "\n")

# =============================================================================
# Stage C — subset splitting
# =============================================================================
# Four subsets per cohort. Omics membership is decided by a sentinel feature —
# the most significant Olink (L3) and NMR (L4) feature in the manifest: a
# participant belongs to an omics subset iff that sentinel is non-missing, which
# is equivalent to having been assayed on the platform.

cat("\n", strrep("=", 70), "\n[C] subset splitting\n", strrep("=", 70), "\n", sep = "")

manifest_all <- rbindlist(list(manifest, derived_manifest), fill = TRUE)
l3_dt <- manifest_all[layer == "L3" & !is.na(pval)][order(pval)]
l4_dt <- manifest_all[layer == "L4" & !is.na(pval)][order(pval)]
SENTINEL_OLINK <- if (nrow(l3_dt) > 0) l3_dt$feature_id[1] else NA_character_
SENTINEL_NMR   <- if (nrow(l4_dt) > 0) l4_dt$feature_id[1] else NA_character_
cat("  sentinels: Olink =", SENTINEL_OLINK, "| NMR =", SENTINEL_NMR, "\n")

split_summary <- list()
for (grp in GROUPS) {
  cat("\n  --- splitting ", grp, " ---\n", sep = "")
  derived_path <- file.path(VALUES_DIR, sprintf("master_values_%s_derived.rds", grp))
  plain_path   <- file.path(VALUES_DIR, sprintf("master_values_%s.rds", grp))
  rds_path <- if (file.exists(derived_path)) derived_path else plain_path
  if (!file.exists(rds_path)) { cat("    [WARN] not found; skipping\n"); next }

  dt <- readRDS(rds_path)
  cat("    loaded:", basename(rds_path), "-", nrow(dt), "x", ncol(dt), "\n")
  if (!is.na(SENTINEL_OLINK) && !SENTINEL_OLINK %in% names(dt)) {
    cat("    [ERROR] Olink sentinel absent; skipping cohort\n"); next
  }
  if (!is.na(SENTINEL_NMR) && !SENTINEL_NMR %in% names(dt)) {
    cat("    [ERROR] NMR sentinel absent; skipping cohort\n"); next
  }

  dt[, has_olink := FALSE][, has_nmr := FALSE]
  if (!is.na(SENTINEL_OLINK)) dt[, has_olink := !is.na(get(SENTINEL_OLINK))]
  if (!is.na(SENTINEL_NMR))   dt[, has_nmr   := !is.na(get(SENTINEL_NMR))]

  n_total <- nrow(dt); n_olink <- sum(dt$has_olink); n_nmr <- sum(dt$has_nmr)
  n_olink_nmr <- sum(dt$has_olink & dt$has_nmr)
  cat(sprintf("    subset sizes: full %s | olink %s | nmr %s | olink_nmr %s\n",
              formatC(n_total, big.mark = ","), formatC(n_olink, big.mark = ","),
              formatC(n_nmr, big.mark = ","), formatC(n_olink_nmr, big.mark = ",")))
  if (n_olink == 0) cat("    [WARN] olink subset empty\n")
  if (n_nmr == 0)   cat("    [WARN] nmr subset empty\n")

  cf <- co <- cn <- con <- NA_integer_
  if (OUTCOME %in% names(dt)) {
    cf  <- sum(dt[[OUTCOME]] == 1, na.rm = TRUE)
    co  <- sum(dt[has_olink == TRUE][[OUTCOME]] == 1, na.rm = TRUE)
    cn  <- sum(dt[has_nmr == TRUE][[OUTCOME]] == 1, na.rm = TRUE)
    con <- sum(dt[has_olink & has_nmr][[OUTCOME]] == 1, na.rm = TRUE)
    cat(sprintf("    cases: full %d | olink %d | nmr %d | olink_nmr %d\n", cf, co, cn, con))
  }

  save_cols <- setdiff(names(dt), c("has_olink", "has_nmr"))
  subsets <- list(full      = dt[, ..save_cols],
                  olink     = dt[has_olink == TRUE, ..save_cols],
                  nmr       = dt[has_nmr == TRUE, ..save_cols],
                  olink_nmr = dt[has_olink & has_nmr, ..save_cols])
  for (sub_name in names(subsets)) {
    saveRDS(subsets[[sub_name]], file.path(ML_DIR, sprintf("%s_%s.rds", grp, sub_name)))
    cat(sprintf("      -> %-10s %d x %d\n", sub_name,
                nrow(subsets[[sub_name]]), ncol(subsets[[sub_name]])))
  }

  split_summary[[grp]] <- data.table(
    group = grp, n_full = n_total, n_olink = n_olink, n_nmr = n_nmr,
    n_olink_nmr = n_olink_nmr, n_features = length(save_cols) - 1,
    case_full = cf, case_olink = co, case_nmr = cn, case_olink_nmr = con)
  rm(dt, subsets); gc(verbose = FALSE)
}
split_dt <- rbindlist(split_summary)
print(split_dt)
fwrite(split_dt, file.path(REPORT_DIR, "subset_summary.csv"))

# =============================================================================
# Stage D - prune and freeze
# =============================================================================
# The only features removed from the data are the individual protective
# PheCodes, which stage B already collapsed into n_protective_dx; keeping both
# would double-count the same healthcare-utilisation signal. Missingness and
# collinearity exclusions are NOT applied to the data — they are emitted as
# per-model-family lists at stage H, so the matrices stay one frozen artefact.

cat("\n", strrep("=", 70), "\n[D] prune and freeze\n", strrep("=", 70), "\n", sep = "")

prot_features <- fread(file.path(REPORT_DIR, "protective_phecode_list.csv"))$feature_id
cat("  protective PheCode columns to remove:", length(prot_features), "\n")
prune_log <- data.table(feature = prot_features,
                        reason = "protective_phecode_collapsed_to_n_protective_dx", auto = TRUE)

file_log <- list()
for (grp in GROUPS) for (sub in SUBSETS) {
  fpath <- file.path(ML_DIR, sprintf("%s_%s.rds", grp, sub))
  if (!file.exists(fpath)) next
  dt <- readRDS(fpath); n_cols_before <- ncol(dt)
  cols_to_keep <- setdiff(names(dt), intersect(prot_features, names(dt)))
  dt_pruned <- dt[, ..cols_to_keep]; n_cols_after <- ncol(dt_pruned)
  saveRDS(dt_pruned, fpath)
  file_log[[length(file_log) + 1]] <- data.table(
    file = basename(fpath), rows = nrow(dt_pruned),
    cols_before = n_cols_before, cols_after = n_cols_after,
    cols_removed = n_cols_before - n_cols_after)
  cat("    ", basename(fpath), ": ", n_cols_before, " -> ", n_cols_after, "\n", sep = "")
  rm(dt, dt_pruned); gc(verbose = FALSE)
}

manifest_orig <- copy(manifest)
manifest_orig[, pruned := feature_id %in% prot_features]
manifest_orig[, prune_reason := NA_character_]
manifest_orig[pruned == TRUE, prune_reason := "protective_phecode_collapsed_to_n_protective_dx"]

final_orig <- manifest_orig[pruned == FALSE]
final_orig[, c("pruned", "prune_reason") := NULL]
manifest_final <- rbindlist(list(final_orig, derived_manifest), fill = TRUE)


cat(sprintf("  original %d - pruned %d + derived %d = final %d\n",
            nrow(manifest), length(prot_features), nrow(derived_manifest), nrow(manifest_final)))
fwrite(manifest_final,     FINAL_FILE)
fwrite(prune_log,             file.path(REPORT_DIR, "pruning_log.csv"))
fwrite(rbindlist(file_log),   file.path(REPORT_DIR, "pruning_file_log.csv"))
# Export and dx upload to RAP  (ml_ready/*.rds + the final manifest are the
# frozen inputs every downstream model reads)

cat("  by layer:\n");  print(manifest_final[, .N, by = layer][order(layer)])
cat("  by deployability tier:\n")
print(manifest_final[!is.na(deployability_tier), .N, by = deployability_tier][order(deployability_tier)])

# =============================================================================
# Stage E - tiered model definitions
# =============================================================================
# Six nested models, ordered by how hard their inputs are to obtain in practice.
# Genetic PCs are excluded from every tiered model: they need population-level
# genotyping plus reference-panel projection and so are not deployable. A
# PC-inclusion ablation recipe is printed so the exclusion can be tested.

cat("\n", strrep("=", 70), "\n[E] tiered model definitions\n", strrep("=", 70), "\n", sep = "")

# The manifest is already pruned, so the tier lists name n_protective_dx
# directly and need no later reconciliation.
tier_source <- copy(manifest_final)

pc_features <- grep("^PC[0-9]+$", tier_source$feature_id, value = TRUE)
cat("  genetic PCs excluded from all tiered models:", length(pc_features), "\n")
tier_source[feature_id %in% pc_features, var_role := "excluded_non_deployable"]

outcome_features <- tier_source[grepl("outcome", var_role, ignore.case = TRUE), feature_id]
if (length(outcome_features) == 0) outcome_features <- OUTCOME
excluded_features <- c(outcome_features, pc_features)
sensory_features  <- tier_source[grepl("sensory", var_role, ignore.case = TRUE), feature_id]
l0_features       <- tier_source[layer == "L0" & !feature_id %in% excluded_features, feature_id]

get_tier <- function(tv) tier_source[deployability_tier == tv & !feature_id %in% excluded_features, feature_id]
tier_A <- get_tier("A"); tier_A_female <- get_tier("A_female")
tier_B <- get_tier("B"); tier_C <- get_tier("C"); tier_D <- get_tier("D")
d_olink <- tier_source[deployability_tier == "D" & layer == "L3" & !feature_id %in% excluded_features, feature_id]
d_nmr   <- tier_source[deployability_tier == "D" & layer == "L4" & !feature_id %in% excluded_features, feature_id]

cat("  pools: L0", length(l0_features), "| A", length(tier_A), "| A_female", length(tier_A_female),
    "| B", length(tier_B), "| C", length(tier_C), "| D", length(tier_D), "\n")

models <- list(
  list(model_id = "M1_TierA", description = "Universal-access tier (questionnaire + demographics; no genotyping)",
       sample_subset = "full", tiers_label = "A",
       features = unique(c(l0_features, tier_A, tier_A_female, sensory_features))),
  list(model_id = "M2_TierAB", description = "Clinical-records tier (adds PheCode diagnoses + comorbidity indices)",
       sample_subset = "full", tiers_label = "A+B",
       features = unique(c(l0_features, tier_A, tier_A_female, tier_B, sensory_features))),
  list(model_id = "M3_TierABC", description = "Targeted-genotype tier (adds APOE e4 dose)",
       sample_subset = "full", tiers_label = "A+B+C",
       features = unique(c(l0_features, tier_A, tier_A_female, tier_B, tier_C, sensory_features))),
  list(model_id = "M4_TierD_Olink", description = "Olink proteomic tier (research-grade)",
       sample_subset = "olink", tiers_label = "A+B+C+D(Olink)",
       features = unique(c(l0_features, tier_A, tier_A_female, tier_B, tier_C, d_olink, sensory_features))),
  list(model_id = "M5_TierD_NMR", description = "NMR metabolomic tier (research-grade, no Olink)",
       sample_subset = "nmr", tiers_label = "A+B+C+D(NMR)",
       features = unique(c(l0_features, tier_A, tier_A_female, tier_B, tier_C, d_nmr, sensory_features))),
  list(model_id = "M6_TierD_Full", description = "Full omics ceiling (Olink + NMR)",
       sample_subset = "olink_nmr", tiers_label = "A+B+C+D(Full)",
       features = unique(c(l0_features, tier_A, tier_A_female, tier_B, tier_C, tier_D, sensory_features)))
)

model_defs <- list()
for (m in models) {
  feats <- setdiff(m$features, outcome_features)
  model_defs[[m$model_id]] <- data.table(
    model_id = m$model_id, description = m$description, sample_subset = m$sample_subset,
    tiers = m$tiers_label, n_features = length(feats), feature_ids = paste(feats, collapse = ";"))
}
model_defs_dt <- rbindlist(model_defs)
for (i in seq_len(nrow(model_defs_dt)))
  cat(sprintf("    %-16s %-10s %-18s %d features\n", model_defs_dt$model_id[i],
              model_defs_dt$sample_subset[i], model_defs_dt$tiers[i], model_defs_dt$n_features[i]))
fwrite(model_defs_dt, TIER_FINAL)

cat("  incremental additions:\n")
prev_feats <- character(0)
for (m in models) {
  feats <- setdiff(m$features, outcome_features)
  if (length(prev_feats) == 0) {
    cat(sprintf("    %-16s %3d features (base)\n", m$model_id, length(feats)))
  } else {
    cat(sprintf("    %-16s %3d features (+%d new)\n", m$model_id, length(feats),
                length(setdiff(feats, prev_feats))))
  }
  prev_feats <- feats
}
cat("  PC-inclusion ablation: refit M3 with c(M3 features, ",
    paste(pc_features, collapse = ", "), ") on the same CV folds and compare AUC.\n", sep = "")

# Export and dx upload to RAP  (tier_model_definitions_final.csv is what
# the modelling step reads)
# =============================================================================
# Stage F - missingness within each model's own context
# =============================================================================
# Missingness is only interpretable relative to the subset a model trains on:
# an Olink feature is 100% missing in the full cohort but near-complete in the
# Olink subset. Both a model-aware table and a raw all-features table are kept.

cat("\n", strrep("=", 70), "\n[F] missingness\n", strrep("=", 70), "\n", sep = "")

FLAG_OK <- 5; FLAG_MODERATE <- 20; FLAG_REVIEW <- 50
flag_missingness <- function(pct)
  fifelse(pct < FLAG_OK, "OK",
          fifelse(pct < FLAG_MODERATE, "MODERATE",
                  fifelse(pct < FLAG_REVIEW, "REVIEW", "HIGH")))

raw_results <- list()
for (grp in GROUPS) for (sub in SUBSETS) {
  fpath <- file.path(ML_DIR, sprintf("%s_%s.rds", grp, sub))
  if (!file.exists(fpath)) next
  dt <- readRDS(fpath); n <- nrow(dt); feat_cols <- setdiff(names(dt), "eid")
  miss_dt <- data.table(
    feature = feat_cols, group = grp, subset = sub, n_total = n,
    n_missing = vapply(feat_cols, function(f) sum(is.na(dt[[f]])), integer(1)),
    n_valid   = vapply(feat_cols, function(f) sum(!is.na(dt[[f]])), integer(1)))
  miss_dt[, pct_missing := round(n_missing / n_total * 100, 2)]
  miss_dt[, flag := flag_missingness(pct_missing)]
  raw_results[[paste(grp, sub, sep = "_")]] <- miss_dt
  rm(dt); gc(verbose = FALSE)
}
fwrite(rbindlist(raw_results), file.path(REPORT_DIR, "missingness_raw_detail.csv"))

model_results <- list()
for (grp in GROUPS) for (mi in seq_len(nrow(model_defs_dt))) {
  m_id <- model_defs_dt$model_id[mi]; m_sub <- model_defs_dt$sample_subset[mi]
  m_feats <- strsplit(model_defs_dt$feature_ids[mi], ";")[[1]]
  fpath <- file.path(ML_DIR, sprintf("%s_%s.rds", grp, m_sub))
  if (!file.exists(fpath)) next
  dt <- readRDS(fpath); n <- nrow(dt)
  available_feats <- intersect(m_feats, names(dt))
  missing_feats   <- setdiff(m_feats, names(dt))
  if (length(missing_feats) > 0 && grp == REF_GROUP)
    cat("  [WARN]", m_id, "/", grp, ":", length(missing_feats), "features not in the data\n")
  for (feat in available_feats) {
    n_miss <- sum(is.na(dt[[feat]])); pct <- round(n_miss / n * 100, 2)
    model_results[[length(model_results) + 1]] <- data.table(
      model_id = m_id, group = grp, subset = m_sub, feature = feat,
      n_total = n, n_missing = n_miss, pct_missing = pct, flag = flag_missingness(pct))
  }
  rm(dt); gc(verbose = FALSE)
}
model_miss <- rbindlist(model_results)
fwrite(model_miss, file.path(REPORT_DIR, "missingness_by_model.csv"))

ref_summary <- model_miss[group == REF_GROUP, .(
  n_features = .N, n_OK = sum(flag == "OK"), n_MODERATE = sum(flag == "MODERATE"),
  n_REVIEW = sum(flag == "REVIEW"), n_HIGH = sum(flag == "HIGH"),
  max_miss_pct = max(pct_missing), median_miss = round(median(pct_missing), 2)
), by = .(model_id, subset)][order(model_id)]
print(ref_summary)
fwrite(ref_summary, file.path(REPORT_DIR, "missingness_by_model_summary.csv"))

# =============================================================================
# Stage G - collinearity flagging
# =============================================================================
# Pairs are flagged, never auto-deleted: tree models tolerate collinearity, so
# the decision belongs to the per-model-family exclusion lists at stage H.
# Computed on the discovery cohort's full subset; features above 50% missing are
# skipped as too sparse for a stable rank correlation.

cat("\n", strrep("=", 70), "\n[G] collinearity\n", strrep("=", 70), "\n", sep = "")

dt <- readRDS(file.path(ML_DIR, sprintf("%s_full.rds", REF_GROUP)))
cat("  ", REF_GROUP, "_full: ", nrow(dt), " x ", ncol(dt), "\n", sep = "")
feat_cols <- setdiff(names(dt), c("eid", OUTCOME))
miss_pct  <- sapply(feat_cols, function(f) sum(is.na(dt[[f]])) / nrow(dt) * 100)
sparse_feats <- names(miss_pct[miss_pct > 50])
feat_cols <- setdiff(feat_cols, sparse_feats)
cat("  features correlated:", length(feat_cols), "| skipped (>50% missing):", length(sparse_feats), "\n")

mat <- as.matrix(dt[, ..feat_cols])
rank_mat <- apply(mat, 2, function(x) {
  r <- rep(NA_real_, length(x)); valid <- !is.na(x)
  r[valid] <- rank(x[valid], ties.method = "average"); r
})
cor_mat <- cor(rank_mat, use = "pairwise.complete.obs")
rm(mat, rank_mat); gc(verbose = FALSE)

flagged_pairs <- list()
for (i in 1:(ncol(cor_mat) - 1)) for (j in (i + 1):ncol(cor_mat)) {
  rho <- cor_mat[i, j]
  if (!is.na(rho) && abs(rho) > RHO_THRESHOLD) {
    f1 <- feat_cols[i]; f2 <- feat_cols[j]
    m1 <- manifest_final[feature_id == f1]; m2 <- manifest_final[feature_id == f2]
    flagged_pairs[[length(flagged_pairs) + 1]] <- data.table(
      feature_1 = f1, feature_2 = f2, rho = round(rho, 4), abs_rho = round(abs(rho), 4),
      layer_1 = if (nrow(m1) > 0) m1$layer[1] else NA_character_,
      layer_2 = if (nrow(m2) > 0) m2$layer[1] else NA_character_,
      source_1 = if (nrow(m1) > 0) m1$source_analysis[1] else NA_character_,
      source_2 = if (nrow(m2) > 0) m2$source_analysis[1] else NA_character_,
      tier_1 = if (nrow(m1) > 0) m1$deployability_tier[1] else NA_character_,
      tier_2 = if (nrow(m2) > 0) m2$deployability_tier[1] else NA_character_,
      same_source = if (nrow(m1) > 0 && nrow(m2) > 0) m1$source_analysis[1] == m2$source_analysis[1] else NA)
  }
}
corr_path <- file.path(REPORT_DIR, "correlation_flagged_pairs.csv")
if (length(flagged_pairs) > 0) {
  pairs_dt <- rbindlist(flagged_pairs)[order(-abs_rho)]
  cat("  flagged pairs:", nrow(pairs_dt),
      "| within-source:", sum(pairs_dt$same_source == TRUE, na.rm = TRUE),
      "| cross-source:", sum(pairs_dt$same_source == FALSE, na.rm = TRUE), "\n")
  fwrite(pairs_dt, corr_path)
} else {
  cat("  no pairs above the threshold\n")
  fwrite(data.table(), corr_path)
}

max_rho <- data.table(
  feature = feat_cols,
  max_abs_rho = sapply(seq_along(feat_cols), function(i) max(abs(cor_mat[i, -i]), na.rm = TRUE)))
max_rho <- merge(max_rho[order(-max_abs_rho)],
                 manifest_final[, .(feature_id, layer, source_analysis, deployability_tier)],
                 by.x = "feature", by.y = "feature_id", all.x = TRUE)
fwrite(max_rho, file.path(REPORT_DIR, "correlation_matrix_summary.csv"))
rm(dt, cor_mat); gc(verbose = FALSE)

# =============================================================================
# Stage H - feature exclusions
# =============================================================================
# Only the gradient-boosted models are fitted, and they take missing values
# as NA and learn a default split direction, so the only exclusion is
# structural: a feature missing in almost every participant of the reference
# cohort carries no usable signal. The data itself is never modified here.

cat("\n", strrep("=", 70), "\n[H] feature exclusions\n", strrep("=", 70), "\n", sep = "")

# Structural high-missingness that is informative rather than noisy, so
# the models keep it: smell duration/impact are asked only of participants who
# reported a smell change, and the NA pattern itself carries signal.
P1_KEEP <- c("smell_time", "smell_extent")
P2_KEEP <- c("fpq_pca_pc1", "tinnitus", "age_menopause")
FORCE_REMOVE <- character(0); FORCE_KEEP <- character(0)

ref_miss <- model_miss[group == REF_GROUP,
                       .(pct_missing_ref = max(pct_missing, na.rm = TRUE),
                         n_models = uniqueN(model_id)), by = feature]
classified <- merge(
  manifest_final[, .(feature_id, layer, source_analysis, deployability_tier, female_only, var_role, pval)],
  ref_miss, by.x = "feature_id", by.y = "feature", all.x = TRUE)
setnames(classified, "feature_id", "feature")
classified[is.na(pct_missing_ref), pct_missing_ref := 0]
classified[is.na(n_models), n_models := 0L]

classified[, p_tier := fcase(
  pct_missing_ref >= P1_THRESHOLD, "P1_drop_completely",
  pct_missing_ref >= P2_THRESHOLD, "P2_review_for_drop",
  pct_missing_ref >= P3_THRESHOLD, "P3_moderate_high",
  pct_missing_ref >= P4_THRESHOLD, "P4_moderate",
  pct_missing_ref >= P5_THRESHOLD, "P5_moderate",
  default = "P6_ok")]
classified[p_tier == "P2_review_for_drop" & feature %in% P2_KEEP, p_tier := "P2_keep"]
classified[p_tier == "P1_drop_completely" & feature %in% P1_KEEP, p_tier := "P1_keep"]
classified[feature %in% FORCE_REMOVE, p_tier := "FORCE_REMOVE"]
classified[feature %in% FORCE_KEEP,   p_tier := "FORCE_KEEP"]
print(classified[, .N, by = p_tier][order(p_tier)])

manifest_remove <- classified[p_tier %in% c("P1_drop_completely", "P2_review_for_drop", "FORCE_REMOVE")][order(-pct_missing_ref)]
manifest_remove[, removal_reason := fcase(
  p_tier == "P1_drop_completely", sprintf(">=%d%% missing (conditional branch / sparsity)", P1_THRESHOLD),
  p_tier == "P2_review_for_drop", sprintf("%d-%d%% missing (skip pattern / redundant alternative)", P2_THRESHOLD, P1_THRESHOLD),
  p_tier == "FORCE_REMOVE", "manual override")]

cat("  structural exclusions:", nrow(manifest_remove), "\n")

fwrite(classified[order(-pct_missing_ref)], file.path(REPORT_DIR, "tier_classification_full.csv"))
fwrite(manifest_remove[, .(feature, pct_missing_ref, p_tier, layer, source_analysis, deployability_tier, removal_reason)],
       file.path(REPORT_DIR, "manifest_removal_candidates.csv"))
fwrite(manifest_remove[, .(feature, pct_missing_ref, p_tier, layer, source_analysis, deployability_tier)],
       file.path(REPORT_DIR, "excluded_features.csv"))
n_manifest <- nrow(manifest_final)
n_outcome  <- sum(grepl("outcome", manifest_final$var_role, ignore.case = TRUE), na.rm = TRUE)
n_pcs      <- sum(grepl("^PC[0-9]+$", manifest_final$feature_id))
n_features_used <- n_manifest - nrow(manifest_remove) - n_outcome - n_pcs
cat("  features available to the models:", n_features_used, "\n")


cat("\n=== ML-ready assembly complete ===\n")
cat("  frozen matrices : ", ML_DIR, "/{group}_{subset}.rds\n", sep = "")
cat("  final manifest  : ", FINAL_FILE, "\n", sep = "")
cat("  tier definitions: ", TIER_FINAL, "\n", sep = "")
cat("\nNext: 03_prepare_modeling_features.R\n")
