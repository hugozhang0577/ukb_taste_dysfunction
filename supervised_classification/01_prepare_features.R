#!/usr/bin/env Rscript
# =============================================================================
# 01_prepare_features.R  (modelling entry point)
#
# Combine the final manifest, the tier model definitions, and the
# exclusion list into a single long-format recipe consumed by the
# training script. Derives the feature list of each tiered model;
# does NOT modify ml_ready/*.rds or any manifest.
#
# Input:  feature_manifest/master_feature_manifest_final.csv
#         feature_manifest/tier_model_definitions_final.csv
#         feature_reports/{excluded_features, manifest_removal_candidates}.csv
# Output: manifest/model_feature_recipes.csv          (long: model x feature)
#         feature_manifest/model_feature_recipes_summary.csv  (counts per model)
#         model_reports/prepare_features_log.txt
# =============================================================================

rm(list = ls()); gc()
suppressPackageStartupMessages(library(data.table))

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)

MF_DIR <- "output/feature_manifest"
P5_DIR <- "output/feature_reports"
P6_DIR <- "output/model_reports"
dir.create(P6_DIR, showWarnings = FALSE, recursive = TRUE)

while (sink.number() > 0) sink()
sink(file.path(P6_DIR, "prepare_features_log.txt"), split = TRUE)
on.exit(sink(), add = TRUE)

cat("=== prepare modeling features ===\n")

# ---- [1] load inputs --------------------------------------------------------
cat("[1] loading inputs\n")
manifest <- fread(file.path(MF_DIR, "master_feature_manifest_final.csv"))
tiers    <- fread(file.path(MF_DIR, "tier_model_definitions_final.csv"))
stopifnot("feature_id" %in% names(manifest))
stopifnot(all(c("model_id", "description", "sample_subset", "tiers", "n_features", "feature_ids") %in% names(tiers)))
cat("  final manifest:", nrow(manifest), "features | tier definitions:", nrow(tiers), "models\n")

get_fid <- function(path) {
  if (!file.exists(path)) { warning("MISSING: ", path); return(character()) }
  dt <- fread(path)
  col <- intersect(c("feature", "feature_id"), names(dt))[1]
  if (is.na(col)) stop("no feature column in ", path)
  as.character(dt[[col]])
}
excluded     <- get_fid(file.path(P5_DIR, "excluded_features.csv"))
manif_remove  <- get_fid(file.path(P5_DIR, "manifest_removal_candidates.csv"))
cat("  excluded:", length(excluded), "| manifest_removal:", length(manif_remove), "\n")

# ---- [2] always-excluded set ------------------------------------------------
outcome_feats <- manifest[var_role == "outcome_primary", feature_id]
if (length(outcome_feats) == 0) outcome_feats <- intersect("taste_2w_strict", manifest$feature_id)
pc_feats <- grep("^PC[0-9]+$", manifest$feature_id, value = TRUE)
always_excl <- unique(c(outcome_feats, pc_feats, manif_remove))
cat("[2] always-excluded:", length(always_excl), "(outcome", length(outcome_feats),
    "+ PCs", length(pc_feats), "+ manifest_removal", length(manif_remove), ")\n")

# ---- [3] parse tiers + build the per-model feature list ---------------------
cat("[3] parsing tier definitions\n")
master_fids    <- manifest$feature_id
INLINE_DERIVED <- c("years_baseline_to_taste")  # computed at train time
recipe_rows <- list(); summary_rows <- list()

for (i in seq_len(nrow(tiers))) {
  mid <- tiers$model_id[i]; desc <- tiers$description[i]
  subset <- tiers$sample_subset[i]; tlabel <- tiers$tiers[i]; n_decl <- tiers$n_features[i]
  feats <- trimws(strsplit(tiers$feature_ids[i], ";", fixed = TRUE)[[1]])
  feats <- feats[nzchar(feats)]
  if (length(feats) != n_decl) warning(sprintf("[%s] declared %d but parsed %d", mid, n_decl, length(feats)))
  missing_in_master <- setdiff(feats, c(master_fids, INLINE_DERIVED))
  if (length(missing_in_master) > 0)
    warning(sprintf("[%s] %d features not in master: %s", mid, length(missing_in_master),
                    paste(head(missing_in_master, 5), collapse = ", ")))

  base <- setdiff(feats, always_excl)
  features <- setdiff(base, setdiff(excluded, INLINE_DERIVED))

  cat(sprintf("  %-14s  subset=%-9s  base=%3d  used=%3d\n",
              mid, subset, length(base), length(features)))

  recipe_rows[[length(recipe_rows) + 1]] <- data.table(
    model_id = mid, description = desc, sample_subset = subset, tiers = tlabel,
    feature_id = features)
  summary_rows[[length(summary_rows) + 1]] <- data.table(
    model_id = mid, description = desc, sample_subset = subset, tiers = tlabel,
    n_declared = n_decl, n_base = length(base), n_features = length(features),
    n_always_excl_applied = length(intersect(feats, always_excl)))
}
recipe  <- rbindlist(recipe_rows)
summary <- rbindlist(summary_rows)

# ---- [4] write --------------------------------------------------------------
fwrite(recipe,  file.path(MF_DIR, "model_feature_recipes.csv"))
fwrite(summary, file.path(MF_DIR, "model_feature_recipes_summary.csv"))
cat("\n[4] wrote model_feature_recipes.csv (", nrow(recipe), "rows) + summary\n")
print(summary)
cat("\n=== done ===\nNext: 02_train_xgboost.R\n")
