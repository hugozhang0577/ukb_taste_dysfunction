# Per-model feature lists for the six nested models: the model definitions minus the
# features excluded for missingness (the outcome and genetic PCs are already absent).
library(data.table)
setwd(Sys.getenv("PROJECT_DIR", "."))

MF_DIR <- "output/feature_manifest"
models   <- fread(file.path(MF_DIR, "model_definitions.csv"))
excluded <- fread("output/feature_reports/excluded_features.csv")$feature

recipe <- models[, .(feature_id = setdiff(strsplit(feature_ids, ";", fixed = TRUE)[[1]], excluded)),
                 by = .(model_id, description, sample_subset, data_requirement)]
fwrite(recipe, file.path(MF_DIR, "model_feature_recipes.csv"))
