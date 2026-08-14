# Unsupervised multi-omic subtyping

Code for the Methods subsection *Unsupervised multi-omic subtyping*: a case-only
multi-view factor model (MOFA+) followed by k-means clustering into four
subtypes, fitted separately in men and women. Same conventions as the other
`code_script/` modules: no absolute or user-specific paths (everything anchored
on `$PROJECT_DIR`), `data.table` for R I/O, tabular and figure outputs only.

Whether the same subtypes come back is `../subtype_reproducibility/`; what they
contain is `../subtype_characterisation/`. The matrices come from
`../feature_engineering/` and are read-only here.

## Files

The file numbering is the run order; `_subtype_map.R` is sourced, not run.

| File | Role |
|---|---|
| `01_prepare_clustering_features.R` | per-sex MOFA input: five views (Olink, NMR, continuous clinical, binary clinical, PheCode), z-scored Gaussian views, missing values preserved. Takes an `M`/`F` argument |
| `02_mofa_fit.R` | MOFA+ fit (15 factors, slow convergence, fixed seed); writes factors, loadings and the variance decomposition |
| `03_cluster_on_factors.R` | k-means / PAM / HDBSCAN sweep on the active factors; k chosen by silhouette, bootstrap Jaccard and a size floor; `<sex> <k>` forces k |
| `04_cluster_discriminators.R` | cluster-versus-rest discriminators across the manifest (Wilcoxon / Fisher, BH-FDR) — the per-subtype signature the reporting tables and Figure 4 use |
| `05_reporting_tables.R` | per-subtype signature tables and the subtype baseline table |
| `06_umap_embedding.R` | two-dimensional UMAP view of the active factor space, one per sex, for Figure 4B/4C. Display only: the subtypes were found in the full factor space, not in these two dimensions, so separation in the plot is neither necessary nor sufficient evidence of a cluster |
| `_subtype_map.R` | sourced by every subtype-aware script in all three subtype modules: the single source of truth for the cluster-id to subtype-letter mapping, per sex, with a self-consistency assertion |

## Run

```bash
export PROJECT_DIR=/path/to/project
export CODE_DIR=$PWD            # so the scripts find _subtype_map.R

for SEX in M F; do
  for s in 01 02 03 04; do Rscript ${s}_*.R $SEX; done
done
Rscript 05_reporting_tables.R
Rscript 06_umap_embedding.R
```

A MOFA fit takes roughly 30–90 minutes, so `02` is worth running detached.

## Specification (matches Methods)

- **Input**: cases only, split by sex, from `output/ml_ready/{group}_full.rds`.
  The clustering input excludes the outcome, the genetic principal components,
  and `smell_time` / `smell_extent` — the last two are represented instead by
  `smell_2w_strict_OR`, the composite derived in
  `../feature_engineering/02_build_ml_ready.R` that mirrors the primary-outcome
  OR-logic. They remain available for characterisation.
- **Why per sex**: taste dysfunction is itself sex-differentiated, so the two
  pipelines are independent rather than a pooled fit with a sex covariate.
- **Model**: MOFA+ across five views; k-means on the active factors; k = 4.
- **Subtypes**: **A Aging frailty / B Psychosomatic / C Cardiometabolic /
  D Young idiopathic**. Cluster numbering from the clustering step is arbitrary
  and differs between the sexes, which is why `_subtype_map.R` exists and why no
  script hard-codes a letter.
