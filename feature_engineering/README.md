# Feature engineering

The analysis-ready feature set: which variables enter the modelling matrices,
where their values come from, and which of the six tiered models may use each of
them. Same conventions as the other `code_script/` modules: no absolute or
user-specific paths (everything anchored on `$PROJECT_DIR`), `data.table` for R
I/O, tabular outputs only.

The models that consume these matrices live in `../supervised_classification/`; the subtyping
that reads them lives in `../unsupervised_subtyping/` and the two subtype-evidence
modules beside it.

## Files

The file numbering is the run order.

| File | Role |
|---|---|
| `01_select_features.R` | the feature-selection step: takes the primary result of each association scan, keeps the FDR-significant hits, de-duplicates across scans by smallest P, and assigns each feature a deployability tier |
| `02_build_ml_ready.R` | extract → derive → split → prune → tier map → missingness → collinearity → exclusion lists. Writes the frozen matrices and the two tables everything downstream reads |
| `feature_manifest_published.csv` | reference copy of the feature set the paper reports: every feature with a plain description, its data source, its UK Biobank field and the tier it belongs to. Not read by any script |

## Run

```bash
export PROJECT_DIR=/path/to/project
export CODE_DIR=$PWD

Rscript 01_select_features.R      # scan results -> feature set + tier definitions
Rscript 02_build_ml_ready.R       # -> ml_ready/*.rds + the two manifests
```

## What `02` writes, and why there is a manifest at all

The useful shape is: read the scan results and the per-cohort matrices, write the
modelling matrices, then model. `02` does exactly that — but the matrices alone
do not say *which* of their columns each of the six tiered models may use, and
that assignment is a study decision rather than a property of the data. So `02`
writes four things and nothing else:

| Output | What reads it |
|---|---|
| `output/ml_ready/{group}_{subset}.rds` | `../supervised_classification/`, `../unsupervised_subtyping/`, Figure 4 |
| `output/feature_manifest/master_feature_manifest_final.csv` | `../supervised_classification/01`, `../supervised_classification/09`, `../unsupervised_subtyping/`, Figure 6 |
| `output/feature_manifest/tier_model_definitions_final.csv` | `../supervised_classification/01`, `../supervised_classification/02` |
| `output/feature_reports/excluded_features.csv` | `../supervised_classification/01` |

The eight stages inside `02` run once, in order, and none of them depends on a
previous run of the script: extract, derive, split, **prune**, tier map,
missingness, collinearity, the exclusion list. The prune comes before the tier map,
so each model's feature list is written once against the frozen manifest.

`ml_ready/*.rds` are immutable once `02` has written them: nothing downstream
edits them, and no feature is added to them after the fact.

## Specification (matches Methods)

- **Selection**: FDR-significant features from the proteomic, metabolomic, two
  questionnaire and disease scans, plus the declared covariates, the genetic
  principal components, the outcome and the APOE anchor. A feature selected by
  more than one scan is kept once, at its smallest P.
- **Tier models** (`tier_model_definitions_final.csv`): M1_TierA (universal
  access), M2_TierAB (adds clinical records), M3_TierABC (adds APOE e4 dose),
  M4_TierD_Olink / M5_TierD_NMR / M6_TierD_Full (molecular-assay tiers). Genetic
  principal components are excluded from every tiered model: they need
  population-level genotyping plus reference-panel projection, so a model using
  them is not deployable. An ablation recipe is printed so the exclusion can be
  tested.
- **Derived features**: the baseline-to-questionnaire interval, the ApoB/ApoA1
  ratio, four PheCode comorbidity counts, a protective-PheCode count that
  replaces the individual protective columns, and a composite smell-severity
  indicator built with the same OR-logic as the taste outcome. Each is derived
  once, here, and read from the matrices everywhere else.
- **Subsets**: four per cohort (full, Olink-assayed, NMR-assayed, both).
  Membership is decided by a sentinel feature — the most significant Olink and
  NMR feature — which is equivalent to having been assayed on that platform.

## The published manifest

`feature_manifest_published.csv` lists the 332 features of the frozen feature
set, one per row:

| Column | Meaning |
|---|---|
| `feature_id` | the column name in the modelling matrices |
| `description` | what the variable measures, in plain words |
| `variable_type` | `continuous`, `binary`, `ordinal`, `categorical`, or the `derived_` form of the first two for variables computed here rather than read from a field |
| `data_source` | the exposure family, questionnaire or assay platform it comes from |
| `ukb_field` | UK Biobank field id, PheCode, gene symbol or variant, as applicable |
| `selected_by` | the scan that selected it — ExWAS, PWAS, MWAS, DWAS, GWAS — or Covariate / Outcome / Derived. Both exposure scans are labelled ExWAS; `data_source` says which of the two |

The three variables defined in women only carry the `Reproductive (female-only)`
data source, and were tested as their own family with the sex term dropped.

Which of the six tiered models may use each feature is a separate, study-level
assignment, and it is written by `02_build_ml_ready.R` into
`tier_model_definitions_final.csv` rather than repeated here.

It carries no effect estimates, P values, counts or cluster ids: it says what
each feature is and where it comes from, not what it was found to do. Those
results are in the supplementary tables. Nothing in the pipeline reads this
file — `02_build_ml_ready.R` writes its own manifest — so it can be read on its
own without running anything.

Descriptions come from the UK Biobank data dictionary wherever the variable maps
to a single field, and otherwise from the derivation defined in this repository.
Metabolite, protein and diagnosis descriptions match the names used in the
supplementary tables.
