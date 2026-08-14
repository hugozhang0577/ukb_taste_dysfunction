# Supervised classification and feature attribution

Code for the Methods subsection *Supervised classification and feature
attribution*: six nested models ordered by how hard their inputs are to obtain
in practice, trained on the
discovery cohort and read out on the two held-out cross-population cohorts.
Same conventions as the other `code_script/` modules: no absolute or
user-specific paths (everything anchored on `$PROJECT_DIR`), `data.table` for R
I/O, tabular and figure outputs only.

The frozen matrices, the feature manifest and the tier definitions come from
`../feature_engineering/`; nothing here modifies them.

## Files

The file numbering is the run order.

| File                          | Role                                                                                                                         |
| ----------------------------- | ---------------------------------------------------------------------------------------------------------------------------- |
| `01_prepare_features.R`     | per-model feature recipe from the tier definitions and the exclusion list                                                    |
| `02_train_xgboost.R`        | XGBoost on the discovery cohort; five-fold stratified cross-validated grid search; out-of-fold predictions saved with`eid` |
| `03_evaluate.R`             | ROC-AUC with DeLong intervals, PR-AUC, per-model summary                                                                     |
| `04_calibrate.R`            | fits the probability calibrators the held-out and fairness steps apply. Not reported: see the note below                     |
| `05_heldout_validation.R`   | the two held-out cross-population cohorts                                                                                    |
| `06_fairness_audit.R`       | subgroup performance gaps                                                                                                    |
| `07_sensitivity_no_smell.R` | refit without the smell features                                                                                             |
| `08_build_tables.R`         | the publication tables                                                                                                       |

## Run

```bash
export PROJECT_DIR=/path/to/project
export CODE_DIR=$PWD

for s in 01 02 03 04 05 06 07 08; do Rscript ${s}_*.R; done
```

`02` is the long step (six models × a cross-validated grid search); run it
detached.

## Expected inputs (all from `../feature_engineering/02_build_ml_ready.R`)

```
output/
  ml_ready/{group}_{subset}.rds
  feature_manifest/master_feature_manifest_final.csv
  feature_manifest/tier_model_definitions_final.csv
  feature_reports/excluded_features.csv
```

Results go to `output/models/xgboost/` and `output/model_reports/`.

## Specification (matches Methods)

- **Tier models**: M1_TierA (universal access: questionnaire and demographics,
  no genotyping), M2_TierAB (adds clinical records), M3_TierABC (adds APOE e4
  dose), M4_TierD_Olink / M5_TierD_NMR / M6_TierD_Full (molecular-assay tiers).
  The models are nested, so the increment from one tier to the next is what the
  extra data buys.
- **Engine**: XGBoost, fitted on the discovery cohort only. Five-fold
  stratified cross-validated grid search; out-of-fold predictions are what the
  discovery-cohort metrics are computed from. The held-out cohorts are never
  used for selection or tuning.
- **Missing values**: kept as NA. The boosted trees learn a default split
  direction per feature, so missingness is information rather than a problem.
  The only features dropped are the structurally unusable ones — missing in
  almost every participant of the discovery cohort — listed in
  `excluded_features.csv`.
- **Reported outputs**: ROC-AUC with DeLong confidence intervals, PR-AUC,
  a subgroup fairness audit, and a no-smell sensitivity analysis. Calibration
  quality and decision-curve net benefit are not reported.
