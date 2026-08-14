# Metabolome-wide association analysis

Cleaned, review-facing code for the MWAS reported under *Metabolome-wide
association analysis*. Each Nightingale NMR measure is tested against the taste
outcome by logistic regression. No absolute or user-specific paths, all I/O
resolved under `$PROJECT_DIR`, tabular outputs only. Subject-level outputs are
produced only where the user holds authorised access to the corresponding UK
Biobank data.

## Files

The file numbering is the run order; the two unnumbered scripts are sourced or
invoked by the numbered ones, not run by hand.

| File                             | Role                                                                                                                                                                              |
| -------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `01_nmr_qc.R`                  | `ukbnmr` (Ritchie 2023, version-3 algorithm): vendor QC flags, technical-variation removal, outlier-plate removal, extended ratios. Writes `nmr_biomarkers_qc.csv`            |
| `02_preprocess_cohorts.R`      | outlier trim (>4×IQR from the median → NA),`log1p`, then z-scoring on the pooled baseline so the three cohorts share one scale. Splits into `metabolomics_group{1,2,3}.csv` |
| `03_run_primary.R`             | the discovery-cohort scan under the pre-specified 20-covariate primary model                                                                                                      |
| `04_run_batch2_sensitivity.R`  | 18 sensitivity models: one added block at a time (oral health, metabolic, clinical completeness, effect-size brackets)                                                            |
| `05_run_batch3_validation.R`   | 3 models: held-out cross-population validation on the two other cohorts                                                                                                           |
| `06_effective_tests.R`         | correlation structure, PCA, effective number of tests, complete-linkage clustering at\|r\|>0·80. Writes `metabolite_clusters.csv`                                              |
| `07_cluster_representatives.R` | per-cluster representative selection →`cluster_representatives.csv`, which the feature manifest reads                                                                          |
| `mwas_config.R`                | sourced by`03`–`05`: paths, outcome, factor variables, `METHOD = glm`, and the `run_model()` / `summarise_batch()` helpers                                             |
| `mwas_glm_analysis.R`          | the per-feature regression engine, invoked once per model                                                                                                                         |

## Run

```bash
export PROJECT_DIR=/mnt/project          # RAP project mount, or a local checkout
export CODE_DIR=$PWD                     # so the batch scripts find the engine

Rscript 01_nmr_qc.R
Rscript 02_preprocess_cohorts.R
# copy metabolomics_group{1,2,3}.csv from output/mwas/derive/
# into input/analysis_ready/ before the batches

nohup Rscript 03_run_primary.R > mwas_primary.log 2>&1 &
nohup Rscript 04_run_batch2_sensitivity.R > mwas_batch2.log 2>&1 &
nohup Rscript 05_run_batch3_validation.R > mwas_batch3.log 2>&1 &

Rscript 06_effective_tests.R
Rscript 07_cluster_representatives.R
```

Each batch takes roughly 15–30 h on ~28 cores (327 measures × ~91k
participants), so run them detached.

## Expected inputs (under `$PROJECT_DIR/input/`)

```
input/
  metabolomics/
    nmr_merged_all.csv               # raw Nightingale export
  analysis_ready/
    phenotype_group{1,2,3}.csv       # eid, outcome, covariates per cohort
    nmr_extra_covariates.csv         # fasting hours, blood-draw hour, assessment month
    metabolomics_group{1,2,3}.csv    # from 02
```

QC outputs go to `output/metabolomics_qc/`, derived matrices to
`output/mwas/derive/`, results to `output/mwas/{primary,batch2_sensitivity,batch3_validation}/`,
interpretation to `output/mwas/effective_tests/`.

## Dependencies

`ukbnmr`, `data.table`. The scripts check for what they need and stop; they do
not install packages.

## Model specification (matches Methods)

- **Platform**: Nightingale NMR; **327** measures (absolute concentrations,
  ratios and percentages) after QC and preprocessing.
- **Values analysed**: per SD (outlier-trimmed, `log1p`, z-scored), so effects
  are odds ratios per SD and comparable across measures.
- **Engine**: standard `glm(family = binomial)`. Measures with fewer than 50
  samples or 5 cases are skipped. **No Firth fallback**: with ~91k participants
  and >3k cases the penalisation is unnecessary, and the estimator is fixed for
  every feature rather than varying with convergence behaviour.
- **Discovery cohort**: 91,042 participants with NMR and a valid outcome; the
  primary model analyses up to ~87,209 per measure (≈3,156 cases) after
  listwise deletion.
- **Primary model** (20 covariates, pre-specified): age at baseline, sex, assessment
  centre, PC1–10, baseline smoking status, alcohol-risk category, Townsend
  index, cumulative taste-affecting surgery, and three pre-analytical
  covariates (fasting hours, blood-draw hour, assessment month). BMI and the
  metabolic block are deliberately absent: they are mediators, not confounders.
- **Multiple testing**: Benjamini–Hochberg and Bonferroni across all 327
  measures; 222 FDR-significant, 151 Bonferroni-significant.
- **Effective tests and de-redundancy** (interpretation only; does not change
  the 327-measure FDR): Li & Ji effective tests = **54**; complete-linkage
  clustering at \|r\|>0·80 → **92** clusters. De-redundancy of the
  FDR-significant set gives **57** representatives, **13** of them CE-IVD.
- **CE-IVD sub-analysis**: the **37** CE-certified Nightingale biomarkers
  (Julkunen et al. 2021, *eLife* 10:e63033). FDR is applied under the full
  327-measure panel, not within the subset.

## Notes 

- The three cohorts are z-scored on one pooled baseline, so a per-SD effect
  means the same thing in the discovery and the held-out cohorts.
- The per-cluster representative rule is specified a priori on association
  fidelity and clinical deployability: for each correlation cluster, take the
  CE-IVD member whose effect direction agrees with the cluster lead
  (minimum-P), and keep the lead where no certified member agrees. It is
  independent of any downstream analysis.
