# Follow-up questionnaire exposure-wide scan

Cleaned, review-facing code for the ExWAS across seven UK Biobank online
follow-up questionnaires (2011–2021), reported under *Exposure-wide association
across online follow-up questionnaires*. Same generic engine and conventions as
the baseline scan: no absolute or user-specific paths, all I/O resolved under
`$PROJECT_DIR`, tabular outputs only. Subject-level outputs are produced only
where the user holds authorised access to the corresponding UK Biobank data.

## Files

The file numbering is the run order.

| File                             | Role                                                                                                                                                                                                       |
| -------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `field_list.txt`               | the field extract this scan needs: an`eid` line and 948 RAP column names (307 UK Biobank fields), one per line. Both the export specification and the list the derivation steps check themselves against |
| `01_derive_mental_digestive.R` | loads the export and the phenotype files into the session, defines the shared helpers and the field-list gate, derives digestive health and mental health                                                  |
| `02_derive_food_preference.R`  | taste-quality items, then PCA on the remaining food-liking items                                                                                                                                           |
| `03_derive_diet_recall.R`      | nutrient means over up to five recalls, sex-specific implausible-energy exclusion, correlation filter, supplements                                                                                         |
| `04_derive_pain.R`             | BPI interference, DN4, body-site counts, pain conditions                                                                                                                                                   |
| `05_derive_cognitive.R`        | instance selection, transforms, Trail Making B−A                                                                                                                                                          |
| `06_derive_work_environment.R` | across-job array collapse to lifetime ever-exposed indicators; closes the field-list check                                                                                                                 |
| `07_assemble_exposures.R`      | merges the per-family files into per-cohort matrices, builds the unified dictionary and case counts, flags the item-level target-leakage exclusion                                                         |
| `08_run_exwas.R`               | batch runner: 3 cohorts × 6 variants = 18 runs; writes per-run results, exclusion logs and a consolidated report                                                                                          |
| *(engine)*                     | reuses`../exwas_baseline/exwas_common.R` and `../exwas_baseline/exwas_regression.R` — the same study configuration and the same generic hybrid glm/Firth engine                                       |

## Step 1 input: one field extract

`field_list.txt` is what the RAP Table Exporter takes: an `eid` line followed by
the exact column names, spelled out in full, because the exporter matches names
literally and has no wildcard. Export them in one run against the participant
table with `coding_option = REPLACE` — the derivations compare the readable
answer text, so a `RAW` export silently matches nothing (top-level README has the
full rule) — and place the result at
`$PROJECT_DIR/input/raw/followup_exwas_merged.csv`.

How many columns a field contributes depends on how it was collected, and the
list carries exactly the ones the scripts read:

| Family                | Fields        | Columns       | Column shape                                           | Analysis variables |
| --------------------- | ------------- | ------------- | ------------------------------------------------------ | ------------------ |
| Food preferences      | 140           | 140           | `p20600` — one administration                       | 12                 |
| Diet (24-hour recall) | 63            | 315           | `p26002_i0`…`_i4` — up to five recalls, averaged | 41                 |
| Digestive health      | 15            | 15            | `p21024`                                             | 16                 |
| Mental health         | 31            | 31            | `p20514`                                             | 11                 |
| Experience of pain    | 42            | 42            | `p120019` — a single 2019 administration            | 16                 |
| Cognitive function    | 6             | 12            | `p20156_i0`/`_i1`, `p20132_i0_a0`…`_a2`       | 7                  |
| Work environment      | 10            | 393           | `p22606_a0`… — one slot per lifetime job           | 11                 |
| **Total**       | **307** | **948** |                                                        | **114**      |


## Execution model

`01` to `06` share one R session: `01` loads `raw` and `pheno_list` and the
later steps operate on those objects in memory, so they are sourced in file-number
order rather than run independently. `07` and `08` are standalone.

```bash
export PROJECT_DIR=/mnt/project          # RAP project mount, or a local checkout
export CODE_DIR=$PWD                     # so the scripts find field_list.txt

# 1. derivation — one session, in order
Rscript -e 'for (f in sort(list.files(pattern = "^0[1-6]_derive_.*[.]R$"))) source(f)'

# 2. assembly — standalone
Rscript 07_assemble_exposures.R

# 3. put the two engine inputs where the runner looks for them
cp "$PROJECT_DIR"/output/followup_exwas/derive/exwas_followup_group{1,2,3}.csv \
   "$PROJECT_DIR"/output/followup_exwas/derive/followup_exwas_variable_dictionary.csv \
   "$PROJECT_DIR"/input/analysis_ready/

# 4. the 18-run batch
nohup Rscript 08_run_exwas.R > followup_exwas_run.log 2>&1 &
```

## Expected inputs (under `$PROJECT_DIR/input/`)

```
input/
  raw/
    followup_exwas_merged.csv          # the field extract, per field_list.txt
  analysis_ready/
    phenotype_group{1,2,3}.csv         # eid, outcome, covariates per cohort
    exwas_followup_group{1,2,3}.csv    # eid + exposure columns  (from 07)
    followup_exwas_variable_dictionary.csv                       # (from 07)
```

Results, logs and run-time metadata are written under
`$PROJECT_DIR/output/followup_exwas/`.

## Model specification (matches Methods)

- **Outcome**: `taste_2w_strict` (binary).
- **Primary covariates**: age at baseline, sex, Townsend index, assessment
  centre, baseline smoking status, alcohol-risk category. BMI and current smell
  change are added back only as sensitivity variants (they lie on mediating
  paths); overall disease burden is excluded a priori and not reintroduced.
- **Engine / filters / FDR**: identical to the baseline scan (three-tier
  100/10/20 filter, hybrid glm with Firth fallback, BH within each family).
- **Variants**: `primary`, `bmi_sens`, `smell_sens`, `minimal`, and `cell5` /
  `cell15` at a permissive and a strict sparse-cell threshold. Unlike the
  baseline scan there is no female-only sub-batch and no anthropometric family
  to exclude from the BMI variant: the follow-up questionnaires contain neither
  reproductive items nor body-size measurements.
- **Counts**: 114 candidate exposures; 1 excluded as item-level target leakage
  (`gi_life_interference`, flagged in `07`); 113 enter regression.
- **FDR**: Benjamini–Hochberg at 5% within each of the seven families, never
  pooled across them.

## Exposure families and the two reporting layers

| Family                | Questionnaire (UK Biobank category)                                        | Layer                 |
| --------------------- | -------------------------------------------------------------------------- | --------------------- |
| Food preferences      | food-preference questionnaire (1039) — 7 taste probes + 5 food-liking PCs | directional           |
| Diet (24-hour recall) | Oxford WebQ — nutrients (100117), supplements (100112)                    | directional           |
| Work environment      | work history and occupational exposures (130)                              | directional           |
| Digestive health      | digestive-health questionnaire (153)                                       | comorbidity correlate |
| Mental health         | mental-health questionnaire (136) — PHQ-9, GAD-7, AUDIT                   | comorbidity correlate |
| Experience of pain    | experience-of-pain questionnaire (154) — DN4, BPI interference            | comorbidity correlate |
| Cognitive function    | cognitive-function battery (116)                                           | comorbidity correlate |

The family name is the value stored in the dictionary and carried into every
results table, so nothing downstream needs a code-to-name lookup.
