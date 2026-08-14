# Baseline exposure-wide association scan

Cleaned, review-facing code for the baseline ExWAS reported under *Exposure-wide
association across baseline assessment variables*. Adapted from the project
working scripts: no absolute or user-specific paths, all I/O resolved under
`$PROJECT_DIR`, tabular outputs only. Subject-level outputs are produced only
where the user holds authorised access to the corresponding UK Biobank data.

## Files

The file numbering is the run order; the two unnumbered scripts are sourced, not
run.

| File                      | Role                                                                                                                                                                              |
| ------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `field_list.txt`        | the 161 RAP column names this scan needs (158 UK Biobank fields), one per line: the export specification, and the list`01` checks its field registry against                    |
| `01_derive_exposures.R` | cleaning, recoding and dictionary build, then assembly of one analysis-ready exposure matrix per cohort                                                                           |
| `02_run_exwas.R`        | batch runner: 3 cohorts × 7 variants = 21 runs; writes per-run results, exclusion logs, and a consolidated text report                                                           |
| `exwas_regression.R`    | sourced. Generic hybrid glm/Firth regression engine with grouped BH FDR. Pure algorithm + I/O, no study config; shared with the follow-up ExWAS and the disease phenome-wide scan |
| `exwas_common.R`        | sourced. Study-specific configuration: outcome, covariate sets, filter thresholds, exposure families, female handling, and`$PROJECT_DIR`-based paths                            |

`01` writes the exposure matrices and dictionary to
`output/baseline_exwas/derive/`; set `EXPOSURE_DIR` to that path before running
`02`, or copy the files into `input/analysis_ready/`.

## Step 1 input: one field extract

`field_list.txt` is what the RAP Table Exporter takes: an `eid` line followed by
the exact column names, spelled out in full, because the exporter matches names
literally and has no wildcard. Only baseline-instance columns are listed: `01` keeps the
baseline instance, so later instances are not needed. Export them in one Table
Exporter run against the `participant` entity with `coding_option = REPLACE`.
`REPLACE` and not `RAW`: the derivations below compare the readable answer text
(`== "Never"`, `== "Excellent"`), so a `RAW` export would leave every one of them
unmatched without erroring. See the top-level README for the full rule.

`01` parses the field ids back out of those column names and stops if they
disagree with its internal field registry, naming the ids that appear in only
one of them, so the export specification and the fields the script expects
cannot drift apart.

## Expected inputs (under `$PROJECT_DIR/input/`)

```
input/
  raw/
    baseline_exwas_table_export.csv   # the Table Exporter output, per field_list.txt
  analysis_ready/
    phenotype_group{1,2,3}.csv      # eid, outcome, covariates per cohort
    exwas_baseline_group{1,2,3}.csv   # eid + exposure columns  (from 01)
    baseline_exwas_variable_dictionary.csv           # var_name, var_type, source, female_only, ... (from 01)
```

Results, logs, and run-time metadata are written under
`$PROJECT_DIR/output/baseline_exwas/`.

## How to run

Local or RAP (Swiss Army Knife / Cloud Workstation job):

```bash
export PROJECT_DIR=/mnt/project          # RAP project mount
export CODE_DIR=$PWD                     # so field_list.txt and the sourced
                                         # scripts are found

Rscript 01_derive_exposures.R
export EXPOSURE_DIR="$PROJECT_DIR/output/baseline_exwas/derive"
nohup Rscript 02_run_exwas.R > baseline_exwas_run.log 2>&1 &
tail -f baseline_exwas_run.log
```

`02` first builds three engine metadata files from the dictionary (main /
BMI-variant / female), subsets the female phenotypes, then loops the 21 runs. To
fan out as independent RAP jobs instead of one job, replace the `call_exwas` call
in `02_run_exwas.R` with a per-variant `dx run` submission.

## Model specification (matches Methods)

- **Outcome**: `taste_2w_strict` (binary).
- **Primary covariates**: age at baseline, sex, Townsend deprivation index, UK
  Biobank assessment centre, non-current smoking (`smoking`), and alcohol-risk
  category (`drink`). Defined once; see the Covariate definitions section.
- **Three-tier inclusion / estimation filter**:
  - L1 — drop exposures with < 100 outcome-complete cases in the cohort;
  - L2 — drop binary exposures with any 2×2 outcome-by-exposure cell < 10;
  - L3 — fit by standard ML logistic regression, with a Firth penalised-likelihood
    fallback triggered by a sparse cell (min cell < 20), non-convergence, or a
    separation warning.
- **FDR**: Benjamini-Hochberg at 5% **within each of the 11 exposure families**.
- **Counts**: 240 registered exposures; 235 enter the primary discovery-cohort
  (G1) model after filtering.
- **Sensitivity variants**: BMI-adjusted (the anthropometric family excluded to
  avoid self-adjustment); current-smell-change-adjusted (taste-specific upper
  bound); minimal age+sex model; Firth sparse-cell threshold at 5 and 15. The
  female-only group is refitted without sex. Group 3 additionally adjusts for
  self-reported ethnicity subgroup.

## Exposure families (the FDR groups)

The family name is the value stored in the variable dictionary and carried into
every results table, so the Methods text, the results tables and the figure
legends all use the same string with no code-to-name lookup anywhere.

| Family                           | Notes                                                               |
| -------------------------------- | ------------------------------------------------------------------- |
| Demographics & SES               |                                                                     |
| Lifestyle & sleep                |                                                                     |
| Dietary intake & preferences     |                                                                     |
| Anthropometric & physiological   | excluded from the BMI-adjusted variant                              |
| Oral health                      |                                                                     |
| General & mental health          |                                                                     |
| Sensory function & pain          |                                                                     |
| Clinical screening & treatment   |                                                                     |
| Blood biochemistry & haematology |                                                                     |
| Environmental exposures          | one variable, so FDR within this family is not meaningful           |
| Handedness & laterality          | three variables                                                     |
| Reproductive (female-only)       | separate family, refitted on female participants without a sex term |
