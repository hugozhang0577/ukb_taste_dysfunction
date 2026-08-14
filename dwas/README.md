# Disease phenome-wide association scan

Cleaned, review-facing code for the disease-wide association scan reported under *Disease
phenome-wide association scan of prior ICD-coded diagnoses*. Each PheCode is
tested as a binary exposure against the taste outcome. Same conventions as the
exposure-wide scans: no absolute or user-specific paths, all I/O resolved under
`$PROJECT_DIR`, tabular outputs only. Subject-level outputs are produced only
where the user holds authorised access to the corresponding UK Biobank data.

## Files

The file numbering is the run order; the two unnumbered scripts are invoked by
`04`, not run by hand.

| File                            | Role                                                                                                                                                                                         |
| ------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `field_list.txt`              | the field extract this scan needs: an`eid` line and 2262 RAP column names, one per line. Both the export specification and the list `01` checks the First-Occurrences dictionary against |
| `01_clean_first_occurrence.R` | merges the export batches into a long table; cleans the UK Biobank sentinel dates; decodes the report source (objective vs self-report)                                                      |
| `02_build_phecode_matrix.R`   | keeps diagnoses dated before the taste assessment; maps ICD-10 to PheCodes (PheCode Map v1.2); writes one case/control matrix per cohort plus per-PheCode summaries                          |
| `03_build_diagnostic_count.R` | `total_dx_count` — distinct ICD-10 codes per participant — merged into the phenotype files as a primary-model covariate                                                                  |
| `04_run_dwas.R`               | batch runner: five covariate models on the discovery cohort, one on each held-out cohort                                                                                                     |
| `dwas_regression.R`           | invoked by`04`, once per cohort × model. Merges the PheCode matrix and the phenotype at run time, fits glm with a Firth fallback, applies BH and Bonferroni across the analysable set     |
| `summarise_results.R`         | invoked by`04` once the models finish; writes the per-model comparison table and the labelled FDR-significant table that Figure 3 uses for PheCode labels                                  |

## Step 1 input: one field extract

`field_list.txt` is what the RAP Table Exporter takes: an `eid` line followed by
the exact column names, spelled out in full, because the exporter matches names
literally and has no wildcard. The First Occurrences fields (Category 1712) are
single-instance and unarrayed, so each is one column and the list is exactly
2262 names — a `Date <code> first reported` field and a matching
`Source of report of <code>` field for each of 1131 ICD-10 three-character
codes. Export against the participant table; both `RAW` and `REPLACE` work,
since `01` decodes the source values itself.

`01` rebuilds the date/source pairing from `FO_dic.csv`, the First-Occurrences
field dictionary, and stops if that pairing and `field_list.txt` disagree,
naming the fields that appear in only one of them — so the export specification
and the fields the script expects cannot drift apart. Fields listed but absent
from the exported batches are reported as a shortfall rather than dropped
silently.

The batches may be split across several CSVs (`01` merges them by `eid`); this
is only a convenience for exporting ~2000 columns.

## Step 2 second input: the two date anchors

`02` needs one more small extract, `tastetime_assesstime.csv`, with three
columns: `eid`, `p53_i0` (baseline assessment date) and `p28755`
(taste-assessment date). The temporal filter keeps only diagnoses first recorded
before the taste assessment, so every PheCode in the scan is a **prior**
diagnosis. Records whose date the UK Biobank flags as unreliable are kept with
the date set to missing rather than discarded, and are retained by the filter.

## Expected inputs (under `$PROJECT_DIR/input/`)

```
input/
  raw/
    fo_batches/*.csv                 # the field extract, per field_list.txt
    FO_dic.csv                       # First-Occurrences field dictionary
  analysis_ready/
    tastetime_assesstime.csv         # eid, p53_i0, p28755
    phenotype_group{1,2,3}.csv       # eid, outcome, covariates per cohort
    phecode_matrix_group{1,2,3}.rds  # from 02
```

`01` and `02` write under `$PROJECT_DIR/output/dwas/derive/`; `03`
writes `total_dx_count` back into the phenotype files in place (dropping and
recomputing any existing column, so re-running is a no-op). Results and logs go
to `$PROJECT_DIR/output/dwas/results/`.

## How to run

Local or RAP (Swiss Army Knife / Cloud Workstation job):

```bash
export PROJECT_DIR=/mnt/project          # RAP project mount, or a local checkout
export CODE_DIR=$PWD                     # so the scripts find field_list.txt

Rscript 01_clean_first_occurrence.R
Rscript 02_build_phecode_matrix.R
Rscript 03_build_diagnostic_count.R

# 04 reads the matrices from input/analysis_ready/; either copy them there or
# point PHECODE_DIR at 02's output directory
export PHECODE_DIR="$PROJECT_DIR/output/dwas/derive"
nohup Rscript 04_run_dwas.R > dwas_run.log 2>&1 &
```

`N_JOBS` sets the number of cores the engine uses per model (default 10).

## Dependencies

`data.table`, `optparse`, `logistf` (the Firth fallback), and `PheWAS` — the
last supplies `createPhenotypes()` and PheCode Map v1.2, and is not on CRAN.

## Model specification (matches Methods)

- **Outcome**: `taste_2w_strict` (binary); each PheCode is the exposure.
- **Primary covariates**: age at baseline, sex, assessment centre, Townsend
  index, baseline smoking status, alcohol-risk category, BMI, and the total
  count of distinct ICD-coded diagnoses (`total_dx_count`).
- **Engine**: glm with a Firth penalised-likelihood fallback on
  convergence or separation warnings.
- **Inclusion**: PheCodes with at least 100 cases in the discovery cohort; at
  least 50 in the held-out cohorts, which are far smaller — the discovery
  threshold would leave almost nothing testable there.
- **Multiple testing**: Benjamini–Hochberg **and** Bonferroni across the full
  analysable PheCode set (454 PheCodes in the discovery cohort).
- **Sensitivity models** (discovery cohort): `primary+smell` (taste-specific
  upper bound); `no_bmi`; `no_centre` (drops assessment centre and the
  diagnostic-count term together); `minimal` (age, sex and the diagnostic-count
  term).

## Output naming

Result files are named `group{1,2,3}_<model>.csv`, as produced:
`group1_primary.csv`, `group1_primary+smell.csv`, `group1_no_bmi.csv`,
`group1_no_centre.csv`, `group1_minimal.csv`, `group2_primary.csv`,
`group3_primary.csv`. The two held-out files are named `_primary` but are run
under the smell-adjusted covariate set — see the next section.

## Reported tables

`summarise_results.R` writes two tables into the results directory:

| File                               | Contents                                                                                                                                                                                                                        |
| ---------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `dwas_models_comparison.csv`     | one row per cohort × model: hit counts at each threshold, the top PheCode, and the inflation diagnostic                                                                                                                        |
| `dwas_fdr_significant_table.csv` | the FDR-significant PheCodes of the**smell-adjusted** discovery model, with description, clinical category, case and control counts, OR with 95% CI, p, FDR and direction. Figure 3 reads it as a PheCode-to-label lookup |


## Three-layer PheCode annotation (reporting)

PheCode labels used in figures and tables come from the PheWAS catalogue, the
UK Biobank First-Occurrences ICD-10 names, and a cohort-derived ICD-10 trace, in
that order, with unresolved codes placed in an "Other" residual group. That
annotation is applied at the reporting stage, not here: this package writes the
raw PheCode identifiers and the association statistics.
