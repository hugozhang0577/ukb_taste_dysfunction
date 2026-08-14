# Preprocessing — field extract to analysis-ready variables

Everything between the raw UK Biobank export and the variables the association
scans treat as given: the base table, the phenotype definitions, the per-cohort
covariate tables, and the genome-wide sample QC. The file numbering is the run
order.

| Script | Reads | Writes |
|---|---|---|
| `01_build_base_table.R` | `input/raw/<extract>.csv` (per `field_list.txt`), `input/reference/` drug tables | `output/base_table/base_table_{full,covar}.csv` |
| `02_define_phenotypes.R` | `input/analysis_ready/base_table_full.csv` | `output/gwas_phenotypes/`, `output/smell_features/` |
| `03_merge_cohort_covariates.R` | `input/eids/` cohort lists, `input/analysis_ready/base_table_covar.csv` | `output/analysis_ready/*_covariates.csv` |
| `04_gwas_sample_qc.R` | `output/gwas_phenotypes/taste_gwas_phenotypes.csv`, `output/smell_features/smell_features.csv`, `output/base_table/base_table_covar.csv`, `input/gwas_qc/` | `output/gwas_sample_qc/` |
| `05_covid_timing_flags.R` | `input/raw/covid19_result_{england,scotland,wales}.csv` (per `covid_field_list.txt`), `input/analysis_ready/tastetime_assesstime.csv` | `output/covid/covid_temporal_flags.csv` |

`01` writes the base table under `output/`; copy or link `base_table_full.csv` and
`base_table_covar.csv` into `input/analysis_ready/` so that `02`, `03` and the
figure code find them at the fixed paths (see the top-level README for why derived
matrices live under `input/analysis_ready/`).

`04` and `05` are the steps whose input is not the participant field extract: `04` reads
the UK Biobank central genotype QC fields from `input/gwas_qc/`. That extract
must carry, besides the QC flags:

| Field | Used as |
|---|---|
| `p22001` | genetic sex; also the `sex_plink` covariate (PLINK coding 1 male / 2 female) |
| `p22006` | genetic ethnic grouping; `Caucasian` defines the White subset |
| `p22009_a1`–`p22009_a40` | genetic principal components |
| `p22010`, `p22018`, `p22019`, `p22027` | the sample-QC exclusion flags |

The genotyping batch is **not** in that extract: field 22000 is in the
participant field list, so `01` carries it through as the `batch_number` factor
and it reaches the GWAS with the other covariates.

`04` writes the SAIGE phenotype/covariate file with every column the genome-wide
covariate models name: it joins `age`, `batch_number`, `smoking`, `drink`,
`BMI`, `surg_taste_affecting_full` from the covariate table and `smell_any` from
the smell features, so the GWAS uses the same derived covariates as the panel
scans rather than re-deriving them. It stops if any of those columns is absent,
and reports per-covariate missingness rather than letting the mixed model drop
rows silently. Copy the cohort's `*_gwas_pheno.txt` into
`$COHORT_DIR/SAIGE/input/` before running `../gwas/saige/`.

## Step 1 input: one field extract

`field_list.txt` lists the 444 columns `01_build_base_table.R` needs, one field
name per line. It is both the export specification and the presence check the
script runs, so the two cannot drift apart: if a column is absent the script stops
and names it.

On DNAnexus RAP the extract is one Table Exporter run against the participant
table of the project dataset:

| Table Exporter input | value |
|---|---|
| `entity` | `participant` |
| `field_names_file_txt` | `field_list.txt` |
| `coding_option` | `RAW` (the script does its own recoding) |
| `output` | `ukb_base_fields` |

`coding_option = REPLACE` also works. The two categorical recodes
(`recode_labelled()` for fields 6150/6153/6177, `yes_no()` for field 2443) accept
the numeric UK Biobank coding or the replaced value label and return the numeric
coding either way. Likewise, a header that carries the entity prefix
(`participant.p31`) is normalised to the bare field name on read, so either header
style produces the same base table.

The extract is wide — the procedure-date and medication-code array blocks are 318
of the 444 columns — but neither output carries them: each block is reduced to its
covariates and the source columns are dropped immediately afterwards.

Place the export under `input/raw/`. Under the UK Biobank data-access agreement,
neither the input extract nor the outputs may be redistributed with this code.
Participant identifiers are used only as join keys and are never printed.

## Step 5 input: the SARS-CoV-2 test records

The test results are **not participant fields**, and this is the only extract in
the package that is not. Every other field list in every other module names
`p`-prefixed participant field IDs and is exported against the `participant`
entity. The COVID results are record-level tables — one row per test, a
participant appearing once per test or not at all — and each nation is its own
entity. `covid_field_list.txt` accordingly holds plain column names, not field
IDs; it is not malformed.

Three Table Exporter runs, one per nation, all with the same field list:

| Table Exporter input | value |
|---|---|
| `entity` | `covid19_result_england`, then `covid19_result_scotland`, then `covid19_result_wales` |
| `field_names_file_txt` | `covid_field_list.txt` — `eid`, `specdate`, `result` |
| `coding_option` | `RAW` (`result` is wanted as 0/1, and `05` checks that it is) |
| `output` | `covid19_result_<nation>` |

Export all three columns **in one run per nation**. Exporting the dates
separately and pairing them with the results by row order would attach each date
to whichever test happened to sit on the same line, and `eid` cannot fix it
afterwards because a participant has many tests. `05` stops if a file is missing
any of the three columns, for exactly this reason.

`05` also needs `input/analysis_ready/tastetime_assesstime.csv` for the
questionnaire completion date (field 28755) — the same small extract
`../dwas/02_build_phecode_matrix.R` uses. Every flag it writes is defined
relative to that date, because a positive test recorded after the questionnaire
cannot explain an answer given before it.

Its output goes to `output/covid/`; copy it into `input/analysis_ready/` for
`../subtype_characterisation/05_covid_association.R`.

## Step 1 second input: two published lookup tables

The medication covariate needs two tables under `input/reference/`. Both are
drug-level, contain no participant data, and are the material behind the
manuscript's supplementary medication table:

| File | Columns | Role |
|---|---|---|
| `taste_drugs_full_list.csv` | `category`, `drug_name`, `ukb_code`, `ukb_med_name` | UK Biobank medication codes for every drug in the literature review |
| `taste_affecting_drugs_literature.csv` | `drug_name`, `category`, `high_risk` | which of those drugs have a reported dysgeusia incidence of at least 1% |

The covariate uses the high-risk subset: the `ukb_code` values of the drugs flagged
`high_risk`, matched on `drug_name`. Both medication variables use that same code
set and differ only in which instances they scan (all instances vs instance 0).

## Field groups (step 1)

| Group | Fields | Derived columns |
|---|---|---|
| Taste / smell severity | p28612–p28617 | `taste_change`, `taste_time`, `taste_extent`, `smell_change`, `smell_time`, `smell_extent` |
| Demographics | p31, p34, p21003 | `sex`, `age`, `age_baseline` |
| Genotyping batch | p22000 | `batch_number` (factor; the GWAS batch covariate) |
| Anthropometry / spirometry | p48, p49, p50, p74, p21001, p20150, p20151, p20258 | `Waist_circumference`, `Hip_circumference`, `Standing_height`, `Fasting_time`, `BMI`, `FEV1`, `FVC`, `FEV1_FVC_ratio` |
| Pulse and blood pressure | p102, p4079, p4080 (arrays a0, a1) | `paulse`, `DBP`, `SBP` (mean of the two readings) |
| Blood count, urine, biochemistry | p30000–p30890, baseline instance | 55 named biomarkers, values unchanged |
| Smoking | p1239, p1249, p2644 | `smoking_bin`, `smoking_unknown`, `smoking`, `smoking_sens` |
| Alcohol | p1568, p1578, p1588, p1598, p1608, p5364 | `alcohol_unit`, `drink` |
| Hypertension | p6150 (instances 0–2), p6153, p6177, plus measured BP | `n6150_i0-i2`, `n6153_i0`, `n6177_i0`, `hypertension` |
| Diabetes | p2443 (instances 0–3) | `diabetes` |
| Prevalent cancer | p40008 (all instances) | `cancer` |
| Oral health | p6149 (instances 0–3) | 6 problems × baseline/cumulative, `periodontal_indicator`, `oral_problem_count`, `any_oral_problem` (each also `_baseline`) |
| Oral / dental surgery | p41272, p41282 (arrays), p53 | 15 procedure groups × `_full`/`_baseline`, plus `surg_taste_affecting_*` and `surg_non_taste_affecting_*` |
| Taste-affecting medication | p20003 (4 instances × 48 arrays) | `taste_drug_user_numeric`, `taste_drug_user_baseline_numeric` |

Taste / smell coding: p28615 is 1 yes / 0 no / -1 do not know / -3 prefer not to
answer, and the base table is restricted to the definite 1/0 answers. p28616 and
p28613 are duration bands 0–3; p28617 and p28614 are daily-life impact 0/1. The
five fields other than p28615 keep their source coding in the base table, negative
codes included; `02_define_phenotypes.R` decides how to treat them.

Smoking is coded 1 = non-current (never, or smoked in the past and no longer does),
0 = current. `smoking` is the primary covariate and groups the participants who
cannot be classified from these three fields with the current smokers;
`smoking_sens` is the sensitivity variant and groups them with the non-current
smokers instead. `smoking_bin` keeps them as missing.

Note that the descriptive baseline table counts an unclassifiable participant as
non-current, whereas the `smoking` covariate counts them as current, so the two
figures are not interchangeable.

`drink` is 1 for high-risk drinking, using a sex-specific weekly-unit threshold on
`alcohol_unit`.

Field 6149 is a multiple-response field: an instance holds the selected codes
pipe-separated (`4|5|6`), with `-7` for none of the above. The periodontal
indicator is bleeding gums, painful gums **or** loose teeth — any one of the three,
not a count of them.

Field 41272 holds a participant's OPCS-4 procedure codes as one pipe-separated
list, and `p41282_a0…a125` hold the corresponding dates positionally aligned to
that list. The baseline assessment date (p53) splits lifetime exposure (`_full`)
from pre-baseline exposure (`_baseline`). Procedures are grouped by whether the
operative field plausibly involves the taste pathway (lip, tongue, palate, jaw,
tonsil, salivary gland) or not (dental and other oral procedures).

## Step 2: phenotype definitions

The primary outcome, `taste_2w_strict`, is a case if the participant reported a
loss or change in taste that either lasted at least two weeks **or** affected daily
life — either criterion suffices, not both. Controls reported neither a taste
change nor a smell change; participants with an isolated smell change are excluded
rather than counted as controls. The four-week and relaxed-control variants are
sensitivity definitions and never redefine the primary outcome.

The same script derives the three machine-learning smell features (`smell_any`,
`smell_time`, `smell_extent`). The taste and smell sections work on separate copies
of the input, because they treat a negative `smell_change` code differently: the
taste definitions read it as "no smell change", so the participant can still serve
as a control, whereas the smell features read it as missing, so an unknown smell
status does not become a `smell_any` of 0.

## Running

```bash
export PROJECT_DIR=/path/to/project      # or the RAP project mount
export CODE_DIR=$PWD                     # so field_list.txt is found
export UKB_EXTRACT=ukb_base_fields.csv   # optional, this is the default

Rscript 01_build_base_table.R
# copy output/base_table/base_table_{full,covar}.csv -> input/analysis_ready/
Rscript 02_define_phenotypes.R
Rscript 03_merge_cohort_covariates.R
Rscript 04_gwas_sample_qc.R
```

No absolute machine paths are hard-coded.
