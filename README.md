# Analysis code for UKB taste-dysfunction research
Review-facing code for the UK Biobank taste-dysfunction study: cohort and phenotype derivation, association scans, and the downstream feature selection, supervised classification and unsupervised subtyping.
Conventions throughout:

- no absolute or machine-specific paths; every script resolves its I/O under
  `$PROJECT_DIR`, which you set to your RAP project;
- R with `data.table`. Shell scripts appear only where an external command-line
  tool does the work (PLINK 2, SAIGE, PRSice-2);
- tabular and figure outputs only. Participant identifiers are used as join keys
  and are never printed;
- data-export points are marked `# Export and dx upload to RAP`. On a DNAnexus
  worker a bare `fwrite` writes to ephemeral local disk, so the upload is what
  makes the result persist. Pure-computation steps need no such step.

Neither the UK Biobank inputs nor the participant-level outputs may be
redistributed with this code. Running it requires your own approved UK Biobank
application.

## Directory layout

```text
code_script/
  preprocessing/        field extract -> base table, every derived covariate,
                        the phenotype definitions and genome-wide sample QC
                        and the SARS-CoV-2 test timing flags
                        (+ field_list.txt and covid_field_list.txt, the two
                        export specifications)
  exwas_baseline/       baseline questionnaire exposure-wide scan
  exwas_followup/       follow-up questionnaire exposure-wide scan
  dwas/                 disease phenome-wide scan (PheCodes as exposures)
  pwas/                 proteome-wide scan (Olink) + pathway enrichment
  mwas/                 metabolome-wide scan (Nightingale NMR) + de-redundancy
  gwas/                 genome-wide scan (SAIGE), summary-statistic merge and
                        diagnostics, FUMA, APOE, PRS
  feature_engineering/  which variables enter the modelling matrices, and the
                        matrices themselves
  supervised_classification/
                        the tiered case-identification benchmark (XGBoost)
  unsupervised_subtyping/
                        MOFA+ factor model and the k-means subtypes
  subtype_reproducibility/
                        cross-sex concordance and seed re-runs
  subtype_characterisation/
                        proteomic, metabolomic, diagnostic, severity,
                        SARS-CoV-2, cross-population and single-locus genetic
                        characterisation of the subtypes
  figures/              main-figure plotting code (Figures 2-6)
  metadata/             covariate dictionary
```

Each directory has its own README covering inputs, run order and the reported
model specification.

The five downstream directories correspond one-to-one with the Methods
subsections, so a reader can go from a paragraph to the code that produced it:

| Methods subsection                                                   | Directory                      |
| -------------------------------------------------------------------- | ------------------------------ |
| Feature selection and preparation for machine learning and subtyping | `feature_engineering/`       |
| Supervised classification and feature attribution                    | `supervised_classification/` |
| Unsupervised multi-omic subtyping                                    | `unsupervised_subtyping/`    |
| Subtype reproducibility                                              | `subtype_reproducibility/`   |
| Subtype characterisation                                             | `subtype_characterisation/`  |

Reproducibility and characterisation are separate directories because they carry
separate claims: reproducibility asks whether the same subtypes come back, and
characterisation describes what they contain.

One part of the characterisation section is only partly included. The variant
reported per subtype was found by a subtype-stratified genome-wide scan and a
subset-based meta-analysis that run outside this package; what is here is the
per-subtype association *at that given variant*
(`subtype_characterisation/06_hla_variant_extract.sh` and `07_hla_persubtype.R`),
which is the estimate the figure shows. The search that selected it is not
reproduced.

## Exporting from RAP: entity, and RAW versus REPLACE

Every extract in the package is a Table Exporter run, and each module carries the
field list for its own. Two of the exporter's inputs are worth stating once here,
because getting either wrong produces a file that loads cleanly and then makes
the derivation code silently wrong.

### entity

All but one extract run against the **`participant`** entity, and their field
lists hold `p`-prefixed field IDs. The exception is the SARS-CoV-2 test results:
those are record-level tables with one row per test, each nation its own entity
(`covid19_result_england`, `covid19_result_scotland`, `covid19_result_wales`), so
`preprocessing/covid_field_list.txt` holds plain column names instead.

### coding_option

The option decides what a categorical answer looks like in the exported file:

- **`RAW`** — the numeric UK Biobank coding. A yes/no field arrives as `1`, `0`,
  `-1` (do not know) or `-3` (prefer not to answer).
- **`REPLACE`** — the data-coding applied, so the answer arrives as the text the
  participant was shown: `Yes`, `No`, `Never`, `Once or more daily`.

| Extract                                                              | Option      | Why                                                               |
| -------------------------------------------------------------------- | ----------- | ----------------------------------------------------------------- |
| `exwas_baseline/field_list.txt`, `exwas_followup/field_list.txt` | `REPLACE` | the derivation code matches the readable answer text              |
| everything else                                                      | `RAW`     | the derivation code does its own recoding from the numeric coding |

The two questionnaire scans contain hundreds of categorical questions, each with its own response scale. For the ExWAS analysis, we use the `REPLACE` export format, in which each response is stored as the label shown to the participant. This makes the derivation code easy to read and verify against the original questionnaire—for example, `== "Never"` or `== "Once or more daily"`. With the `RAW` format, the same comparisons would use numeric codes such as `== 0` or `== 5`. Because these codes vary across fields, they must be checked against a separate coding reference, making errors much harder to detect. For this reason, `exwas_baseline/01_derive_exposures.R` compares response labels rather than numeric codes.

For the other analyses, all recoding rules are defined explicitly in the scripts. The numeric values are therefore appropriate, so these scripts use the `RAW` export format. The script `preprocessing/01_build_base_table.R` can handle either format for its two categorical variables: it accepts both the raw numeric codes and the corresponding response labels, and always returns the numeric coding.

Importantly, using the wrong export format does not usually produce an error. Instead, the comparisons fail to match, causing the derived variables to be empty or constant. To make this problem visible, each derivation script reports the number of values successfully recoded. A count of zero therefore signals a likely export-format mismatch rather than a questionnaire scan containing no relevant responses.

## Expected file layout

```text
$PROJECT_DIR/
  input/
    raw/              UK Biobank field extracts, as exported from RAP
    eids/             cohort membership lists (group{1,2,3}_{full,olink}.csv)
    gwas_qc/          UK Biobank central genotype QC fields
    analysis_ready/   per-cohort matrices produced by the derivation steps
    assoc_results/    the primary result CSV of each scan, for the manifest step
    reference/        published lookup tables (medication code list, PheCode map)
  output/             everything the analysis produces
```

`input/analysis_ready/` holds derived rather than supplied data: the derivation
steps write there because the matrices are inputs to every later stage. The file
names are fixed, so the scans and the manifest step find them without
configuration:

```text
input/analysis_ready/
  base_table_full.csv           base table (severity fields + covariate panel)
  base_table_covar.csv          covariate-reduced base table
  proteomics_{group}.csv        QC'd Olink matrix
  metabolomics_{group}.csv      QC'd, preprocessed NMR matrix
  phecode_matrix_{group}.rds    PheCode matrix
  exwas_baseline_{group}.csv    baseline questionnaire exposures
  exwas_followup_{group}.csv    follow-up questionnaire exposures
  phenotype_{group}.csv         outcome + covariates
  apoe_{group}.csv              APOE diplotype call
  *_dictionary.csv              variable dictionaries for the two ExWAS scans
```

`group` is `group1` (discovery), `group2` and `group3` (held-out
cross-population cohorts).

## Running

```bash
export PROJECT_DIR=/path/to/project      # RAP project mount
export CODE_DIR=$PWD/<module>            # the module you are running
```

`CODE_DIR` lets a script find its siblings (a regression engine, a shared config)
without assuming a working directory. It defaults to `.`, so running a script
from inside its own directory needs no setting.

Order: `preprocessing/` first, numbered 01 to 04 — one field extract in, base
table, covariates, phenotypes and the genome-wide sample list out. Then the
per-source analysis-ready matrices (each scan's `derive_*/` subdirectory), then
the scans, then `feature_engineering/`.

## Outcome definition

The primary outcome, `taste_2w_strict`, is a case if the participant reported a
loss or change in taste that either lasted at least two weeks **or** affected
daily life — either criterion suffices, not both. Controls reported neither a
taste change nor a smell change; participants with an isolated smell change are
excluded rather than counted as controls, since they are not a clean comparison
for a taste-specific question.
