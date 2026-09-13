# Analysis code — UK Biobank taste dysfunction

Cohort and phenotype derivation, six association scans, feature selection,
supervised case identification, multi-omic subtyping, and the three main figures.

- Every script resolves its I/O under `$PROJECT_DIR`. Scripts that source a
  sibling also need `CODE_DIR` set to that sibling's directory.
- R with `data.table`. Shell only where an external tool does the work
  (PLINK 2, SAIGE, PRSice-2).
- Tabular and figure outputs only; participant identifiers are join keys and are
  never printed.
- `# Export and dx upload to RAP` marks an export point. A bare `fwrite` writes
  to ephemeral worker disk — upload it.

UK Biobank inputs and participant-level outputs may not be redistributed with
this code. Running it requires your own approved UK Biobank application.

## Layout and run order

```text
preprocessing/              01 base table  02 phenotypes  03 covariates
                            04 genome-wide sample QC  05 SARS-CoV-2 test timing
exwas_baseline/             01 derive exposures  02 run scan
exwas_followup/             01-06 derive by questionnaire block
                            07 assemble  08 run scan
dwas/                       01-03 PheCode matrix + diagnosis count  04 run scan
pwas/                       01 Olink QC  02 run scan  03 pathway enrichment
mwas/                       01 NMR QC  02 preprocess  03 primary
                            04 sensitivity  05 held-out  06 effective tests
                            07 cluster representatives
gwas/                       01 sample QC  02 variant QC
                            saige/ 01-04    fuma/ 01-02
                            apoe_diplotype/ 01-03    prs/ 01-06
feature_engineering/        01 select features  02 build modelling matrices
supervised_classification/  01 prepare  02 train  03 evaluate  04 calibrate
                            05 held-out  06 fairness  07 no-smell sensitivity
                            08 SHAP  09 tables
unsupervised_subtyping/     01 prepare  02 MOFA+  03 k-means  04 discriminators
                            05 tables  06 embedding
subtype_reproducibility/    01 cross-sex  02 seed re-runs  03 design  04 ARI
subtype_characterisation/   01 omics  02 diagnoses  03 severity
                            04 cross-population  05 SARS-CoV-2  06-07 HLA
figures/                    fig1/ 01-03    fig2/ a-d + assemble
                            fig3/ a-e + assemble
```

Numbered files run in order. Unnumbered `.R` files are sourced, not run.
The scans depend on `preprocessing/`; `feature_engineering/` depends on the
scans; everything after it depends on `feature_engineering/`.

Start these detached (`nohup Rscript <step> > <step>.log 2>&1 &`); a dropped
session otherwise loses the job:

| Step | Runtime |
| --- | --- |
| `mwas/03`, `mwas/04`, `mwas/05` | 15-30 h each on ~28 cores (327 measures x ~91k participants) |
| `gwas/prs/04` | hours, across the 22 chromosomes |
| `unsupervised_subtyping/02`, `subtype_reproducibility/02`, `03` | 30-90 min per sex |
| the four scan runners and `supervised_classification/02` | long enough to want a log |

`N_JOBS` caps the worker count in the scans (default 16); set it to the instance.

The five downstream directories are named after the Methods subsections:

| Methods subsection | Directory |
| --- | --- |
| Feature selection and preparation for machine learning and subtyping | `feature_engineering/` |
| Supervised classification and feature attribution | `supervised_classification/` |
| Unsupervised multi-omic subtyping | `unsupervised_subtyping/` |
| Subtype reproducibility | `subtype_reproducibility/` |
| Subtype characterisation | `subtype_characterisation/` |

### Figures

```bash
export PROJECT_DIR=/path/to/project

Rscript figures/fig1/01_collect_scan_results.R      # harmonise the six scans
Rscript figures/fig1/02_exwas_per_sd.R              # per-SD scaling table
Rscript figures/fig1/03_exposure_forest.R           # Figure 1

Rscript figures/fig2/b_apoe_forest.R                # a-d, then assemble
Rscript figures/fig2/c_protein_forest.R
Rscript figures/fig2/d_metabolite_forest.R
Rscript figures/fig2/assemble.R                     # builds a, caches it, Figure 2

Rscript figures/fig3/a_feature_sample_heatmap.R
Rscript figures/fig3/b_signature_forest.R
Rscript figures/fig3/cd_omics_heatmaps.R
Rscript figures/fig3/e_disease_heatmap.R
Rscript figures/_common/fig3_layout.R               # writes the geometry CSV
python  figures/fig3/assemble.py                    # Figure 3
```

Panel letters are added at assembly; the panel scripts do not draw them.
Appendix-figure code is not part of this release.

## Table Exporter settings

Each module carries the field list for its own extract. Two exporter inputs
matter, because getting either wrong yields a file that loads cleanly and then
makes the derivation silently wrong.

**entity** — all extracts run against `participant` and their field lists hold
`p`-prefixed field IDs. The exception is the SARS-CoV-2 test results: record-level
tables, one entity per nation (`covid19_result_england`, `_scotland`, `_wales`),
so `preprocessing/covid_field_list.txt` holds plain column names.

**coding_option** — `RAW` exports the numeric UK Biobank coding (`1`, `0`, `-1`,
`-3`); `REPLACE` exports the label the participant was shown (`Yes`, `Never`).

| Extract | Option |
| --- | --- |
| `exwas_baseline/field_list.txt`, `exwas_followup/field_list.txt` | `REPLACE` |
| everything else | `RAW` |

The two questionnaire scans hold hundreds of categorical questions, each on its
own response scale, so their derivation code compares labels. Everywhere else the
recoding rules are explicit in the script and the numeric coding is what they
expect. `preprocessing/01_build_base_table.R` accepts either and returns the
numeric coding.

The wrong option does not raise an error — comparisons simply stop matching and
the derived variable comes out empty or constant. Each derivation script
therefore reports how many values it recoded; a count of zero means a format
mismatch, not a questionnaire with no relevant answers.

## Environment

R with `data.table`. Packages, by where they come from:

- CRAN: `data.table ggplot2 ggrepel patchwork scales cowplot gridExtra ggplotify
  RColorBrewer circlize scattermore pheatmap uwot dbscan cluster fpc clue mclust
  xgboost pROC PRROC logistf broom optparse stringr arrow R.utils msigdbr
  systemfonts ukbnmr`
- Bioconductor: `MOFA2 ComplexHeatmap clusterProfiler fgsea org.Hs.eg.db
  AnnotationDbi`
- GitHub: `PheWAS/PheWAS`
- Python, for MOFA2: `mofapy2`

External tools, called from the shell steps only: PLINK 2, SAIGE 1.5.0,
PRSice-2, bgenix.

On RAP, install all of it into the image: workers have no internet, so anything
resolved at run time fails. Two consequences already handled in the code —
MOFA2's basilisk provisions a conda environment on first use, so set
`MOFA_BASILISK=FALSE` once `mofapy2` is importable; and Arial is absent from a
bare Linux worker, so install `ttf-mscorefonts-installer` or set `FIG_FONT`.

`$PROJECT_DIR` must be **writable** — every step creates directories under
`output/`. The `/mnt/project` mount is read-only, so work in a scratch directory
and link the inputs in:

```bash
mkdir -p ~/work && ln -s /mnt/project/input ~/work/input
export PROJECT_DIR=~/work
```

Worker storage is ephemeral: `# Export and dx upload to RAP` marks the points
where results have to be uploaded before the job ends.

## Expected file layout

```text
$PROJECT_DIR/
  input/
    raw/              UK Biobank field extracts, as exported from RAP
                      (incl. the SARS-CoV-2 result tables)
    eids/             cohort membership (group{1,2,3}_{full,olink}.csv)
    gwas_qc/          UK Biobank central genotype QC fields
    gwas_genotype/    PLINK --export A dosages: apoe_genotypes.raw,
                      hla_lead_genotypes.raw
    metabolomics/     nmr_merged_all.csv, the raw Nightingale extract
    fuma/             the job directory downloaded from the FUMA web service
    analysis_ready/   frozen matrices, copied here by the step that builds
                      them (see the "#> ship" lines); also
                      tastetime_assesstime.csv, supplied: eid + field 53
                      (baseline date) + field 28755 (taste questionnaire date)
    assoc_results/    each scan's primary result, renamed on the way in
    reference/        published lookup tables (PheCode map, medication list,
                      Olink panel composition, NMR annotation, variable labels)
  output/             everything the analysis produces
```

## Declared inputs and outputs

Each step declares its own I/O in a block at the top of the file, so the chain is
read from the code rather than from this file:

```text
#> in   input/analysis_ready/base_table_full.csv       consumed by this step
#> opt  input/assoc_results/pwas_{group}_primary.csv   skipped if absent
#> out  output/gwas_phenotypes/taste_gwas_phenotypes.csv   consumed downstream
#> ship output/base_table/base_table_full.csv -> input/analysis_ready/
```

`ship` is the copy that moves a frozen matrix or a scan result under `input/`,
where later stages look for it; an arrow ending in a filename also renames.
Only cross-step files are declared — a step's own reports and logs are not.
`{group}` is `group1` (discovery), `group2` and `group3` (the held-out
cross-population cohorts); other braces are likewise placeholders.

```bash
Rscript check_io.R          # every declared input has a producer
```

## Outcome

`taste_2w_strict` is a case if the participant reported a loss or change in taste
that either lasted at least two weeks **or** affected daily life — either
criterion suffices, not both. Controls reported neither a taste nor a smell
change; participants with an isolated smell change are excluded rather than
counted as controls, as they are not a clean comparison for a taste-specific
question.
