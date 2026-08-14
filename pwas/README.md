# Proteome-wide association analysis

Cleaned, review-facing code for the PWAS reported under *Proteome-wide
association analysis*. Each Olink protein is tested against the taste outcome by
Firth penalised-likelihood logistic regression. No absolute or user-specific
paths, all I/O resolved under `$PROJECT_DIR`, tabular outputs only.
Subject-level outputs are produced only where the user holds authorised access
to the corresponding UK Biobank data.

## Files

The file numbering is the run order; the unnumbered script is invoked by
`02`, not run by hand.

| File                      | Role                                                                                                                                                                                                  |
| ------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `01_olink_qc.R`         | QC of the Olink long-format export: filter to the cohort, drop proteins >30% and samples >50% missing, separate below-LOD from true missing, stratified imputation. Writes`proteomics_<cohort>.csv` |
| `02_run_pwas.R`         | runs the pre-specified primary model over the QC'd matrix                                                                                                                                             |
| `03_gsea.R`             | pre-ranked pathway enrichment (fgsea) on the primary result, restricted to the assayed panel as background                                                                                            |
| `pwas_firth_analysis.R` | invoked by`02`. Firth regression per protein, parallelised, profile-likelihood CIs, BH and Bonferroni                                                                                               |

## Run

```bash
export PROJECT_DIR=/mnt/project          # RAP project mount, or a local checkout
export CODE_DIR=$PWD

# QC, one cohort per run
COHORT=group1 Rscript 01_olink_qc.R
COHORT=group2 Rscript 01_olink_qc.R
COHORT=group3 Rscript 01_olink_qc.R

# 02 reads the QC'd matrix from input/analysis_ready/; either copy it there or
# point OLINK_DIR at 01's output directory
export OLINK_DIR="$PROJECT_DIR/output/pwas/olink_qc"
nohup Rscript 02_run_pwas.R > pwas_run.log 2>&1 &

Rscript 03_gsea.R
```

`N_JOBS` sets the number of cores per model (default 16).

## Expected inputs (under `$PROJECT_DIR/input/`)

```
input/
  proteomics/
    olink_full_dataset.csv        # long-format Olink export (eid, protein_id, result, ins_index, Batch)
  analysis_ready/
    phenotype_group{1,2,3}.csv    # eid, outcome, covariates per cohort
    extra_covariates.csv          # Townsend, joined by eid
    proteomics_group1.csv         # from 01 (or set OLINK_DIR)
  assoc_results/
    pwas_primary.csv              # the primary result, read by 03
```

QC outputs go to `$PROJECT_DIR/output/pwas/olink_qc/`, association results to
`output/pwas/results/`, enrichment to `output/pwas_gsea/panel_restricted/`.

## Dependencies

`data.table`, `optparse`, `logistf`, `parallel`, `ggplot2`, `gridExtra`; `03`
additionally needs `clusterProfiler`, `org.Hs.eg.db`, `msigdbr` and `fgsea`.
The scripts check for what they need and stop; they do not install packages.

## Model specification (matches Methods)

- **Platform**: Olink Explore 3072; 2,920 proteins after QC (NPX, baseline
  instance).
- **Engine**: Firth (`logistf`, maxit 250) for every protein, not as a
  fallback: with ~17k assayed participants and few cases per protein, maximum
  likelihood is biased and can fail to converge. Proteins with fewer than 50
  samples or 5 cases are skipped.
- **Primary covariates** (17, specified in advance from the study design):
  age at baseline, sex, Olink batch, PC1–10, baseline smoking status,
  alcohol-risk category, cumulative taste-affecting surgery, Townsend index.
  BMI and the wider metabolic block are deliberately absent: they lie on the
  causal path between the protein panel and the outcome rather than confounding
  it, so adjusting for them would remove signal by construction.
- **Cross-omics alignment**: the metabolomic primary model is this set with the
  assay-specific terms swapped — no Olink batch (the NMR assay has no
  equivalent), plus assessment centre and three pre-analytical covariates — so
  the two omics results are directly comparable.
- **Discovery cohort**: median 15,110 participants (561 cases) per protein.
- **Multiple testing**: Benjamini–Hochberg and Bonferroni across all 2,920
  proteins; 17 FDR-significant, 4 Bonferroni-significant.
- **Held-out validation**: the same engine on the two held-out
  cross-population cohorts, with their own phenotype inputs.

## Notes

- QC runs one cohort at a time. Imputation is fitted within the cohort being
  processed, so the cohorts are deliberately not pooled; every output file
  carries the cohort in its name.
- Below-LOD values are imputed with MinProb (left-censored, missing not at
  random) and true missing values with the protein-wise median, rather than
  treating the two as one missingness mechanism.
- `03_gsea.R` uses the assayed panel as the enrichment background rather than
  the whole proteome, so the result reflects enrichment within what was
  measured. Its output table supplies the pathway panel of Figure 2.
