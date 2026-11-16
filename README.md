# Multiomic Associations and Subtypes of Self-reported Taste Dysfunction in UK Biobank

Analysis code for the article of the same title (Journal of Oral Rehabilitation).

## Data

The analyses use UK Biobank data under an approved application. Individual-level data cannot be
shared; access is available through the UK Biobank Access Management System
(https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access). The `field_list.txt` files
list the UK Biobank fields used by each step.

## Contents

| Directory | Analysis |
| --- | --- |
| `preprocessing/` | Case and control definitions, covariates, genotype sample QC, SARS-CoV-2 test timing |
| `exwas_baseline/` | Association scan of baseline assessment variables |
| `exwas_followup/` | Association scan of online follow-up questionnaire variables |
| `dwas/` | Association scan of prior diagnoses (three-character ICD-10 codes) |
| `pwas/` | Plasma protein (Olink) associations and pathway enrichment |
| `mwas/` | Metabolic biomarker (Nightingale NMR) associations |
| `gwas/` | GWAS (SAIGE), FUMA input, APOE genotypes, polygenic risk score, subtype GWAS (ASSET) |
| `feature_engineering/` | Feature matrix for classification and subtyping |
| `supervised_classification/` | Nested XGBoost classification models and SHAP |
| `unsupervised_subtyping/` | MOFA+ factors and k-means subtypes, by sex |
| `subtype_reproducibility/` | Cross-sex matching, re-initialization, and design sensitivity |
| `subtype_characterization/` | Molecular, diagnostic, symptom, and SARS-CoV-2 comparisons of subtypes |

## Running

Scripts are numbered in run order within each directory and are run from the repository root,
with `PROJECT_DIR` pointing to the data directory:

```bash
export PROJECT_DIR=/path/to/data
Rscript preprocessing/01_build_base_table.R
```

The association scans use the outputs of `preprocessing/`; `feature_engineering/` uses the scan
results; the classification and subtyping directories use `feature_engineering/`. In
`exwas_followup/`, scripts 01-06 run in one R session (see the header of `01`).

## Requirements

R with the packages loaded at the top of each script:

- CRAN: data.table, xgboost, pROC, PRROC, logistf, cluster, fpc, clue, mclust, dbscan, uwot,
  ukbnmr, msigdbr, stringr, optparse, broom, ggplot2, ggrepel, scales, scattermore
- Bioconductor: MOFA2 (with the Python package mofapy2), ASSET, fgsea, clusterProfiler,
  org.Hs.eg.db, AnnotationDbi

External tools: PLINK 2, SAIGE 1.5.0, PRSice-2, and bgenix.
