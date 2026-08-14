# Subtype characterisation

Code for the Methods subsection *Subtype characterisation*: what the four
subtypes contain, across the full proteomic and metabolomic panels, the disease
phenome, chemosensory severity, and the two held-out cross-population cohorts.
Same conventions as the other `code_script/` modules: no absolute or
user-specific paths (everything anchored on `$PROJECT_DIR`), `data.table` for R
I/O, tabular and figure outputs only.

Whether the subtypes reproduce is `../subtype_reproducibility/`; the fit itself
is `../unsupervised_subtyping/`.

## Files

Each script stands alone; the numbering is the order they appear in Methods.

| File                                | Role                                                                                                                                                                                                                               |
| ----------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `01_fullpanel_omics.R`            | full Olink and NMR panels, each subtype against the other three pooled; the effect is the subtype term of`biomarker ~ is_subtype + age_baseline + sex + centre`, BH-FDR within subtype and platform. Writes what Figure 5B draws |
| `02_phecode_enrichment.R`         | per-subtype enrichment across the full PheCode phenome, age- and sex-adjusted (glm with a Firth fallback), BH-FDR within subtype. The PheCodes used as MOFA inputs are removed from the discovery universe                         |
| `03_severity_orthogonality.R`     | whether the subtypes are orthogonal to chemosensory severity rather than gradations of one severity axis                                                                                                                           |
| `04_crosspopulation_projection.R` | projects the two held-out cohorts into the discovery factor space; nearest-centroid assignment plus an independent k-means concordance. Takes an`M`/`F` argument                                                               |
| `05_covid_association.R`          | pre-assessment SARS-CoV-2 positivity, one subtype versus the rest, crude and age-adjusted (and a spline-age variant). Writes what Figure 5E draws on the left                                                                      |
| `06_hla_variant_extract.sh`       | one Swiss Army Knife job: PLINK2 exports additive dosages at the reported lead variant from the QC-filtered chromosome 6 bgen                                                                                                      |
| `07_hla_persubtype.R`             | that variant's effect in each subtype's cases against the shared controls, with a heterogeneity test across the four. Writes what Figure 5D draws                                                                                  |

## Run

```bash
export PROJECT_DIR=/path/to/project
export CODE_DIR=$PWD/../unsupervised_subtyping    # holds _subtype_map.R

Rscript 01_fullpanel_omics.R
Rscript 02_phecode_enrichment.R
Rscript 03_severity_orthogonality.R
Rscript 04_crosspopulation_projection.R M
Rscript 04_crosspopulation_projection.R F
Rscript 05_covid_association.R
bash    06_hla_variant_extract.sh
Rscript 07_hla_persubtype.R
```
