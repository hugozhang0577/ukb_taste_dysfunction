# Main-figure code (Figures 2-6)

Plotting code for the five main figures. Supplementary-figure code is not part of
this release.

Figure 1 is a study schematic drawn by hand and has no code.

## Layout

```text
figures/
  fig2_multiomics_assoc.R    Fig 2, six panels in one script
  fig3/                      Fig 3, two preparation steps + one figure script
    tier_flagged_results.R
    exwas_continuous_sd.R
    exwas_dwas.R
  fig4/                      Fig 4, three panel scripts + assembler
    panelA_feature_sample_heatmap.R
    panelBC_umap.R
    panelD_signature_forest.R
    assemble.R
  fig5/                      Fig 5, four panel scripts + assembler
    panelA_reproducibility_radial.R
    panelB_omics_heatmap.R
    panelC_disease_heatmap.R
    panelDE_genetics_covid.R
    assemble.R
  fig6/                      Fig 6, attribution step + six panels
    shap_prep.R
    ml_panel.R
```

## Running

```bash
export PROJECT_DIR=/path/to/project

Rscript figures/fig2_multiomics_assoc.R
Rscript figures/fig3/tier_flagged_results.R
Rscript figures/fig3/exwas_continuous_sd.R
Rscript figures/fig3/exwas_dwas.R

CODE_DIR=figures/fig4 Rscript figures/fig4/panelA_feature_sample_heatmap.R
CODE_DIR=figures/fig4 Rscript figures/fig4/panelBC_umap.R
CODE_DIR=figures/fig4 Rscript figures/fig4/panelD_signature_forest.R
CODE_DIR=figures/fig4 Rscript figures/fig4/assemble.R

CODE_DIR=figures/fig5 Rscript figures/fig5/panelA_reproducibility_radial.R
CODE_DIR=figures/fig5 Rscript figures/fig5/panelB_omics_heatmap.R
CODE_DIR=figures/fig5 Rscript figures/fig5/panelC_disease_heatmap.R
CODE_DIR=figures/fig5 Rscript figures/fig5/panelDE_genetics_covid.R
CODE_DIR=figures/fig5 Rscript figures/fig5/assemble.R

Rscript figures/fig6/shap_prep.R
Rscript figures/fig6/ml_panel.R
```

## What each figure shows

| Figure      | Content                                                                                                                                                                                         |
| ----------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **2** | Multi-omic association overview: GWAS Manhattan and QQ, APOE diplotype forest, Olink volcano, pathway enrichment, and the non-redundant NMR panel across the discovery and two held-out cohorts |
| **3** | Exposure-wide and disease-wide scans: paired volcano and signal-count panels per exposure family                                                                                                |
| **4** | Subtype structure: feature-by-participant heatmap, the per-sex factor embedding, and the per-subtype signature forest                                                                           |
| **5** | Subtype evidence: reproducibility across re-runs and design choices, molecular signature, disease enrichment, and the genetic and COVID panels                                                  |
| **6** | Prediction: ROC and cross-cohort AUC for the six tier models, with and without the smell features, plus SHAP importance for the universal-access tier                                           |

## Input provenance

Produced by the analysis code in this repository:

| Input                                                                                    | Produced by                                                                                                        |
| ---------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------ |
| `output/feature_manifest/master_feature_manifest_final.csv`                            | `feature_engineering/02_build_ml_ready.R`                                                                        |
| `output/ml_ready/group1_full.rds`                                                      | same                                                                                                               |
| `output/models/xgboost/*_oof_preds.rds`, `*_booster.rds`                             | `supervised_classification/02_train_xgboost.R`                                                                   |
| `output/model_reports/eval/xgb_eval_summary.csv`                                       | `supervised_classification/03_evaluate.R`                                                                        |
| `output/model_reports/eval/heldout_val_summary.csv`                                    | `supervised_classification/05_heldout_validation.R`                                                              |
| `output/model_reports/eval/sensitivity_no_smell_summary.csv`                           | `supervised_classification/07_sensitivity_no_smell.R`                                                            |
| `output/model_reports/eval/shap_*_mean_abs.csv`                                        | `figures/fig6/shap_prep.R`                                                                                       |
| `output/subtyping/inputs/clustering_input_g1_{m,f}.rds`                                | `unsupervised_subtyping/01_prepare_clustering_features.R`                                                        |
| `output/subtyping/mofa/mofa_loadings_g1_{m,f}.rds`                                     | `unsupervised_subtyping/02_mofa_fit.R`                                                                           |
| `output/subtyping/clusters/cluster_assignments_g1_{m,f}_k4.rds`                        | `unsupervised_subtyping/03_cluster_on_factors.R`                                                                 |
| `output/subtyping/reports/cluster_profiles_{m,f}_k4.csv`                               | `unsupervised_subtyping/04_cluster_discriminators.R`                                                             |
| `output/subtyping/evidence_disease/phecode_enrichment_g1_by_subtype.csv`               | `subtype_characterisation/02_phecode_enrichment.R`                                                               |
| `output/subtyping/evidence_omics/omics_signature_{nmr,olink}_plotdata.csv`             | `subtype_characterisation/01_fullpanel_omics.R`                                                                  |
| `output/subtyping/evidence_covid/covid_onevsrest_OR.csv`, `covid_subtype_labels.csv` | `subtype_characterisation/05_covid_association.R`                                                                |
| `input/analysis_ready/covid_temporal_flags.csv` | `preprocessing/05_covid_timing_flags.R` |
| `output/subtyping/evidence_genetic/hla_persubtype_effects.csv`                         | `subtype_characterisation/06_hla_variant_extract.sh` then `07_hla_persubtype.R`                                |
| `output/subtyping/reports/fig_tables/umap_coords_{m,f}_k4.csv`                         | `unsupervised_subtyping/06_umap_embedding.R`                                                                     |
| `output/evidence_tiering/tier_flagged_results.csv`                                     | `figures/fig3/tier_flagged_results.R`                                                                     |
| `output/evidence_tiering/exwas_continuous_SD.csv`                                      | `figures/fig3/exwas_continuous_sd.R`                                                                      |
| `output/dwas/results/dwas_fdr_significant_table.csv`                                 | `dwas/summarise_results.R`                                                                                       |
| `output/pwas_gsea/panel_restricted/tables/GSEA_all_significant.csv`                    | `pwas/03_gsea.R`                                                                                                 |
| `output/subtyping/reports/fig4_row_order.csv`                                          | `fig4/panelA_feature_sample_heatmap.R` (written, then reused by supplementary figures to keep row order aligned) |
| `input/analysis_ready/base_table_full.csv`                                             | `preprocessing/01_build_base_table.R`                                                                            |

The primary result of each association scan, renamed into `input/assoc_results/`
as the manifest step also expects:

| Input                                                            | Scan and model                                                                                                    |
| ---------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------- |
| `gwas_primary_sumstats.tsv.gz`, `gwas_primary_lead_snps.csv` | `gwas/saige/03_merge_sumstats.sh` (model M1) and the top-hit table from `gwas/saige/04_manhattan_qq_lambda.R` |
| `apoe_diplotype_association.csv`                               | `gwas/apoe_diplotype/02_apoe_association.R`                                                                     |
| `pwas_primary.csv`                                             | `pwas/02_run_pwas.R` (the primary model)                                                                        |
| `mwas_primary.csv`                                             | `mwas/03_run_primary.R`, the primary model                                                                      |
| `mwas_group{2,3}_primary.csv`                                  | `mwas/05_run_batch3_validation.R`, V1 and V3                                                                    |

Not produced by this repository; supply them alongside the inputs above:

| Input                                                                            | What it is                                                                                                                                                                                                                                   |
| -------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `input/reference/phecode_map_v1_2_icd10.csv`, `phecode_definitions_v1_2.csv` | PheCode Map v1.2, published by its authors                                                                                                                                                                                                   |
| `input/reference/olink_panel_mapping.csv`                                      | Olink Explore panel composition, from the assay documentation                                                                                                                                                                                |
| `input/reference/nmr_label_map.csv`                                            | Nightingale biomarker code to published abbreviation                                                                                                                                                                                         |
| `input/reference/exwas_{baseline,followup}_labels.csv`                         | variable code to display label for the questionnaire exposures                                                                                                                                                                               |

| `input/reference/nmr_metabolite_annotation.csv`                                | Nightingale biomarker annotation used by the evidence-tiering step                                                                                                                                                                           |

## Conventions

- Subtype colours are fixed once and shared by Figures 4 and 5, so a subtype
  keeps one colour throughout the paper.
- Panels carry no in-figure title.
- Cluster numbering from the clustering step is arbitrary. Every script that
  needs a subtype letter reads `unsupervised_subtyping/_subtype_map.R`
  rather than hard-coding the mapping.
