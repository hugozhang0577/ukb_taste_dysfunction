# Genome-wide association analysis

Code for the taste-dysfunction GWAS, its phenotype-definition sensitivity runs,
functional mapping, the APOE diplotype analysis and the polygenic risk score.

Phenotype, covariate and sample-QC derivation live in `../preprocessing/`
(`02_define_phenotypes.R`, `04_gwas_sample_qc.R`); this module covers the
genetic-analysis engines that consume those exports. Containers are
pinned (SAIGE 1.5.0), and paths resolve under `$PROJECT_DIR` and `$COHORT_DIR`.

## Layout

```text
gwas/
  01_sample_qc.R      individual-level QC and ancestry split
  02_variant_qc.sh    genotype prep: analysis-ready imputed bgen + LD-pruned GRM set
  saige/              the primary single-variant GWAS, summary-statistic merge,
                      and the Manhattan / QQ / lambda diagnostics

  fuma/               FUMA functional mapping
  apoe_diplotype/     APOE e2/e3/e4 diplotype association
  prs/                polygenic risk score (clumping and thresholding, PRSice-2)
```

## The primary covariate model

`sex_plink`, PC1-PC10, `age` (at the taste assessment), `batch_number`,
`smoking`, `drink` and `surg_taste_affecting_full`, with relatedness carried by
the sparse GRM. `saige/README.md` gives the full specification and which terms
are declared categorical.

Every one of these columns is written by `../preprocessing/04_gwas_sample_qc.R`
into the phenotype/covariate file the SAIGE scripts read, so the genome-wide
analysis has no covariate input this package does not produce. `age`,
`batch_number`, `smoking`, `drink`, `BMI`, `surg_taste_affecting_full` and
`smell_any` are the same derived columns the panel scans use, taken from the
covariate table rather than re-derived, so a covariate means the same thing in
every analysis.

`batch_number` (genotype measurement batch, field 22000) enters as a **factor**.
It is a batch label, not a quantity: a numeric batch code would impose an
ordering the assay does not have.

`sex_plink` is the exception and is named apart on purpose: it is genetic sex
(field 22001) in PLINK coding, 1 male / 2 female, whereas the `sex` the panel
scans use is self-reported (field 31), coded 1 male / 0 female. Different source
field, different coding, different name.

Genomic-control lambda is reported for the genome-wide analysis. It is not
reported for the panel scans (proteomic, metabolomic, questionnaire, disease):
those test outcome-anchored, health-relevant panels in which most variables are
genuinely associated, so a high lambda there reflects real signal rather than
inflation and would be misread as the opposite.

## Phenotype-definition sensitivity cohorts

Two further GWAS re-run the primary model unchanged on cohorts defined by a
different taste phenotype, testing whether the chr19/APOE signal depends on where
the case boundary is drawn rather than on the biology.

| Cohort | Outcome | Definition | Relative to primary |
|---|---|---|---|
| `cohort_sens_4w` | stricter | `taste_4w_strict`: taste change lasting at least 4 weeks **or** affecting daily life | narrower |
| `cohort_primary` | primary | `taste_2w_strict`: at least 2 weeks **or** impact | — |
| `cohort_sens_any` | broader | `taste_basic_strict`: any reported taste change | broader |

All three use smell-free controls: a participant reporting no taste change and no
smell change. The two sensitivity cohorts re-run the primary covariate model
unchanged; only the phenotype, and therefore the control set, differs.

## Run order

```bash
export PROJECT_DIR=/path/to/project
export COHORT_DIR=$PROJECT_DIR/gwas/cohort_primary

Rscript 01_sample_qc.R                  # ../preprocessing/02 and /04 make the exports
#       02_variant_qc.sh                # SAK command lines: imputed, array,
#                                       # then merge, then GRM prune
#       saige/01_make_sparse_grm.sh     # SAK command lines
#       saige/02_run_saige_primary.sh   # SAK command lines: step 1, then step 2
#       saige/03_merge_sumstats.sh      # SAK: 22 chromosomes -> one file
Rscript saige/04_manhattan_qq_lambda.R  # Manhattan, QQ, lambda
Rscript fuma/01_prepare_fuma_input.R    # upload to FUMA, then 02_summarise_fuma.R
#       apoe_diplotype/01_extract_apoe_genotype.sh    # SAK; then 02 and 03
Rscript prs/01_split_cohort.R           # then the discovery GWAS, then prs/02-06
                                        # (JupyterLab workspace, see prs/README.md)
```

## What this package ships, and what it does not

`saige/` holds the **primary** scan only: one step-1 command line and one step-2
command line. The Supplementary Methods report several variants of that same
fit, and each is the same two commands with one thing changed — no separate
code, and none is included here:

| Reported variant | What changes |
|---|---|
| behaviour-free baseline | drop `surg_taste_affecting_full` from `--covarColList` |
| metabolic sensitivity | add `BMI` to `--covarColList` |
| smell-adjusted sensitivity | add `smell_any` to `--covarColList` and `--qCovarColList` |
| stricter phenotype (`taste_4w_strict`) | `--phenoFile=cohort6_4w_strict_gwas_pheno.txt` |
| broader phenotype (`taste_basic_strict`) | `--phenoFile=cohort1_basic_strict_gwas_pheno.txt` |

Each needs its own `--outputPrefix`, and the phenotype variants need their own
genotype QC because the control set differs. `saige/04_manhattan_qq_lambda.R`
takes `SUMSTATS` and `RUN_TAG` from the environment, so the diagnostics are the
same script pointed at whichever merged file a run produced.

The smell-adjusted variant is a lower bound on the taste-specific signal rather
than an alternative primary estimate: adjusting for concurrent smell change
removes the chemosensory component the two outcomes share.

## Variant-QC thresholds

Applied in `02_variant_qc.sh` and reported in the Supplementary Methods:

- imputed bgen: INFO at least 0.8, MAF at least 0.001, call rate at least 0.99,
  Hardy-Weinberg P above 1e-6 in controls, duplicates removed;
- hard-call array: MAF at least 0.01, missingness at most 0.01, Hardy-Weinberg
  1e-6;
- GRM LD-pruning: MAF at least 0.05, `--indep-pairwise 200 50 0.2`.

## RAP persistence

Data-export points are marked `# Export and dx upload to RAP`. A bare `fwrite` on
a DNAnexus worker writes to ephemeral local disk; the upload is what makes the
result persist. Pure-computation steps need no such step.

## Genotype steps are Swiss Army Knife command lines

Every step that touches genotype data — `02_variant_qc.sh`,
`saige/01_make_sparse_grm.sh`, `apoe_diplotype/01_extract_apoe_genotype.sh` and
`../subtype_characterisation/06_hla_variant_extract.sh` — is a list of command
lines to paste into the Swiss Army Knife terminal, with that step's inputs
attached to the job. They are not local pipelines and do not orchestrate
anything.

Two things follow, and both are visible in the commands:

- **Bare filenames, not `$PROJECT_DIR` paths.** The job's inputs arrive in the
  working directory and whatever the command leaves there is saved. This is the
  only part of the package where paths are not written relative to
  `$PROJECT_DIR`, because inside a SAK job there is no project mount to be
  relative to.
- **`--threads "$(nproc)"` and the `--memory` expression read the instance the
  job is on**, so the same line works on any instance type without editing.

The per-chromosome steps are given twice: once for a single chromosome, to run
as 22 independent jobs, and once as a loop for a single job with all 22
chromosomes attached. The merge and GRM steps consume the per-chromosome
outputs, so those have to finish first.
