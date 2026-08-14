# SAIGE single-variant GWAS

The primary genome-wide scan: a logistic mixed model with a saddlepoint
approximation, which is what makes the test usable at this case-control ratio
(5,753 cases against 152,617 controls, roughly 1:27) while carrying relatedness
through a sparse genetic relationship matrix rather than by dropping relatives.

SAIGE is installed on the RAP platform, so `01` and `02` are command lines to
paste into the Swiss Army Knife terminal with SAIGE selected as the job's
application. The job's inputs arrive in the working directory, so the commands
use bare filenames and no paths — this and `../02_variant_qc.sh` are the only
parts of the package not written relative to `$PROJECT_DIR`, because inside a SAK
job there is no project mount to be relative to. `04` is ordinary R and does use
`$PROJECT_DIR`.

| File | Role |
|---|---|
| `01_make_sparse_grm.sh` | the two step-1 prerequisites: 50,000 seeded random markers for the variance ratio, and the sparse GRM at relatedness cutoff 0.125 |
| `02_run_saige_primary.sh` | step 1 (null model) and step 2 (single-variant test, chromosomes 1-22) |
| `03_merge_sumstats.sh` | the 22 per-chromosome files into one `gwas_sumstats.tsv.gz` |
| `04_manhattan_qq_lambda.R` | Manhattan, QQ, genomic-control lambda, top hits |

## The primary model

| Term | Source | Declared |
|---|---|---|
| `sex_plink` | genetic sex, field 22001, PLINK coding 1 male / 2 female | categorical |
| `PC1`-`PC10` | genetic principal components, field 22009 | continuous |
| `age` | age at the taste assessment | continuous |
| `batch_number` | genotyping batch, field 22000 | categorical |
| `smoking`, `drink` | non-current smoking, alcohol risk | categorical |
| `surg_taste_affecting_full` | cumulative taste-affecting surgery | categorical |

Two of these are easy to get wrong and are worth stating plainly:

- **`batch_number` is declared categorical.** It is a batch label, not a
  quantity; entering it as a number would impose an ordering the assay does not
  have. Nothing about how the column is stored enforces this — the
  `--qCovarColList` declaration is what keeps it out of the model as a number.
- **`sex_plink` is not `sex`.** It is genetic sex (field 22001, coded 1 male /
  2 female), whereas the `sex` the panel scans use is self-reported (field 31,
  coded 1 male / 0 female). Different field, different coding, so it carries a
  different name on purpose.

Every one of these columns is written by `../../preprocessing/04_gwas_sample_qc.R`
into the phenotype file, so the genome-wide analysis has no covariate input this
package does not produce, and a covariate means the same thing here as in the
panel scans.

Age is the one covariate that differs from the panel scans by design: the GWAS
uses age at the taste assessment, the panel scans use age at recruitment, because
their exposures were measured at recruitment. Genotype has no measurement date.

## Engine settings, and why

| Setting | Value | Reason |
|---|---|---|
| `--useSparseGRMtoFitNULL` | TRUE | relatedness is carried by the sparse GRM rather than by excluding relatives |
| `--minMAC` | 20 | the association-level variant filter; below this the test is unstable whatever the model |
| `--is_Firth_beta` / `--pCutoffforFirth` | TRUE / 0.05 | the saddlepoint approximation gives a calibrated P value but a biased effect size under strong imbalance, so variants reaching P < 0.05 get a Firth-corrected estimate |
| `--LOCO` | FALSE | leave-one-chromosome-out is off: with a sparse GRM built from LD-pruned markers the proximal contamination it guards against is small, and it would cost one null model per chromosome |
| `--SPAcutoff` | 0.5 | where the saddlepoint approximation takes over from the normal approximation |
| `--AlleleOrder` | ref-first | matches how `../02_variant_qc.sh` exported the bgen; a mismatch would flip every effect direction silently |

Genome-wide significance is the conventional 5×10⁻⁸ and suggestive 1×10⁻⁵.

## Outputs

```text
step1_taste_change.rda                  null model
step1_taste_change.varianceRatio.txt    variance ratio -> step 2
taste_change_chr<N>.txt                 one file per chromosome   (02)
gwas_sumstats.tsv.gz                    merged, genome-wide       (03)
output/gwas_plots/                      Manhattan, QQ, lambda     (04)
```

## Run

```text
01_make_sparse_grm.sh     once per cohort
02_run_saige_primary.sh   step 1, then step 2 for chromosomes 1-22
03_merge_sumstats.sh      after all 22 chromosomes exist
```

Then `Rscript 04_manhattan_qq_lambda.R`, which defaults to the primary run and
takes `SUMSTATS` and `RUN_TAG` from the environment for any other.

Genomic-control lambda is reported for this analysis. It is deliberately not
reported for the panel scans: those test outcome-anchored panels in which most
variables are genuinely associated, so a high lambda there reflects real signal
and would be misread as inflation.

## Only the primary scan is here

The Supplementary Methods report several variants of this fit — three covariate
specifications and two alternative phenotype definitions. Each is these same two
command lines with one argument changed, so none is shipped as separate code;
`../README.md` lists exactly which argument each one alters.
