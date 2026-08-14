# Polygenic risk score (clumping and thresholding, PRSice-2)

A C+T score built from a discovery GWAS in 70% of the cohort and evaluated in the
held-out 30%. The reported result is that the score does not discriminate, which
is why it is a short paragraph rather than a figure.

This module runs in a **JupyterLab workspace on RAP**, not as Swiss Army Knife
jobs. PRSice-2 is a precompiled binary with no platform package, and the scoring
loop wants a shell that persists across chromosomes, which is what a notebook
session gives. The `.sh` files here are command lines to paste into the
JupyterLab terminal; `dx download` brings inputs into the workspace and
`dx upload` puts results back, because workspace storage does not survive the
session.

Only the C+T pipeline is here. The separate PGS-Catalog scoring workflow is not
part of this package.

| File | Role |
|---|---|
| `01_split_cohort.R` | one stratified 70/30 split; writes the discovery phenotype file and the held-out subset |
| `02_setup_prsice.sh` | unzip the PRSice-2 release into the workspace and make it executable |
| `03_format_base.sh` | discovery SAIGE summary statistics into PRSice base format, split by chromosome |
| `04_score_by_chr.sh` | PRSice-2 scores at every threshold, one chromosome at a time |
| `05_merge_scores.R` | sum the per-chromosome scores into the genome-wide score |
| `06_regression.R` | score against outcome, in the held-out subset only |

## The discovery GWAS is not in this directory

Between `01` and `03` there is a GWAS, and it is not a separate piece of code: it
is `../saige/02_run_saige_primary.sh` with

```text
--phenoFile=discovery_70_gwas_pheno.txt
--outputPrefix=step1_taste_change_discovery70
```

followed by `../saige/03_merge_sumstats.sh`. Only the sample set changes; the
covariate model, the QC and the engine settings are the primary ones.

## The one thing this design exists to prevent

Weights estimated in the same people the score is evaluated in make the score
partly a record of those people's own outcomes. The result is an R-squared near 1
and an AUC near 1 — which reads as a spectacular finding and is in fact the
signature of the error.

Three things keep that from happening quietly:

- the split is drawn once, by `01`, from a fixed seed, and written to disk;
- the two subsets are checked for overlap before they are written;
- `06` refuses to run when the phenotype file it is given is much larger than the
  held-out subset should be, rather than trusting the caller to pass the right
  file.

Scoring, unlike evaluation, is run on everybody. That is safe — the weights still
come only from the discovery subset — and it avoids scoring the cohort twice.

## Parameters

| Setting | Value |
|---|---|
| clumping | `--clump-kb 250 --clump-r2 0.1 --clump-p 1` |
| thresholds | 1e-5, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1 |
| ambiguous variants | kept (`--keep-ambig`) |
| regression at scoring time | off (`--no-regress`) |
| covariates in the evaluation | none |

Clumping is what stops one strong locus contributing once per correlated variant.
The base file is pre-filtered at P ≤ 0.5, which is a speed measure only: those
variants cannot enter any reported threshold.

`--keep-ambig` keeps strand-ambiguous A/T and C/G variants. Base and target are
the same UK Biobank imputation here, so strand cannot disagree between them and
dropping those variants would only lose information; the flag would be the wrong
choice for a base file from an external study.

The score is evaluated unadjusted. Adding ancestry components and age would if
anything flatter a score whose honest performance is the finding.

## Column mapping

SAIGE reports the effect on `Allele2`, so the base file sets `A1 = Allele2` and
`A2 = Allele1`. Reversing this flips every weight and yields a score of identical
magnitude pointing the wrong way, with nothing in the output to reveal it.

## Practical notes

- `unzip` does not reliably preserve the executable bit, so `chmod +x` is part of
  `02` rather than an afterthought.
- Whole-chromosome bgen files are read from workspace disk, not across the
  read-only project mount, which is why `04` downloads each chromosome and
  deletes it after scoring.
- `04` over all 22 chromosomes runs for hours, so it is started under `nohup` and
  survives the notebook session closing.
