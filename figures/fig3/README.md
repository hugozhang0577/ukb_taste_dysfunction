# Figure 3 — exposure-wide and disease-wide associations

Four panels: a volcano and a signal-count bar for the questionnaire exposures,
and the same pair for the disease phenome.

The figure is one script, but it is preceded by two preparation steps, because
the three scans it draws together do not arrive in a comparable state. Those two
steps live here rather than in a module of their own: nothing else in the package
reads their output.

| File | Role |
|---|---|
| `tier_flagged_results.R` | reads the primary result of each scan for the discovery and the two held-out cohorts, harmonises differently-shaped outputs into one long table with a common schema, and flags each row on four tiers |
| `exwas_continuous_sd.R` | per-exposure standard deviation in the discovery cohort, so continuous effects can be drawn per SD rather than per raw unit |
| `exwas_dwas.R` | the four panels |

## Run

```bash
export PROJECT_DIR=/path/to/project

Rscript figures/fig3/tier_flagged_results.R
Rscript figures/fig3/exwas_continuous_sd.R
Rscript figures/fig3/exwas_dwas.R
```

The two preparation steps are cheap and only need re-running when a scan is
re-run.

## Why the preparation steps exist

The two questionnaire scans and the disease scan report different columns, apply
FDR within their own families, and keep their held-out results in separate files.
Panels A–D put all three on one pair of axes. Harmonising once means the figure
selects and colours rows rather than re-deriving significance, and the
supplementary tables built from the same table cannot disagree with the figure
about what passed.

Figure 2 needs none of this: its panels are one scan each, so it reads each
scan's own summary statistics directly.

## Inputs

Each scan writes under `output/<scan>/`; copy its primary result into
`input/assoc_results/` under the fixed name below, the same staging convention
`input/analysis_ready/` uses for the derived matrices.

| Expected file | Comes from |
|---|---|
| `pwas_primary.csv` | `../../pwas/02_run_pwas.R` |
| `mwas_primary.csv` | `../../mwas/` |
| `exwas_baseline_primary.csv` | `../../exwas_baseline/02_run_exwas.R` |
| `exwas_followup_primary.csv` | `../../exwas_followup/08_run_exwas.R` |
| `dwas_phecode_primary.csv` | `../../dwas/04_run_dwas.R` |
| `*_{group2,group3}_primary.csv` | the same scans on the held-out cohorts |

A held-out file that is absent is reported and skipped, not fatal: the
replication flags for that layer are then simply not set.

`exwas_continuous_sd.R` additionally reads the two ExWAS analysis-ready matrices
and their variable dictionaries from `input/analysis_ready/`.

Outputs land in `output/evidence_tiering/`, and the figure in
`output/figures/fig3/`.

## The four tiers, and what they are not

| Flag | Meaning |
|---|---|
| `pass_bonf` | within-domain Bonferroni q < 0.05 |
| `pass_fdr` | within-domain Benjamini-Hochberg q < 0.05 |
| `pass_effect_tier` | \|beta\| at or above a layer-specific floor (odds ratio 1.10 for the exposure and protein scans, 0.05 SD for the metabolite scan) |
| `pass_g2_replication` / `pass_g3_replication` | same direction in a held-out cohort, and nominally significant there when available |

`evidence_score` sums four of these and `headline_flag` marks a score of at least
three. Two cautions on reading them:

- **The flags are not independent tests.** Bonferroni implies FDR, so the score
  is an ordering device for reading a long table, not a probability, and a score
  of 4 is not "four separate confirmations".
- **`pass_effect_tier` is scale-dependent.** It is applied to each layer's beta
  in that layer's own units, so for a continuous exposure whether it passes
  depends on the unit the exposure was measured in. `exwas_continuous_sd.R` is
  what puts those exposures on a comparable per-SD axis, and it is applied when
  the figure is drawn — after these flags, not before them.

FDR is computed within domain and never pooled across domains or layers, matching
how each scan reports its own results.

The genome-wide scan is not tiered. Ten million variants do not belong in a
row-per-variable evidence table, and genome-wide significance is its own
threshold; Figure 2 reads the lead variants directly.

## Panel selection

Both panel pairs select on the `layer` column — `ExWAS-baseline` and
`ExWAS-followup` for A and B, `DWAS` for C and D — which is the same column the
load step filters on, so the panels cannot end up selecting on a vocabulary the
tiering step has moved on from.
