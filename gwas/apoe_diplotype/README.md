# APOE diplotype — taste-dysfunction association (Fig 3d)

The chr19 genome-wide significant signal tags the APOE locus. This module
resolves the **APOE e2/e3/e4 diplotype** from the two defining variants and
tests it against the primary taste outcome — the data behind main-text Fig 3d.

## Files

| File | Role |
|---|---|
| `01_extract_apoe_genotype.sh` | PLINK2 `--export A` of rs429358 + rs7412 (chr19 imputed bgen) → additive `.raw` dosages. |
| `02_apoe_association.R` | Infer the 6-level diplotype; test e4-carrier (binary), per-e4-allele dose (trend), and 6-level diplotype (e3/e3 reference) against the primary outcome with the primary covariate set; export the forest + frequency tables. |
| `03_forest_plot.R` | Publication forest plot (left table + right forest) from the association output. |

## Variants & encoding

- rs429358 (19:45411941) + rs7412 (19:45412079), GRCh37/hg19.
- **ALT-allele assumption** (as used for the reported
  OR): in `apoe_genotypes.raw`, rs429358 is coded to the **T** allele and rs7412
  to the **C** allele → `e4_count = 2 - rs429358_dose`, `e2_count = 2 -
  rs7412_dose`. Verify against the `.raw` header (column suffix = ALT) before
  reuse.
- Covariates = the primary set **without the genotyping batch**: `age, sex_plink,
  PC1–10, smoking, drink, surg_taste_affecting_full`. The batch term is a
  single-variant artefact control that a two-variant diplotype model taken from
  the same genotype file does not need.
- Reported direction: clean additive e4 dose effect (per-allele OR ≈ 1·15);
  e2 shows no signal.

## Run

```bash
export PROJECT_DIR=/path/to/project
export COHORT_DIR=$PROJECT_DIR/gwas/cohort_primary
bash    01_extract_apoe_genotype.sh   # -> apoe_genotypes.raw
Rscript 02_apoe_association.R          # -> Results/*.csv  (dx upload on RAP)
Rscript 03_forest_plot.R              # -> Figures/APOE_forest_plot.{png,pdf}
```

> Labels are ASCII (`e3/e4`, not the Greek glyph) to stay cairo_pdf-safe and to
> match the locked Fig 3 convention.
