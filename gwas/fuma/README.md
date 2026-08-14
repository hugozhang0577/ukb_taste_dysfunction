# Functional mapping (FUMA)

Post-GWAS functional annotation of the primary signal, reported under
*Functional mapping*. SNP2GENE is run on the **FUMA web platform**; this folder
holds the input-preparation and result-summarisation code. No absolute paths
beyond the `PROJECT_DIR` default.

## Files

| File | Role |
|---|---|
| `01_prepare_fuma_input.R` | SAIGE sumstats → FUMA upload format (A1 = SAIGE effect allele = Allele2; N = N_case + N_ctrl); writes `gwas_for_fuma.txt.gz`. |
| `02_summarise_fuma.R` | Reads the FUMA SNP2GENE output, integrates lead SNPs / loci / mapped genes / GWAS-Catalog overlap, writes the main-text Table 1 + Supplementary Tables S1–S5. |

## FUMA SNP2GENE configuration (v1.5.2; verified 2026-05-16)

- **Reference panel**: UK Biobank release 2b, White-British 10k (`WBrits_10k`).
- **Build**: GRCh37 / hg19.
- **Gene mapping**: positional + eQTL (`posMap` + `eqtlMap`); chromatin-interaction
  mapping OFF.
- **MHC**: excluded (`exMHC = 1`).
- **Locus merge distance**: 250 kb; **MAGMA** gene/gene-set test ON.
- **Two jobs** (one input, two lead-SNP thresholds): suggestive `1e-5` and
  genome-wide `5e-8`.

> If FUMA returns `ERROR 003`, it is usually scientific-notation BP — re-export
> with integer positions (`01_prepare_fuma_input.R` already forces this via
> `scipen`).

## Run

```bash
export PROJECT_DIR=/path/to/project
Rscript 01_prepare_fuma_input.R     # -> gwas_for_fuma.txt.gz  (upload to FUMA)
# ... run SNP2GENE on https://fuma.ctglab.nl/ , download results into FUMA_1e5/ ...
Rscript 02_summarise_fuma.R         # -> FUMA_summary/Table1 + S1-S5
```
