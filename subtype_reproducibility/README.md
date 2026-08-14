# Subtype reproducibility

Code for the Methods subsection *Subtype reproducibility*: whether the same four
subtypes come back when something that should not matter is changed. Same
conventions as the other `code_script/` modules: no absolute or user-specific
paths (everything anchored on `$PROJECT_DIR`), `data.table` for R I/O, tabular
and figure outputs only.

This module asks whether the solution is stable. What the subtypes contain is
`../subtype_characterisation/`; the fit itself is `../unsupervised_subtyping/`.

## Files

| File                           | Role                                                                                                                                                                           |
| ------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `01_cross_sex_concordance.R` | cross-sex concordance in shared clinical space: within-sex z-scoring, cosine similarity between subtype centroids, greedy matching                                             |
| `02_reseed_stability.R`      | re-fits MOFA changing only the random seed, re-clusters, and reports the Adjusted Rand Index against the main solution. Takes`<sex> <seed>`                                  |
| `03_design_sensitivity.R`    | re-fits under a different MOFA design choice (view weighting) and re-clusters at the same k. Takes`<sex> [k]`                                                                |
| `04_design_ari.R`            | puts every comparison on one ARI scale, with a permutation null and a clustering-step bootstrap ceiling, and writes the per-participant retention rates Figure 5's panel shows |

## Run

```bash
export PROJECT_DIR=/path/to/project
# none of these need _subtype_map.R: they work on cluster ids, not letters

Rscript 01_cross_sex_concordance.R

# repeat across seeds; each run appends one row to the summary
for SEED in 10001 10002 10003; do
  Rscript 02_reseed_stability.R M $SEED
  Rscript 02_reseed_stability.R F $SEED
done

Rscript 03_design_sensitivity.R M
Rscript 03_design_sensitivity.R F

Rscript 04_design_ari.R
```

`02` and `03` each re-run a whole MOFA fit, so they cost one fit per call.

## The design grid

The manuscript reports a 2×2 of two orthogonal design choices — view weighting
and feature breadth:

|                                 | default weighting             | balanced weighting          |
| ------------------------------- | ----------------------------- | --------------------------- |
| **main feature set**      | the main analysis (reference) | `03_design_sensitivity.R` |
| **broadened feature set** | separate fit                  | separate fit                |

k is fixed rather than re-selected in each cell: the comparison is between
partitions at the same granularity. The natural k of the alternative design is a
separate question and is written to the sweep table, not used for the ARI.

The broadened-feature cells relax the feature-inclusion threshold and add a
view. They are fitted separately; `04_design_ari.R` reads their cluster
assignments from `output/subtyping/clusters/` and reports each cell it finds.
