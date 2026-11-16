#!/bin/bash
# Marker set for the SAIGE variance ratio (50,000 random pruned array variants) and the sparse
# genetic relationship matrix (5,000 random markers, relatedness cutoff 0.125), SAIGE 1.5.0.
THREADS=16
plink2 --bfile grm_for_saige_taste_2w_strict_white --thin-count 50000 --seed 1766121598 --make-bed \
  --out random_50000_markers
createSparseGRM.R --plinkFile=grm_for_saige_taste_2w_strict_white --nThreads=$THREADS --outputPrefix=sparseGRM \
  --numRandomMarkerforSparseKin=5000 --relatednessCutoff=0.125
