#!/bin/bash
plink2 --bfile grm_for_saige_taste_2w_strict_white --thin-count 50000 --seed 1766121598 --make-bed --threads "$(nproc)" --memory "$(( $(awk '/MemTotal/{print int($2*0.8/1024)}' /proc/meminfo) ))" --out random_50000_markers

# [2] sparse GRM  (run with the SAIGE 1.5.0 image selected for the job)

createSparseGRM.R --plinkFile=grm_for_saige_taste_2w_strict_white --nThreads="$(nproc)" --outputPrefix=sparseGRM --numRandomMarkerforSparseKin=5000 --relatednessCutoff=0.125
