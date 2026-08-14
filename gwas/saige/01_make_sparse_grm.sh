#!/bin/bash
# =============================================================================
# SAIGE — sparse GRM and variance-ratio markers (step-1 prerequisites)
# =============================================================================
#
# The two inputs SAIGE step 1 needs, both built from the LD-pruned hard-call
# PLINK set that ../02_variant_qc.sh produces (grm_for_saige_<cohort>):
#
#   1. a random 50,000-marker PLINK set for variance-ratio estimation, and
#   2. a sparse genetic relationship matrix at relatedness cutoff 0.125.
#
# Both are cohort-invariant given the sample set, so step 1 for every covariate
# model reuses the same pair. Run once per cohort.
#
# Run in the Swiss Army Knife terminal with the PLINK set attached as job inputs;
# the files arrive in the working directory, so the commands use bare filenames.
# The second command needs the SAIGE image, which the job runs in rather than the
# default one, so that the version is fixed by a file rather than by whatever a
# registry currently serves.
#
# The marker thinning is seeded, so the 50,000 markers are the same on every run.
#
# Inputs:  grm_for_saige_cohort3_2w_strict.{bed,bim,fam}
# Outputs: random_50000_markers.{bed,bim,fam}
#          sparseGRM_relatednessCutoff_0.125_5000_randomMarkersUsed.sparseGRM.mtx
#          ...sparseGRM.mtx.sampleIDs.txt
# =============================================================================


# -----------------------------------------------------------------------------
# [1] 50,000 random markers for the variance ratio
# -----------------------------------------------------------------------------

plink2 --bfile grm_for_saige_cohort3_2w_strict --thin-count 50000 --seed 1766121598 --make-bed --threads "$(nproc)" --memory "$(( $(awk '/MemTotal/{print int($2*0.8/1024)}' /proc/meminfo) ))" --out random_50000_markers


# -----------------------------------------------------------------------------
# [2] sparse GRM  (run with the SAIGE 1.5.0 image selected for the job)
# -----------------------------------------------------------------------------

createSparseGRM.R --plinkFile=grm_for_saige_cohort3_2w_strict --nThreads="$(nproc)" --outputPrefix=sparseGRM --numRandomMarkerforSparseKin=5000 --relatednessCutoff=0.125
