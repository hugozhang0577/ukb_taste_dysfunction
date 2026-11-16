#!/bin/bash
# Per-subtype GWAS (each subtype vs the shared controls) in SAIGE 1.5.0, with the primary GWAS covariates, GRM
# and step-2 settings; step 2 reads the minor-allele BGEN files written by saige/02_run_saige_primary.sh.
THREADS=16
GRM=sparseGRM_relatednessCutoff_0.125_5000_randomMarkersUsed.sparseGRM.mtx
COVAR=sex_plink,batch_number,smoking,drink,surg_taste_affecting_full,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,age
QCOVAR=sex_plink,smoking,drink,surg_taste_affecting_full

for S in A B C D; do
  step1_fitNULLGLMM.R --sparseGRMFile=$GRM --sparseGRMSampleIDFile=$GRM.sampleIDs.txt --useSparseGRMtoFitNULL=TRUE \
    --plinkFile=random_50000_markers --phenoFile=pheno_subtype_${S}.txt --phenoCol=pheno \
    --sampleIDColinphenoFile=IID --traitType=binary --covarColList=$COVAR --qCovarColList=$QCOVAR \
    --skipVarianceRatioEstimation=FALSE --IsOverwriteVarianceRatioFile=TRUE --nThreads=$THREADS \
    --outputPrefix=step1_subtype_${S}

  # Step 2 is run with Allele2 as the minor allele.
  for CHR in {1..22}; do
    step2_SPAtests.R --bgenFile=ukb22828_c${CHR}_b0_v3_minor.bgen --bgenFileIndex=ukb22828_c${CHR}_b0_v3_minor.bgen.bgi \
      --sampleFile=ukb22828_c${CHR}_b0_v3_minor.sample --AlleleOrder=ref-first --chrom=${CHR} \
      --GMMATmodelFile=step1_subtype_${S}.rda --varianceRatioFile=step1_subtype_${S}.varianceRatio.txt \
      --SAIGEOutputFile=subtype_${S}_chr${CHR}.txt --minMAC=20 --LOCO=FALSE \
      --is_Firth_beta=TRUE --pCutoffforFirth=0.05 --is_fastTest=FALSE --SPAcutoff=0.5
  done
done
