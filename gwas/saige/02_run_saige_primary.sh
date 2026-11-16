#!/bin/bash
# Primary GWAS in SAIGE 1.5.0: step 1 null logistic mixed model on the sparse GRM; step 2 single-variant
# tests per chromosome (MAC >= 20, SPA, Firth correction for P < 0.05).
THREADS=16
GRM=sparseGRM_relatednessCutoff_0.125_5000_randomMarkersUsed.sparseGRM.mtx
COVAR=sex_plink,batch_number,smoking,drink,surg_taste_affecting_full,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,age
QCOVAR=sex_plink,smoking,drink,surg_taste_affecting_full

step1_fitNULLGLMM.R --sparseGRMFile=$GRM --sparseGRMSampleIDFile=$GRM.sampleIDs.txt --useSparseGRMtoFitNULL=TRUE \
  --plinkFile=random_50000_markers --phenoFile=taste_2w_strict_white_gwas_pheno.txt --phenoCol=pheno \
  --sampleIDColinphenoFile=IID --traitType=binary --covarColList=$COVAR --qCovarColList=$QCOVAR \
  --skipVarianceRatioEstimation=FALSE --IsOverwriteVarianceRatioFile=TRUE --nThreads=$THREADS \
  --outputPrefix=step1_taste_change

# Step 2 is run with Allele2 as the minor allele.
for CHR in {1..22}; do
  plink2 --bgen ukb22828_c${CHR}_b0_v3_qc_filtered.bgen ref-first --sample ukb22828_c${CHR}_b0_v3_qc_filtered.sample \
    --maj-ref force --export bgen-1.2 ref-first bits=8 --out ukb22828_c${CHR}_b0_v3_minor
  bgenix -g ukb22828_c${CHR}_b0_v3_minor.bgen -index -clobber
done

for CHR in {1..22}; do
  step2_SPAtests.R --bgenFile=ukb22828_c${CHR}_b0_v3_minor.bgen --bgenFileIndex=ukb22828_c${CHR}_b0_v3_minor.bgen.bgi \
    --sampleFile=ukb22828_c${CHR}_b0_v3_minor.sample --AlleleOrder=ref-first --chrom=${CHR} \
    --GMMATmodelFile=step1_taste_change.rda --varianceRatioFile=step1_taste_change.varianceRatio.txt \
    --SAIGEOutputFile=taste_change_chr${CHR}.txt --minMAC=20 --LOCO=FALSE \
    --is_Firth_beta=TRUE --pCutoffforFirth=0.05 --is_fastTest=FALSE --SPAcutoff=0.5
done
