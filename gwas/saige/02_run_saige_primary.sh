#!/bin/bash
step1_fitNULLGLMM.R --sparseGRMFile=sparseGRM_relatednessCutoff_0.125_5000_randomMarkersUsed.sparseGRM.mtx --sparseGRMSampleIDFile=sparseGRM_relatednessCutoff_0.125_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt --useSparseGRMtoFitNULL=TRUE --plinkFile=random_50000_markers --phenoFile=taste_2w_strict_white_gwas_pheno.txt --phenoCol=pheno --sampleIDColinphenoFile=IID --traitType=binary --covarColList=sex_plink,batch_number,smoking,drink,surg_taste_affecting_full,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,age --qCovarColList=sex_plink,batch_number,smoking,drink,surg_taste_affecting_full --skipVarianceRatioEstimation=FALSE --IsOverwriteVarianceRatioFile=TRUE --nThreads="$(nproc)" --outputPrefix=step1_taste_change

# Step 2 — single-variant association, chromosomes 1 to 22

for CHR in {1..22}; do step2_SPAtests.R --bgenFile=ukb22828_c${CHR}_b0_v3_qc_filtered.bgen --bgenFileIndex=ukb22828_c${CHR}_b0_v3_qc_filtered.bgen.bgi --sampleFile=ukb22828_c${CHR}_b0_v3_qc_filtered.sample --GMMATmodelFile=step1_taste_change.rda --varianceRatioFile=step1_taste_change.varianceRatio.txt --SAIGEOutputFile=taste_change_chr${CHR}.txt --chrom=${CHR} --AlleleOrder=ref-first --minMAC=20 --LOCO=FALSE --is_Firth_beta=TRUE --pCutoffforFirth=0.05 --is_fastTest=FALSE --SPAcutoff=0.5; done

# One chromosome at a time, if step 2 is split across jobs (CHR=22 shown):

step2_SPAtests.R --bgenFile=ukb22828_c22_b0_v3_qc_filtered.bgen --bgenFileIndex=ukb22828_c22_b0_v3_qc_filtered.bgen.bgi --sampleFile=ukb22828_c22_b0_v3_qc_filtered.sample --GMMATmodelFile=step1_taste_change.rda --varianceRatioFile=step1_taste_change.varianceRatio.txt --SAIGEOutputFile=taste_change_chr22.txt --chrom=22 --AlleleOrder=ref-first --minMAC=20 --LOCO=FALSE --is_Firth_beta=TRUE --pCutoffforFirth=0.05 --is_fastTest=FALSE --SPAcutoff=0.5
