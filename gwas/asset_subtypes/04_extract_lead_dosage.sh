#!/bin/bash
# Dosages of the lead variants of the suggestive (P < 1e-6) ASSET loci outside APOE, from the QC-filtered imputed data.
while read -r CHR LEAD; do
  plink2 --bgen ukb22828_c${CHR}_b0_v3_qc_filtered.bgen ref-first --sample ukb22828_c${CHR}_b0_v3_qc_filtered.sample \
    --snp ${LEAD} --export A --out lead_dosage_${LEAD}
done << EOF
2 rs13028944
20 rs7269148
6 rs59668897
EOF
