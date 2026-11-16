#!/bin/bash
# PRSice-2 base files, one per chromosome, from the discovery GWAS (all tested variants; A1 = SAIGE Allele2).
for CHR in {1..22}; do
  awk 'BEGIN{OFS="\t"} NR==1{for(i=1;i<=NF;i++)h[$i]=i; print "SNP","CHR","BP","A1","A2","BETA","P"; next}
       {print $h["MarkerID"],$h["CHR"],$h["POS"],$h["Allele2"],$h["Allele1"],$h["BETA"],$h["p.value"]}' \
    discovery_chr${CHR}.txt > base_chr${CHR}.txt
done
