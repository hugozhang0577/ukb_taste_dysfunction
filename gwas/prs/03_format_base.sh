#!/bin/bash
zcat gwas_sumstats_discovery70.tsv.gz | awk 'BEGIN{OFS="\t"} NR==1{for(i=1;i<=NF;i++)h[$i]=i; print "SNP","CHR","BP","A1","A2","BETA","SE","P"; next} $h["p.value"]!="NA" && $h["BETA"]!="NA" && $h["SE"]!="NA" && $h["p.value"]+0<=0.5 {print $h["MarkerID"],$h["CHR"],$h["POS"],$h["Allele2"],$h["Allele1"],$h["BETA"],$h["SE"],$h["p.value"]}' > gwas_base.txt && wc -l gwas_base.txt

# Split by chromosome, so that scoring can go one chromosome at a time

for CHR in {1..22}; do awk -v c=${CHR} 'NR==1{print; next} $2==c' gwas_base.txt > gwas_base_chr${CHR}.txt; echo "chr${CHR}: $(( $(wc -l < gwas_base_chr${CHR}.txt) - 1 )) variants"; done
