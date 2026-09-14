#!/bin/bash
awk '$8 >= 0.8 {print $2}' ukb22828_c22_b0_v3.mfi.txt > c22_info08.snplist && plink2 --bgen ukb22828_c22_b0_v3.bgen ref-first --sample ukb22828_c22_b0_v3.sample --keep taste_2w_strict_white_keep_ids.txt --pheno taste_2w_strict_white_gwas_pheno.txt --pheno-name pheno --extract c22_info08.snplist --maf 0.001 --geno 0.01 --hwe 1e-6 --rm-dup exclude-all --export bgen-1.2 bits=8 --threads "$(nproc)" --memory "$(( $(awk '/MemTotal/{print int($2*0.8/1024)}' /proc/meminfo) ))" --out ukb22828_c22_b0_v3_qc_filtered && bgenix -g ukb22828_c22_b0_v3_qc_filtered.bgen -index -clobber

# All 22 in one job (only if every chromosome is attached to that job):

for CHR in {1..22}; do awk '$8 >= 0.8 {print $2}' ukb22828_c${CHR}_b0_v3.mfi.txt > c${CHR}_info08.snplist && plink2 --bgen ukb22828_c${CHR}_b0_v3.bgen ref-first --sample ukb22828_c${CHR}_b0_v3.sample --keep taste_2w_strict_white_keep_ids.txt --pheno taste_2w_strict_white_gwas_pheno.txt --pheno-name pheno --extract c${CHR}_info08.snplist --maf 0.001 --geno 0.01 --hwe 1e-6 --rm-dup exclude-all --export bgen-1.2 bits=8 --threads "$(nproc)" --memory "$(( $(awk '/MemTotal/{print int($2*0.8/1024)}' /proc/meminfo) ))" --out ukb22828_c${CHR}_b0_v3_qc_filtered && bgenix -g ukb22828_c${CHR}_b0_v3_qc_filtered.bgen -index -clobber; done

plink2 --bfile ukb22418_c22_b0_v2 --keep taste_2w_strict_white_keep_ids.txt --pheno taste_2w_strict_white_gwas_pheno.txt --pheno-name pheno --no-psam-pheno --maf 0.01 --geno 0.01 --hwe 1e-6 --rm-dup exclude-all --make-bed --threads "$(nproc)" --memory "$(( $(awk '/MemTotal/{print int($2*0.8/1024)}' /proc/meminfo) ))" --out ukb22418_c22_b0_v2_filt

# All 22 in one job:

for CHR in {1..22}; do plink2 --bfile ukb22418_c${CHR}_b0_v2 --keep taste_2w_strict_white_keep_ids.txt --pheno taste_2w_strict_white_gwas_pheno.txt --pheno-name pheno --no-psam-pheno --maf 0.01 --geno 0.01 --hwe 1e-6 --rm-dup exclude-all --make-bed --threads "$(nproc)" --memory "$(( $(awk '/MemTotal/{print int($2*0.8/1024)}' /proc/meminfo) ))" --out ukb22418_c${CHR}_b0_v2_filt; done

for c in $(seq 2 22); do echo ukb22418_c${c}_b0_v2_filt; done > merge_list.txt && plink2 --bfile ukb22418_c1_b0_v2_filt --merge-list merge_list.txt --make-bed --threads "$(nproc)" --memory "$(( $(awk '/MemTotal/{print int($2*0.8/1024)}' /proc/meminfo) ))" --out maf_flt_22chroms

plink2 --bfile maf_flt_22chroms --maf 0.05 --indep-pairwise 200 50 0.2 --threads "$(nproc)" --memory "$(( $(awk '/MemTotal/{print int($2*0.8/1024)}' /proc/meminfo) ))" --out grm_prune && plink2 --bfile maf_flt_22chroms --extract grm_prune.prune.in --keep taste_2w_strict_white_keep_ids.txt --make-bed --threads "$(nproc)" --memory "$(( $(awk '/MemTotal/{print int($2*0.8/1024)}' /proc/meminfo) ))" --out grm_for_saige_taste_2w_strict_white
