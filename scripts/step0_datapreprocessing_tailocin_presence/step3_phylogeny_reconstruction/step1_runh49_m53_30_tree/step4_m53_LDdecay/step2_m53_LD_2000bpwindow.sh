#!/bin/bash -l
#$ -l tmem=16G
#$ -l h_vmem=16G
#$ -l h_rt=8:00:0
#$ -cwd
#$ -V
#$ -N LD_2kbwidnow
#$ -hold_jid m53_merge_bivcf
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/

# ============================================================

source ~/miniconda3/bin/activate phylogeny_snp

# ------------------------------------------------------------
# Directories and file paths
# ------------------------------------------------------------
TOP_DIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step3_phylogeny_tree"
RES_DIR="${TOP_DIR}/results_m53_LDdecay"
VCF="${RES_DIR}/biallelic_fullinfo_m53samples_GTfixed.vcf.gz"

OUT_MAIN="${RES_DIR}/LDcalculation"
OUT_B="${OUT_MAIN}/fullLDdecay_2kbwindow"
mkdir -p  "${OUT_B}"

CHR="utg000001l_p25.C2"

# ------------------------------------------------------------
# PLAN B : Genome-wide LD decay using 20 % subset (~40 k SNP)
# ------------------------------------------------------------
echo "===== PLAN B : 20% subset LD-decay calculation ====="

TOTAL_SNPS=$(bcftools index -n "${VCF}")


# ------------------------------------------------------------
# vcftools --geno-r2    → calculates r² between all SNP pairs
# --ld-window-bp 1e6    → only SNPs ≤1 Mb apart
# --ld-window 999999    → remove SNP-number restriction
# Output .geno.ld has: CHR POS1 POS2 N_INDV R2

# Compute LD within 2 Mb window
vcftools --gzvcf "${VCF}" \
         --hap-r2 \
         --ld-window 999999 \
         --ld-window-bp 2000 \
         --out "${OUT_B}/LD_2kb"


