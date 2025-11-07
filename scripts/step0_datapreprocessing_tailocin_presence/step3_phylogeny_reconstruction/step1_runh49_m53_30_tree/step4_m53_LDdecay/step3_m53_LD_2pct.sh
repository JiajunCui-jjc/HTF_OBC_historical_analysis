#!/bin/bash -l
#$ -l tmem=16G
#$ -l h_vmem=16G
#$ -l h_rt=8:00:0
#$ -cwd
#$ -V
#$ -N LD_subset2
#$ -hold_jid m53_merge_bivcf
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/

# ============================================================
# Combined LD pipeline for m53 modern isolates
#   PLAN A : Focused LD between HTF ±100 kb and LPS ±100 kb
#   PLAN B : Genome-wide LD decay using 20 % random SNP subset
# ============================================================

source ~/miniconda3/bin/activate phylogeny_snp

# ------------------------------------------------------------
# Directories and file paths
# ------------------------------------------------------------
TOP_DIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step3_phylogeny_tree"
RES_DIR="${TOP_DIR}/results_m53_LDdecay"
VCF="${RES_DIR}/biallelic_fullinfo_m53samples_GTfixed.vcf.gz"

OUT_MAIN="${RES_DIR}/LDcalculation"
OUT_B="${OUT_MAIN}/subset2_LDdecay2mb_bg"
mkdir -p  "${OUT_B}"

CHR="utg000001l_p25.C2"

# ------------------------------------------------------------
# PLAN B : Genome-wide LD decay using 2 % subset (~4 k SNP)
# ------------------------------------------------------------
echo "===== PLAN B : 20% subset LD-decay calculation ====="

TOTAL_SNPS=$(bcftools index -n "${VCF}")
SUBSET_COUNT=$(( TOTAL_SNPS / 50 ))  # 2%

echo "→ Selecting ${SUBSET_COUNT} random SNPs (~5%)..."
bcftools query -f '%CHROM\t%POS\n' "${VCF}" | shuf -n ${SUBSET_COUNT} > "${OUT_B}/subset_positions.txt"

# Subset and index
bcftools view -T "${OUT_B}/subset_positions.txt" "${VCF}" -Oz -o "${OUT_B}/subset_2pct.vcf.gz"
bcftools index -f "${OUT_B}/subset_2pct.vcf.gz"
# ------------------------------------------------------------
# vcftools --geno-r2    → calculates r² between all SNP pairs
# --ld-window-bp 1e6    → only SNPs ≤1 Mb apart
# --ld-window 999999    → remove SNP-number restriction
# Output .geno.ld has: CHR POS1 POS2 N_INDV R2

# Compute LD within 2 Mb window
vcftools --gzvcf "${OUT_B}/subset_2pct.vcf.gz" \
         --hap-r2 \
         --ld-window 999999 \
         --ld-window-bp 2000000 \
         --out "${OUT_B}/LD_subset_2pct_2mb"


