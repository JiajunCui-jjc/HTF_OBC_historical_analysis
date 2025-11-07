#!/bin/bash -l
#$ -l tmem=16G
#$ -l h_vmem=16G
#$ -l h_rt=8:00:0
#$ -cwd
#$ -V
#$ -N fullLD_tailocin
#$ -hold_jid m53_merge_bivcf
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/

# ============================================================
# Combined LD pipeline for m53 modern isolates
#   PLAN A : Focused LD between tailocin ±100 kb and LPS ±100 kb
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
OUT_A="${OUT_MAIN}/tailocin_only_region"
mkdir -p "${OUT_A}"

CHR="utg000001l_p25.C2"

# ------------------------------------------------------------
# PLAN A : tailocin–LPS linkage (real coordinates ±100 kb)
# ------------------------------------------------------------
echo "===== PLAN A : tailocin–LPS LD calculation ====="

# --- tailocin (2696041-2716969) 20928bp ---
tailocin_START=2696041
tailocin_END=2716969


# Extract ±100 kb around both regions
bcftools view -r ${CHR}:${tailocin_START}-${tailocin_END} "${VCF}" -Oz -o "${OUT_A}/tailocin_region.vcf.gz"
bcftools index -f "${OUT_A}/tailocin_region.vcf.gz"

# Compute LD within 2 kb window
vcftools --gzvcf "${OUT_A}/tailocin_region.vcf.gz" \
         --hap-r2 \
         --ld-window 999999 \
         --ld-window-bp 20928 \
         --out "${OUT_A}/LD_tailocin_20kb"



