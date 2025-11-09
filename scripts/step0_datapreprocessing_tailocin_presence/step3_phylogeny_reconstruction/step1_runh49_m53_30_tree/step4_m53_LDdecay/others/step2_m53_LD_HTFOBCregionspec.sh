#!/bin/bash -l
#$ -l tmem=16G
#$ -l h_vmem=16G
#$ -l h_rt=8:00:0
#$ -cwd
#$ -V
#$ -N LD_HTF_LPS
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
OUT_A="${OUT_MAIN}/HTF_LPS_only_1kbbuffer"
mkdir -p "${OUT_A}"

CHR="utg000001l_p25.C2"

# ------------------------------------------------------------
# PLAN A : HTF–LPS linkage (real coordinates ±100 kb)
# ------------------------------------------------------------
echo "===== PLAN A : HTF–LPS LD calculation ====="

# --- HTF (2712429–2714231) ±1 kb ---
HTF_START=2711429
HTF_END=2715231

# --- LPS cluster genes range (1012835–1032652) ±1 kb ---
LPS_START=1011835
LPS_END=1033652

# Extract ±100 kb around both regions
bcftools view -r ${CHR}:${HTF_START}-${HTF_END} "${VCF}" -Oz -o "${OUT_A}/HTF_region.vcf.gz"
bcftools view -r ${CHR}:${LPS_START}-${LPS_END} "${VCF}" -Oz -o "${OUT_A}/LPS_region.vcf.gz"
bcftools index -f "${OUT_A}/HTF_region.vcf.gz"
bcftools index -f "${OUT_A}/LPS_region.vcf.gz"
# Merge and index
bcftools concat -a "${OUT_A}/HTF_region.vcf.gz" "${OUT_A}/LPS_region.vcf.gz" \
    -Oz -o "${OUT_A}/HTF_LPS_combined.vcf.gz"
bcftools index -f "${OUT_A}/HTF_LPS_combined.vcf.gz"

# Compute LD within 2 Mb window
vcftools --gzvcf "${OUT_A}/HTF_LPS_combined.vcf.gz" \
         --hap-r2 \
         --ld-window 999999 \
         --ld-window-bp 2000000 \
         --out "${OUT_A}/LD_HTF_LPS_2Mb"

echo "✅ Focused HTF–LPS LD file: ${OUT_A}/LD_HTF_LPS_22Mb.geno.ld"


