#!/bin/bash -l
#$ -l tmem=56G
#$ -l h_vmem=56G
#$ -l h_rt=6:00:00
#$ -cwd
#$ -V
#$ -N LDtailocin_subset10
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/

# ============================================================
# ============================================================

echo "[$(date)] Starting LD_HTF_LPS_vs_subset10 job ..."
source ~/miniconda3/bin/activate phylogeny_snp
set -euo pipefail

# --- Define paths ---
R_SCRIPT="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/scripts/step0_datapreprocessing_tailocin_presence/step3_phylogeny_reconstruction/step1_runh49_m53_30_tree/step4_m53_LDdecay/Rvisual/Rvisual_tailocin_withwithoutHTF.R"

WORKDIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step3_phylogeny_tree/results_m53_LDdecay/LDcalculation"

cd "$WORKDIR"

echo "[$(date)] Running Rscript ..."
Rscript "$R_SCRIPT"

echo "[$(date)] Job finished successfully."
