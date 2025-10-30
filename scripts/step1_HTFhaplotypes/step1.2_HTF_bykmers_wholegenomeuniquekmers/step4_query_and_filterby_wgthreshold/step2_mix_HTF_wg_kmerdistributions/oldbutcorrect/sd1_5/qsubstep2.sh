#!/bin/bash -l
#$ -l tmem=8G
#$ -l h_vmem=8G
#$ -l h_rt=1:00:0
#$ -N plot15sd_HTF_WG
#$ -j y
#$ -V
#$ -cwd
#$ -t 1-106
#$ -e /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step4_query_and_filterby_wgthreshold/step2_mix_HTF_wg_kmerdistributions/logs
#$ -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step4_query_and_filterby_wgthreshold/step2_mix_HTF_wg_kmerdistributions/logs

# ---------------------------------------------------------------
# Array job: each task plots one sample (modern or historical)
# ---------------------------------------------------------------

source ~/miniconda3/bin/activate phylogeny_snp
set -euo pipefail

# Get sample list (unique sample names)
sample_list=$(awk 'NR>1 {print $1}' /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers/kmer_distribution_mix/HTF_kmer_depth_distributions.tsv | sort -u)
samples=($sample_list)
sample_name=${samples[$((SGE_TASK_ID-1))]}

echo "[$(date)] Processing sample: $sample_name"

Rscript /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step4_query_and_filterby_wgthreshold/step2_mix_HTF_wg_kmerdistributions/sd1_5/step2_mix_HTF_wg_kmerdistributions_arraysub.R "$sample_name"

echo "[$(date)] Done: $sample_name"

