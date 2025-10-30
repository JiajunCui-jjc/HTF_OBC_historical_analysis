#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=24:0:0
#$ -S /bin/bash
#$ -N dumpmixplot97
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/
#$ -wd /SAN/ugi/plant_genom/jiajucui/

source ~/miniconda3/bin/activate phylogeny_snp

Rscript /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step4_query_and_filterby_wgthreshold/step2_mix_HTF_wg_kmerdistributions.R
