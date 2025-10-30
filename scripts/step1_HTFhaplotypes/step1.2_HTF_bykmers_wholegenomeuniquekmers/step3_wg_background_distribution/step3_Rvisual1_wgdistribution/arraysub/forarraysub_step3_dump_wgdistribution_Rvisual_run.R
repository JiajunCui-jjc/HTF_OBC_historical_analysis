#$ -l tmem=4G
#$ -l h_vmem=4G
#$ -l h_rt=2:0:0
#$ -S /bin/bash
#$ -N dumpplot106
#$ -t 1-106
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/
#$ -wd /SAN/ugi/plant_genom/jiajucui/

echo "Task id is $SGE_TASK_ID"
source ~/miniconda3/bin/activate phylogeny_snp

Rscript /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step3_wg_background_distribution/step3_Rvisual1_wgdistribution/arraysub/forarraysub_step3_dump_wgdistribution_Rvisual.R $SGE_TASK_ID
