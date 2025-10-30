#!/bin/bash -l
#$ -l tmem=32G              # total memory available per core
#$ -l h_vmem=32G            # hard limit memory
#$ -l h_rt=1:00:00         # maximum runtime (12 hours)
#$ -cwd                     # run in current working directory
#$ -j y                     # merge stdout and stderr
#$ -V                       # export environment
#$ -N violin_h49            # job name
#$ -o violin_h49.log        # log file

echo "===== Job started: $(date) ====="
echo "Working directory: $(pwd)"

# Load R module if needed (depends on your UCL/SAN environment)
source ~/miniconda3/bin/activate phylogeny_snp

cd /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step3_wg_background_distribution/step4_violinplotsfigS1
# Run the R script
Rscript step4_dump_violinwgdistribution_Rvisualfaster.R > error.log 2>&1

echo "===== Job finished: $(date) ====="
