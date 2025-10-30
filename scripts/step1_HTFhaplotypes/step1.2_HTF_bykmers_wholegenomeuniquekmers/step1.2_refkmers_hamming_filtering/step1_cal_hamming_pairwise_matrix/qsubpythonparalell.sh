#!/bin/bash -l
#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=5:00:00
#$ -cwd
#$ -pe smp 2
#$ -N h1amming_calc
#$ -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step1_cal_hamming_pairwise_matrix/hamming_matrix.log
#$ -e /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step1_cal_hamming_pairwise_matrix/hamming_matrix.err
source ~/miniconda3/bin/activate phylogeny_snp
wd=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step1.2_HTF_bykmers_wholegenomeuniquekmers/step1.2_refkmers_hamming_filtering/step1_cal_hamming_pairwise_matrix
kmer_unique=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/intermediate_kmersfiltering/kmers_unique

mkdir -p $kmer_unique/../step1_cal_hamming_pairwise_matrix/
#python $inputkmerspath $output $threshold
python $wd/calculate_hamming_matrix.py $kmer_unique $kmer_unique/../step1_cal_hamming_pairwise_matrix/pairwise_hamming_distance_longformatbypy.tsv 2

mkdir -p /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/kmers_unique_hamming2/
mv /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/intermediate_kmersfiltering/ /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/kmers_unique_hamming2/
