#!/bin/bash -l

#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=10:00:0
#$ -wd /SAN/ugi/plant_genom/jiajucui/
#$ -V
#$ -N pankmer
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/

source ~/miniconda3/bin/activate phylogeny_snp
wd=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2026_PNAS_revision/results/step4_pankmer
mkdir -p $wd
## 1.1. Create indices
pankmer index -t 10 --rounds 10 -g /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2026_PNAS_revision/data/pankmer/pankmer_input_fastas/only_HTF -o $wd/only_HTF.tar

## 1.2. Create adjacency matrices
pankmer adj-matrix -i $wd/only_HTF.tar -o $wd/only_HTF.adjmatrix.tsv

## 1.3. Create Heatmap with Jacard distances
pankmer clustermap -i $wd/only_HTF.adjmatrix.tsv -o $wd/only_HTF.adjmatrix.jaccard.svg --metric jaccard --width 20 --height 20
