#!/bin/bash -l
#$ -l tmem=20G
#$ -l h_vmem=20G
#$ -l h_rt=10:00:0
#$ -N generateHTFdist
#$ -V
#$ -e /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step4_query_and_filterby_wgthreshold/logs
#$ -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step4_query_and_filterby_wgthreshold/logs

source ~/miniconda3/bin/activate phylogeny_snp
set -euo pipefail

# --- Configuration ---fff
wd=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step2_isolate_kmers/
out=$wd/../step3_query_kmers
mkdir -p $out $out/kmer_distribution_mix/
out_file=$out/kmer_distribution_mix/HTF_kmer_depth_distributions.tsv
kmerref=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/kmers_unique_hamming2
echo -e "sample\tref\tkmer\tdepth" > "$out_file"

for jf in "$wd"/isolate_jf/*.jf; do
    iso_name=$(basename "$jf" .jf)
    echo "Processing $iso_name..."
	for ref_kmer in $kmerref/*_unique_filtered_iterative2.txt; do
        ref_name=$(basename "$ref_kmer" _unique_filtered_iterative2.txt)

        while read -r kmer; do
            count=$(jellyfish query "$jf" "$kmer" 2>/dev/null | awk '{print $2}')
            if [[ -z "$count" ]]; then
                count=0
            fi
            echo -e "${iso_name}\t${ref_name}\t${kmer}\t${count}" >> "$out_file"
        done < "$ref_kmer"

    done
done

echo "✅ Done! Output saved to: $out_file"
