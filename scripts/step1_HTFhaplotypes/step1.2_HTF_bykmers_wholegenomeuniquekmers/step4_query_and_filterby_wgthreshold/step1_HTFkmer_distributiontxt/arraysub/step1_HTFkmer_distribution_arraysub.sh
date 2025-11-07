#!/bin/bash -l
#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=2:00:0
#$ -N hm106HTFdisttxt_array
#$ -j y
#$ -V
#$ -cwd
#$ -t 1-106
#$ -e /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step4_query_and_filterby_wgthreshold/logs
#$ -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step4_query_and_filterby_wgthreshold/logs

# ---------------------------------------------------------------
# Each array task processes one isolate (.jf) file.
# The script queries all unique kmers from all HTF references.
# ---------------------------------------------------------------

source ~/miniconda3/bin/activate phylogeny_snp
set -euo pipefail

# === Paths ===
wd="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step2_isolate_kmers"
out="$wd/../step3_query_kmers"
mkdir -p "$out" "$out/kmer_distribution_mix"
kmerref="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/kmers_unique_hamming2"

# === List all .jf files and pick based on SGE_TASK_ID ===
jf_list=($wd/isolate_jf/*.jf)
jf="${jf_list[$((SGE_TASK_ID-1))]}"
iso_name=$(basename "$jf" .jf)

echo "[$(date)] Processing sample: $iso_name (Task ID: $SGE_TASK_ID)"

# --- Output file for this isolate ---
out_file_individual="$out/kmer_distribution_mix/${iso_name}_HTF_kmer_depth.tsv"
echo -e "sample\tref\tkmer\tdepth" > "$out_file_individual"

# --- Query all HTF reference k-mers ---
for ref_kmer in "$kmerref"/*_unique_filtered_iterative2.txt; do
    ref_name=$(basename "$ref_kmer" _unique_filtered_iterative2.txt)

    while read -r kmer; do
        # Query the count of each k-mer in this isolate’s Jellyfish db
        count=$(jellyfish query "$jf" "$kmer" 2>/dev/null | awk '{print $2}')
        if [[ -z "$count" ]]; then
            count=0
        fi
        echo -e "${iso_name}\t${ref_name}\t${kmer}\t${count}" >> "$out_file_individual"
    done < "$ref_kmer"
done

echo "[$(date)] done! Output saved to: $out_file_individual"

