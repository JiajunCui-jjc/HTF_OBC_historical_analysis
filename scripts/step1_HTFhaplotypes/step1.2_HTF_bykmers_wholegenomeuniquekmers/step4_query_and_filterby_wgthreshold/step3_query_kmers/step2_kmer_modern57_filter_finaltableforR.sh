#!/bin/bash -l
#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=5:00:0
#$ -wd /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers/
#$ -N m57query
#$ -V
#$ -e /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers/logs
#$ -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers/logs
source ~/miniconda3/bin/activate phylogeny_snp
set -euo pipefail


wd=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step2_isolate_kmers/
out=$wd/../step3_query_kmers
summary_out=${out}/summarymodern57_tailocin_kmeravgwgdepth_filter.tsv
kmerref=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/kmers_unique_hamming2
cd $out

# --- Configuration ---fff
KMER=31
JFSIZE=300000000


dump_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step2_isolate_kmers/isolate_jf/dump"
thresholds_file="${dump_dir}/threshold_m57.tsv"
#use 65% threshold for historical




summary_out=${out}/summarymodern57_tailocin_kmeravgwgdepth_filter.tsv

echo "[5] Counting isolate kmers and comparing to reference alleles"
echo -e "Isolate\tReference\tTotal_Allele_Kmers\tMatched_Kmers\tProp" > $summary_out

  
# === Loop over modern .jf files starting with p

# Convert threshold file to associative array (requires Bash 4+)
declare -A cutoff_map
while IFS=$'\t' read -r sample cutoff; do
  cutoff_map["$sample"]=$cutoff
done < <(tail -n +2 "$thresholds_file")  # Skip header

# Loop through each isolate .jf
for jf in "$wd"/isolate_jf/p*.jf; do
    iso_name=$(basename "$jf" .jf)
    echo "Processing $iso_name..."

    # Get cutoff from map
    cutoff=${cutoff_map[$iso_name]}
    if [[ -z "$cutoff" ]]; then
        echo "  ⚠️  Warning: No cutoff found for $iso_name — skipping"
        continue
    fi
    echo "  Using cutoff: $cutoff"
    for ref_kmer in $kmerref/*_unique_filtered_iterative2.txt; do
        ref_name=$(basename "$ref_kmer" _unique_filtered_iterative2.txt)

        # Collect k-mer depths
        counts=()
        while read -r kmer; do
            count=$(jellyfish query "$jf" "$kmer" 2>/dev/null | awk '{print $2}')
            [[ -z "$count" ]] && count=0
            counts+=("$count")
        done < "$ref_kmer"

        total_kmers=${#counts[@]}
        if [[ "$total_kmers" -eq 0 ]]; then
            echo -e "${iso_name}\t${ref_name}\t0\t0\t0.0" >> "$summary_out"
            continue
        fi

        # Count matched kmers with depth ≥ cutoff
        matched_kmers=$(printf "%s\n" "${counts[@]}" | awk -v c="$cutoff" '$1 >= c' | wc -l | awk '{print $1}')
        prop=$(python3 -c "print(round($matched_kmers / $total_kmers, 4))")

        echo -e "${iso_name}\t${ref_name}\t${total_kmers}\t${matched_kmers}\t${prop}" >> "$summary_out"
    done
done

echo "✅ Done! Output written to $summary_out"

