#!/bin/bash -l
#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=5:00:0
#$ -N h49query_1sd
#$ -V
#$ -cwd
#$ -e /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers/logs
#$ -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers/logs

source ~/miniconda3/bin/activate phylogeny_snp
set -euo pipefail

# --- Paths ---
wd="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step2_isolate_kmers"
out="$wd/../step3_query_kmers"
kmerref="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/kmers_unique_hamming2"
dump_dir="$wd/isolate_jf/dump"
mkdir -p "$out"

# --- Input thresholds ---
thresholds_file="${dump_dir}/threshold_h49_sd_twosides_0.5sd.tsv"
#>=mean-0.5sd
#PL0027 no kmer mapped to any HTF kmer
summary_out="${out}/summaryhistorical49_tailocin_05sd.tsv"

echo "[5] Counting isolate kmers using two-sided thresholds (mean-SD)"
echo -e "Isolate\tReference\tTotal_Allele_Kmers\tMatched_Kmers\tProp" > "$summary_out"

# --- Read thresholds into two associative arrays ---
declare -A cutoff_low_map
declare -A cutoff_high_map
while read -r sample low high; do
    # skip header or empty lines
    [[ "$sample" == "sample" || -z "$sample" ]] && continue
    cutoff_low_map["$sample"]="$low"
    cutoff_high_map["$sample"]="$high"
done < "$thresholds_file"

# --- Loop over historical isolates (.jf not starting with p) ---
for jf in "$wd"/isolate_jf/*.jf; do
    iso_name=$(basename "$jf" .jf)
    [[ "$iso_name" == p* ]] && continue

    low_cut=${cutoff_low_map[$iso_name]:-NA}
    high_cut=${cutoff_high_map[$iso_name]:-NA}

    if [[ "$low_cut" == "NA" || "$high_cut" == "NA" ]]; then
        echo "Warning: thresholds not found for $iso_name — skipping"
        continue
    fi
    echo "Processing $iso_name  (range: >=$low_cut)"

    for ref_kmer in "$kmerref"/*_unique_filtered_iterative2.txt; do
        ref_name=$(basename "$ref_kmer" _unique_filtered_iterative2.txt)

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

        # Count k-mers within [low_cut, high_cut]
        matched_kmers=$(printf "%s\n" "${counts[@]}" | awk -v low="$low_cut" -v high="$high_cut" '$1>=low ' | wc -l | awk '{print $1}')
        prop=$(python3 -c "print(round($matched_kmers / $total_kmers, 4))")

        echo -e "${iso_name}\t${ref_name}\t${total_kmers}\t${matched_kmers}\t${prop}" >> "$summary_out"
    done
done

echo "✅ Done! Output written to $summary_out"

