#!/bin/bash
set -euo pipefail

# Define paths
BINARY_FILE="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/historical49_modern57_bybamcov/coverage_h49_m57/../final_fivegenes_binary_matrix.tsv"
ESP_LIST="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/checkespE2_57m/final_afterexpasy/list_present39.txt"
OUT_SIX="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/checkespE2_57m/final_sixgenes_withm57epsE2rescued_binary_matrix.tsv"
# Copy the five-gene matrix to new file
cp "$BINARY_FILE" "$OUT_SIX"

# Read header (get samples)
header=$(head -n1 "$BINARY_FILE")
samples=($(echo "$header" | cut -f2-))

# Load espE2 present samples into a hash (set)
declare -A esp_present
while read -r s; do
    esp_present["$s"]=1
done < "$ESP_LIST"

# Start constructing espE2 row
esp_line="m57_espE2_rescued"

for sample in "${samples[@]}"; do
    if [[ "$sample" != p* && "$sample" != "DC3000" ]]; then
        esp_line+="\tNA"
    elif [[ ${esp_present[$sample]+_} ]]; then
        esp_line+="\t1"
    else
        esp_line+="\t0"
    fi
done

# Check column count
esp_cols=$(echo -e "$esp_line" | awk -F'\t' '{print NF}')
header_cols=$(echo "$header" | awk -F'\t' '{print NF}')
if [[ $esp_cols -ne $header_cols ]]; then
    echo "❌ Column mismatch! espE2 has $esp_cols cols, header has $header_cols"
    echo "espE2 line: $esp_line"
    exit 1
fi

# Append espE2 row
echo -e "$esp_line" >> "$OUT_SIX"
echo "✅ espE2 row added to: $OUT_SIX"
