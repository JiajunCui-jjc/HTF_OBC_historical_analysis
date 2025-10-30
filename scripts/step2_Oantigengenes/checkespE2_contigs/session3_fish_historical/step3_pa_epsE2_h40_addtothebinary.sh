#!/bin/bash
set -euo pipefail

# ================= USER CONFIG =================
output_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/fish_historical/using39m_espE2fasta"

# Define paths
BINARY_FILE="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/checkespE2_57m/final_sixgenes_withm57epsE2rescued_binary_matrix.tsv"
OUT_SIX="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/fish_historical/final_sixgenes_withm57_h49_epsE2rescued_binary_matrix.tsv"


ESPh_nocontig="${output_dir}/mappings/nocontigs.txt"
#ESPh_nocontig="${output_dir}/mappings/nocontigs_dotname.txt"

# Clear or create output file
#> "$ESPh_nocontig"

# Replace the first underscore "_" with a dot "." in each line
#while read -r line; do
#    converted=$(echo "$line" | sed 's/_/./')
#    echo "$converted" >> "$ESPh_nocontig"
#done < "$ESPh_nocontig_original"


# ================= READ HEADER =================
header=$(head -n1 "$BINARY_FILE")
samples=($(echo "$header" | cut -f2-))
echo "$header" > "$OUT_SIX"

# ================= LOAD ABSENT SAMPLES =================
declare -A esp_absent
while read -r s; do
    [[ -n "$s" ]] && esp_absent["$s"]=1
done < "$ESPh_nocontig"

# ================= PROCESS EACH ROW =================
tail -n +2 "$BINARY_FILE" | while IFS=$'\t' read -r gene_name rest; do
    if [[ "$gene_name" == "m57_espE2_rescued" ]]; then
        echo " Updating row: $gene_name → m57_h36_espE2_rescued"

        IFS=$'\t' read -r -a values <<< "$rest"
        new_line="espE2"

        for i in "${!samples[@]}"; do
            sample="${samples[$i]}"
            if [[ "$sample" != p* && "$sample" != "DC3000" ]]; then
                if [[ ${esp_absent[$sample]+_} ]]; then
                    new_line+="\t0"
                else
                    new_line+="\t1"
                fi
            else
                new_line+="\t${values[$i]}"
            fi
        done

        echo -e "$new_line" >> "$OUT_SIX"
    else
        echo -e "$gene_name\t$rest" >> "$OUT_SIX"
    fi
done

echo "final matrix written to: $OUT_SIX"
