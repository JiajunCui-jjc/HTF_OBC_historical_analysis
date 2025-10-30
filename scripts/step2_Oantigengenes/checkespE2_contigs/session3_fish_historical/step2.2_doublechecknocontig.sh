#!/bin/bash
set -euo pipefail

# --- Directories ---
base_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/fish_historical/using39m_espE2fasta"
asm_dir="${base_dir}/assemblies"
map_dir="${base_dir}/mappings"
out_file="${map_dir}/nocontigs.txt"

mkdir -p "$map_dir"
> "$out_file"   # clear file

echo "🔍 Checking assemblies in: $asm_dir"
echo "📝 Writing missing contigs to: $out_file"

# --- Loop through sample folders ---
for d in "${asm_dir}"/*; do
    [[ -d "$d" ]] || continue
    sample=$(basename "$d")
    contig="${d}/contigs.fasta"

    if [[ ! -s "$contig" ]]; then
        echo "$sample" >> "$out_file"
        echo "❌ No contigs for: $sample"
    else
        echo "✅ Found contigs for: $sample"
    fi
done

echo "--------------------------------------"
echo "✅ Summary written to: $out_file"
echo "Missing samples count: $(wc -l < "$out_file")"

