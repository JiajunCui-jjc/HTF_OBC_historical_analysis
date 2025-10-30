#!/bin/bash

# Set your working directory
indir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step2_isolate_kmers/isolate_jf/dump"

cd "$indir" || { echo "failed to access $indir"; exit 1; }

# Step 1: Make folder to store original modern isolates
mkdir -p modern_original_big

# Step 2: Move all modern (p*-dump.txt) files to backup folder
mv p*_dump.txt modern_original_big/

# Step 3: Loop through each file and generate 20% sampled output
for file in modern_original_big/p*_dump.txt; do
    [[ -f "$file" ]] || continue
    fname=$(basename "$file")
    outfile="${fname%.txt}_20p_dump.txt"

    # Sample 20% lines randomly using awk
    awk 'BEGIN {srand()} rand() < 0.2' "$file" > "$outfile"

    echo "✅ Sampled 20% from $fname → $outfile"
done
#rm modern_original_big
