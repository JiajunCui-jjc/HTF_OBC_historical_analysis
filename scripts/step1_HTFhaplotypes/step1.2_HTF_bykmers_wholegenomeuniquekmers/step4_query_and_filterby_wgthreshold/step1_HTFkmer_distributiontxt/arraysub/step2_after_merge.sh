# --- Paths ---
out=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers/kmer_distribution_mix

# Create a temporary folder for per-sample files
mkdir -p "$out/tmp"

# Move all individual outputs (*.tsv) into tmp/
mv "$out"/*_HTF_kmer_depth.tsv "$out/tmp"

# --- Combine all TSVs into a single master file ---
# 1️⃣ Take the header (first line) from the first file
head -n 1 $(ls "$out/tmp"/*_HTF_kmer_depth.tsv | head -n 1) > "$out/HTF_kmer_depth_distributions.tsv"

# 2️⃣ Append the rest (skip headers from subsequent files)
#    -q : quiet (no file name prefixes)
#    -n +2 : skip first line from each file
tail -q -n +2 "$out/tmp"/*_HTF_kmer_depth.tsv >> "$out/HTF_kmer_depth_distributions.tsv"

# Optional: confirm
wc -l "$out/HTF_kmer_depth_distributions.tsv"
echo "✅ Merged HTF k-mer distributions written to:"
echo "   $out/HTF_kmer_depth_distributions.tsv"

