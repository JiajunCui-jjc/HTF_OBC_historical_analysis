#!/bin/bash

# Define directories
RESULTS_DIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/historical49_modern57_bybamcov/coverage_h49_m57"
OUTPUT_TABLE="${RESULTS_DIR}/covprop_matrix.tsv"

# Find all sample TSV files
TSV_FILES=($(ls ${RESULTS_DIR}/*coverage.tsv))

# Extract gene names from the first file (assuming all samples have the same genes)
GENE_FILE="genes.txt"
awk 'NR>1 {print $2}' "${TSV_FILES[0]}" > "$GENE_FILE"

# Create the header row with "Gene" followed by sample names
HEADER="Gene"
for file in "${TSV_FILES[@]}"; do
    sample=$(awk 'NR==2 {print $1}' "$file")  # Extract sample name
    HEADER+="\t${sample}"
done
echo -e "$HEADER" > "$OUTPUT_TABLE"

# Prepare a temporary file to store all depth values
TEMP_FILE="temp_covprops.txt"
paste <(cat "$GENE_FILE") > "$TEMP_FILE"  # Start with gene names

# Extract depth values for each sample and append them column-wise
for file in "${TSV_FILES[@]}"; do
    depths=$(awk 'NR>1 && $2 != "" {print ($4=="NA") ? 0 : $4}' "$file")

    # Ensure correct number of rows before appending
    if [[ $(wc -l < "$GENE_FILE") -eq $(echo "$depths" | wc -l) ]]; then
        paste "$TEMP_FILE" <(echo "$depths") > "temp_paste.tsv" && mv "temp_paste.tsv" "$TEMP_FILE"
    else
        echo "⚠️ WARNING: Mismatch in number of genes and extracted values for sample $sample. Skipping..."
    fi
done

# Final output: Combine the header and data
cat "$TEMP_FILE" >> "$OUTPUT_TABLE"

# Clean up intermediate files
rm -f "$GENE_FILE" "$TEMP_FILE"

echo "✅ Table generated: $OUTPUT_TABLE"
