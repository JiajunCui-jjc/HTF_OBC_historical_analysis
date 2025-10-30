#!/bin/bash
set -euo pipefail

# =============================================================
# Build relative_depth_matrix.tsv (Gene × Sample)
# using column 5 (Relative_Read_Depth) from *_gene_coverage.tsv
# =============================================================

RESULTS_DIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/historical49_modern57_bybamcov/coverage_h49_m57"
OUT_FILE="${RESULTS_DIR}/relative_depth_matrix.tsv"

TSV_FILES=($(ls ${RESULTS_DIR}/*_gene_coverage.tsv 2>/dev/null))
if [[ ${#TSV_FILES[@]} -eq 0 ]]; then
    echo "No *_gene_coverage.tsv files found in $RESULTS_DIR"
    exit 1
fi

echo "Building relative_depth_matrix.tsv ..."

GENE_FILE=$(mktemp)
awk 'NR>1 {print $2}' "${TSV_FILES[0]}" > "$GENE_FILE"

# header
HEADER="Gene"
for f in "${TSV_FILES[@]}"; do
    sample=$(awk 'NR==2 {print $1}' "$f")
    HEADER+="\t${sample}"
done
echo -e "$HEADER" > "$OUT_FILE"

TEMP_FILE=$(mktemp)
cp "$GENE_FILE" "$TEMP_FILE"

for f in "${TSV_FILES[@]}"; do
    depths=$(awk 'NR>1 {print ($5=="NA")?0:$5}' "$f")   # col5 = Relative_Read_Depth
    if [[ $(wc -l < "$GENE_FILE") -eq $(echo "$depths" | wc -l) ]]; then
        paste "$TEMP_FILE" <(echo "$depths") > tmp && mv tmp "$TEMP_FILE"
    else
        echo "Skipping mismatched file: $f"
    fi
done

cat "$TEMP_FILE" >> "$OUT_FILE"
rm -f "$TEMP_FILE" "$GENE_FILE"

echo "Relative depth matrix written: $OUT_FILE"

