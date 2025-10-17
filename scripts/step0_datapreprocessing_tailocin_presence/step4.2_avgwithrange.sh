OUT_DIR='/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step4_tailocin_presence'
# -------------------------
# Step 5: Detailed avg (mean–min–max)
# -------------------------
avg_detail_out="$OUT_DIR/avg_stats_tailocin_detailed.txt"

# --- regenerate all_samples.tmp by combining group tables if missing ---
if [ ! -f "$OUT_DIR/all_samples.tmp" ]; then
    echo "♻️  Regenerating all_samples.tmp from H40 + M57 tables..."
    cat "$OUT_DIR/h40_samples.tsv" "$OUT_DIR/m57_samples.tsv" > "$OUT_DIR/all_samples.tmp"
fi

echo -e "group\tcov_prop_mean(min-max)\tdepth_mean(min-max)" > "$avg_detail_out"

compute_detailed_avg () {
    file=$1
    label=$2
    n=$(wc -l < "$file")
    if [[ "$n" -eq 0 ]]; then
        echo -e "$label\tNA\tNA"
    else
        cov_mean=$(awk '{sum+=$2} END {printf "%.3f", sum/NR}' "$file")
        cov_min=$(awk 'NR==1{min=$2} $2<min{min=$2} END{printf "%.3f", min}' "$file")
        cov_max=$(awk 'NR==1{max=$2} $2>max{max=$2} END{printf "%.3f", max}' "$file")

        depth_mean=$(awk '{sum+=$3} END {printf "%.1f", sum/NR}' "$file")
        depth_min=$(awk 'NR==1{min=$3} $3<min{min=$3} END{printf "%.1f", min}' "$file")
        depth_max=$(awk 'NR==1{max=$3} $3>max{max=$3} END{printf "%.1f", max}' "$file")

        printf "%s\t%s (%s–%s)\t%s (%s–%s)\n" \
            "$label" "$cov_mean" "$cov_min" "$cov_max" "$depth_mean" "$depth_min" "$depth_max"
    fi
}

compute_detailed_avg "$OUT_DIR/h40_samples.tsv" "H40" >> "$avg_detail_out"
compute_detailed_avg "$OUT_DIR/m57_samples.tsv" "M57" >> "$avg_detail_out"
compute_detailed_avg "$OUT_DIR/all_samples.tmp" "ALL" >> "$avg_detail_out"

echo "📈 Detailed mean(min–max) summary saved to: $avg_detail_out"

# optional cleanup
mkdir -p "$OUT_DIR/depth"
mv "$OUT_DIR"/*depth.txt "$OUT_DIR/depth" 2>/dev/null || true