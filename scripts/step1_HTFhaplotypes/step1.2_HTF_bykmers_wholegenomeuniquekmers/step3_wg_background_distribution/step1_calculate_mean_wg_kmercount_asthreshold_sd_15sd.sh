#!/usr/bin/env bash
set -euo pipefail
# -------------------------------------------------------------------------
# Calculate per-sample k-mer coverage statistics from Jellyfish dumps.
# This version computes the mean and SD outputs with count >1 both ±1SD and ±0.5SD ranges.
#
# Output:
#   1. wg_kmer_mean.tsv                     → all samples (sample, cutoff)
#   2. threshold_m57.tsv                    → modern samples only
#   3. threshold_h49.tsv                    → historical samples only
#   4. threshold_h49_sd_twosides.tsv        → historical, mean ± SD (after filter)
#   5. threshold_h49_sd_twosides_0.5sd.tsv  → historical, mean ± 0.5×SD (after filter)
# -------------------------------------------------------------------------

# === Input directory ===
dump_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step2_isolate_kmers/isolate_jf/dump/"

# === Output file paths ===
out_combined="${dump_dir}/wg_kmer_mean.tsv"
out_modern="${dump_dir}/threshold_m57.tsv"
out_hist="${dump_dir}/threshold_h49.tsv"
out_hist_sd="${dump_dir}/threshold_h49_sd_twosides.tsv"
out_hist_sd15="${dump_dir}/threshold_h49_sd_twosides_0.5sd.tsv"

# === Headers ===
printf "sample\tcutoff\n" > "$out_combined"
printf "sample\tcutoff\n" > "$out_modern"
printf "sample\tcutoff\n" > "$out_hist"
printf "sample\tmean_minus_sd\tmean_plus_sd\n" > "$out_hist_sd"
printf "sample\tmean_minus_0.5sd\tmean_plus_0.5sd\n" > "$out_hist_sd15"

# Enable nullglob
shopt -s nullglob

# -------------------------------------------------------------------------
# Iterate through all *_dump.txt files (each sample’s k-mer histogram)
# -------------------------------------------------------------------------
for f in "${dump_dir}"/*_dump.txt; do
  [[ -e "$f" ]] || continue

  base=$(basename "$f")
  sample="${base%.*}"
  sample="${sample%_dump}"

  # --- Compute stats with AWK ---
  stats=$(awk '
  NF>=2 {
    c = $2
    if (c >1) { a[++n]=c; sum+=c; sumsq+=c*c }
  }
  END {
    if (n < 1) { print "NA\tNA\tNA\tNA\tNA"; exit }

    # Compute mean and SD using all values (no filtering)
    m  = sum / n
    sd = (n>1) ? sqrt((sumsq - sum*sum/n)/(n-1)) : 0

    # Define cutoffs
    cutoff = m
    mean_minus_sd  = m - sd
    mean_plus_sd   = m + sd
    mean_minus_0_5sd = m - 0.5*sd
    mean_plus_0_5sd  = m + 0.5*sd

    print cutoff "\t" mean_minus_sd "\t" mean_plus_sd "\t" mean_minus_0_5sd "\t" mean_plus_0_5sd
    }
  ' "$f")
  #round to first after dot
  cutoff=$(echo "$stats" | cut -f1 | awk '{printf "%.1f", $1}')
  mean_minus_sd=$(echo "$stats" | cut -f2 | awk '{printf "%.1f", $1}')
  mean_plus_sd=$(echo "$stats" | cut -f3 | awk '{printf "%.1f", $1}')
  mean_minus_0_5sd=$(echo "$stats" | cut -f4 | awk '{printf "%.1f", $1}')
  mean_plus_0_5sd=$(echo "$stats" | cut -f5 | awk '{printf "%.1f", $1}')

  # --- Write outputs ---
  printf "%s\t%s\n" "$sample" "$cutoff" >> "$out_combined"

  if [[ "$sample" =~ ^p ]]; then
    # modern
    printf "%s\t%s\n" "$sample" "$cutoff" >> "$out_modern"
  else
    # historical
    printf "%s\t%s\n" "$sample" "$cutoff" >> "$out_hist"
    printf "%s\t%s\t%s\n" "$sample" "$mean_minus_sd" "$mean_plus_sd" >> "$out_hist_sd"
    printf "%s\t%s\t%s\n" "$sample" "$mean_minus_0_5sd" "$mean_plus_0_5sd" >> "$out_hist_sd15"
  fi
done

# -------------------------------------------------------------------------
# Cleanup for modern sample naming
# -------------------------------------------------------------------------
mv "$dump_dir/threshold_m57.tsv" "$dump_dir/threshold_m57nametmp.tsv"
awk '{sub(/_dump.*/, "", $1); print $1 "\t" $2}' "$dump_dir/threshold_m57nametmp.tsv" \
  > "$dump_dir/threshold_m57.tsv"
rm "$dump_dir/threshold_m57nametmp.tsv"

# -------------------------------------------------------------------------
# Summary
# -------------------------------------------------------------------------
echo "✅ Wrote outputs:"
echo "  - $out_combined"
echo "  - $out_modern"
echo "  - $out_hist"
echo "  - $out_hist_sd"
echo "  - $out_hist_sd15"

