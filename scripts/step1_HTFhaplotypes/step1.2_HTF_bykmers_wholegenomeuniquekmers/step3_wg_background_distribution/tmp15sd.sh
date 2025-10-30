#!/usr/bin/env bash
set -euo pipefail
# -------------------------------------------------------------------------
# Compute historical-sample k-mer coverage thresholds using ±1.5 SD
# after filtering kmers > mean + 3·SD.
#
# Input:  Jellyfish *_dump.txt files
# Output: threshold_h49_sd_twosides_1.5sd.tsv
# -------------------------------------------------------------------------

# === Input directory ===
dump_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step2_isolate_kmers/isolate_jf/dump/"

# === Output ===
out_hist_sd15="${dump_dir}/threshold_h49_sd_twosides_1.5sd.tsv"
printf "sample\tmean_minus_1.5sd\tmean_plus_1.5sd\n" > "$out_hist_sd15"

# Enable nullglob to skip if no files
shopt -s nullglob

for f in "${dump_dir}"/*_dump.txt; do
  [[ -e "$f" ]] || continue

  base=$(basename "$f")
  sample="${base%.*}"
  sample="${sample%_dump}"
  # Skip modern samples (those starting with lowercase p)
  if [[ "$sample" =~ ^p ]]; then
    continue
  fi
  # --- Compute mean ± 1.5 SD (after filter) ---
  stats=$(awk '
    NF>=2 {
      c = $2
      if (c > 2) { a[++n]=c; sum+=c; sumsq+=c*c }
    }
    END {
      if (n < 1) { print "NA\tNA"; exit }

      # 1st pass: mean & SD
      m  = sum / n
      sd = (n>1)?sqrt((sumsq - sum*sum/n)/(n-1)):0
      thr = m + 3*sd

      # Filter kmers within mean + 3SD
      sum2=0; sumsq2=0; n2=0
      for (i=1; i<=n; i++) {
        if (a[i] <= thr) {
          sum2 += a[i]; sumsq2 += a[i]^2; n2++
        }
      }
      if (n2 < 1) { print "NA\tNA"; exit }

      # Recalculate mean & SD on filtered data
      m2  = sum2 / n2
      sd2 = (n2>1)?sqrt((sumsq2 - sum2*sum2/n2)/(n2-1)):0

      minus15 = m2 - 1.5*sd2
      plus15  = m2 + 1.5*sd2
      print minus15 "\t" plus15
    }
  ' "$f")

  mean_minus_1_5sd=$(echo "$stats" | cut -f1)
  mean_plus_1_5sd=$(echo "$stats" | cut -f2)

  printf "%s\t%s\t%s\n" "$sample" "$mean_minus_1_5sd" "$mean_plus_1_5sd" >> "$out_hist_sd15"
done

echo "✅ Done. Wrote: $out_hist_sd15"

