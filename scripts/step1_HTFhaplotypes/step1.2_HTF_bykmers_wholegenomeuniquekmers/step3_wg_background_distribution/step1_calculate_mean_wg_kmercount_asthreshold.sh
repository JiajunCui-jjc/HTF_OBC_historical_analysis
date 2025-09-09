#!/usr/bin/env bash
set -euo pipefail
#with count >2 and < mean +3sd
# === INPUT: set your dump directory ===
dump_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step2_isolate_kmers/isolate_jf/dump/"

# === OUTPUT files ===
out_combined="${dump_dir}/wg_kmer_mean.tsv"         # sample \t cutoff
out_modern="${dump_dir}/threshold_m57.tsv"          # starts with 'p'
out_hist="${dump_dir}/threshold_h40.tsv"            # not starting with 'p'

# Header
printf "sample\tcutoff\n" > "$out_combined"
printf "sample\tcutoff\n" > "$out_modern"
printf "sample\tcutoff\n" > "$out_hist"

shopt -s nullglob

for f in "${dump_dir}"/*_dump.txt; do
  [[ -e "$f" ]] || continue

  base=$(basename "$f")
  sample="${base%.*}"
  sample="${sample%_dump}"

  mean=$(awk '
    NF>=2 {
      c = $2
      if (c > 2) { a[++n]=c; sum+=c; sumsq+=c*c }
    }
    END {
      if (n < 1) { print "NA"; exit }
      m  = sum / n
      sd = (n>1) ? sqrt((sumsq - sum*sum/n)/(n-1)) : 0
      thr = m + 3*sd
      sum2 = 0; n2 = 0
      for (i=1; i<=n; i++) if (a[i] <= thr) { sum2 += a[i]; n2++ }
      if (n2 < 1) { print "NA"; exit }
      m2 = sum2 / n2
      base = int(m2)
      print (m2 - base < 0.5) ? base : base + 1
    }
  ' "$f")

  printf "%s\t%s\n" "$sample" "$mean" >> "$out_combined"

  if [[ "$sample" =~ ^p ]]; then
    printf "%s\t%s\n" "$sample" "$mean" >> "$out_modern"
  else
    printf "%s\t%s\n" "$sample" "$mean" >> "$out_hist"
  fi
done

echo "Wrote:"
echo "  - $out_combined"
echo "  - $out_modern"
echo "  - $out_hist"
mv $dump_dir/threshold_m57.tsv $dump_dir/threshold_m57nametmp.tsv

awk '{sub(/_dump.*/, "", $1); print $1 "\t" $2}' "$dump_dir/threshold_m57nametmp.tsv" \
  > "$dump_dir/threshold_m57.tsv"
rm $dump_dir/threshold_m57nametmp.tsv
