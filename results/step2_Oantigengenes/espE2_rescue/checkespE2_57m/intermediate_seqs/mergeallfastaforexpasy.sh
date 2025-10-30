#!/bin/bash
set -euo pipefail

# Define directories and files
WORKDIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/checkespE2_57m"
INTERMEDIATE="$WORKDIR/intermediate_seqs"
ORIG_FASTA="${INTERMEDIATE}/sub_7full_30_retrieved_epsE2_sequences.fasta"
REPLACE_FASTA="${INTERMEDIATE}/extend200bp_p20D4_p6A10.fasta"
FINAL_FASTA="${INTERMEDIATE}/forexpasy_57espE2.fasta"
ABSENT_LIST="$WORKDIR/minimap_results/sub_20_absent_nomatch.txt"


#!/bin/bash
set -euo pipefail

# Define paths
WORKDIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/checkespE2_57m"
EXTRA_DIR="${INTERMEDIATE}/OantigenP_but_espE2notpresent_p22D1_p25C11_p12F2_p24B5_p5C1usetagH1"

TMP1="${INTERMEDIATE}/tmp_filtered.fasta"
TMP2="${INTERMEDIATE}/tmp_merged.fasta"

# Step 1: Remove p20.D4 and p6.A10 from original 30
awk '
BEGIN{RS=">"; ORS=""}
NR>1 {
  split($1, header, /\s+/);
  id = header[1];
  if (id != "p20.D4:125912-130738" && id != "p6.A10:36463-41289" && id !="p24.B5:1-3856" && id != "p5.C1:17-3909") {
    print ">" $0
  }
}' "$ORIG_FASTA" > "$TMP1"

# Step 2: Add extend200bp_p20D4_p6A10.fasta
cat "$TMP1" "$REPLACE_FASTA" > "$TMP2"

# Step 3: Add 4 espE2 sequences with renamed headers (virtual)
declare -A sample_names=(
    ["p12.F2.espE2.fasta"]="p12.F2"
    ["p22.D1.espE2.fasta"]="p22.D1"
    ["p25.C11.espE2half.fasta"]="p25.C11"
    ["p5.C1.espE2.fasta"]="p5.C1"
)

for fasta_file in "${!sample_names[@]}"; do
    sample="${sample_names[$fasta_file]}"
    input="${EXTRA_DIR}/${fasta_file}"

    awk -v new_header="$sample" 'BEGIN {RS=">"; ORS=""} NR > 1 {
        n = split($0, lines, "\n")
        printf(">%s\n", new_header)
        for (i = 2; i <= n; i++) print lines[i]
        print "\n" 
    }' "$input" >> "$TMP2"
done
# Step 4: Add 18 empty >sample entries (excluding p12.F2 and p5.C1)
grep -v -e "p12.F2" -e "p5.C1" -e "p22.D1" -e "p25.C11" "$ABSENT_LIST" | while read -r sample; do
    echo ">$sample" >> "$TMP2"
done
# Add p24.B5 manually (if not already added)
echo ">p24.B5" >> "$TMP2"

# Step 5: Move to final output
mv "$TMP2" "$FINAL_FASTA"

rm -f "$TMP1"

echo "✅ Final FASTA written to: $FINAL_FASTA"
