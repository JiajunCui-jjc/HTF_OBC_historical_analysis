#!/bin/bash
# -------------------------------------------------------------
# Use depth >= (sample-specific mean - 0.25xsd)
#      covprop >= 65
# to assign presence = 1 else 0
# -------------------------------------------------------------
#
        #1.     Concatenates the depth and covprop files.
        #2.     Adds suffixes _depth and _covprop to gene rows.
        #3.     Filters out unwanted genes (espE_2, tagG_2, tagH_2).
# -------------------------------------------------------------



set -euo pipefail
# -------------------------------------------------------------

RESULTS_DIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/historical49_modern57_bybamcov/coverage_h49_m57"
DEPTH_FILE="${RESULTS_DIR}/relative_depth_matrix.tsv"
THRESH_FILE="${RESULTS_DIR}/depth_thresholds_mean_minus_025_sd.tsv"
COV_FILE="${RESULTS_DIR}/covprop_matrix.tsv"
CAT_FILE="${RESULTS_DIR}/../final_cov_depth_wide_cat.tsv"
OUT_FILE="${RESULTS_DIR}/../final_fivegenes_binary_matrix.tsv"

# -------------------------------------------------------------
# 1. Combine depth and covprop
# -------------------------------------------------------------
echo "🔹 Combining depth and covprop matrices..."
head -n1 "$DEPTH_FILE" > "$CAT_FILE"
tail -n +2 "$DEPTH_FILE" | awk -F'\t' -v OFS='\t' '{ $1=$1 "_depth"; print }' >> "$CAT_FILE"
tail -n +2 "$COV_FILE"   | awk -F'\t' -v OFS='\t' '{ $1=$1 "_covprop"; print }' >> "$CAT_FILE"
echo "✅ Combined: $CAT_FILE"

# -------------------------------------------------------------
# 2. Load thresholds
# -------------------------------------------------------------
#Sample  Mean_total_depth        SD_total_depth  Mean_minus_SD   relative_threshold

declare -A THR
while read -r sample mean sd minus thr ; do
    [[ "$sample" == "Sample" ]] && continue
    THR["$sample"]="$thr"
done < "$THRESH_FILE"
echo "✅ Loaded thresholds for ${#THR[@]} samples"

# -------------------------------------------------------------
# 3. Parse header
# -------------------------------------------------------------
mapfile -t lines < "$CAT_FILE"
header="${lines[0]}"
echo "$header" > "$OUT_FILE"

IFS=$'\t' read -r -a header_fields <<< "$header"
SAMPLES=("${header_fields[@]:1}")   # skip "Gene"
echo "✅ Parsed ${#SAMPLES[@]} sample columns"

# -------------------------------------------------------------
# 4. Separate depth/covprop lines
# -------------------------------------------------------------
declare -A DEPTH_LINE COV_LINE
for line in "${lines[@]:1}"; do
    gene=$(echo "$line" | cut -f1)
    case "$gene" in *tagG_2*|*tagH_2*|*epsE*) continue ;; esac
    if [[ "$gene" == *_depth ]]; then
        DEPTH_LINE["$gene"]="$line"
    elif [[ "$gene" == *_covprop ]]; then
        COV_LINE["$gene"]="$line"
    fi
done

# -------------------------------------------------------------
# 5. Build binary matrix
# -------------------------------------------------------------
for dkey in "${!DEPTH_LINE[@]}"; do
    gbase="${dkey%_depth}"
    ckey="${gbase}_covprop"
    [[ -z "${COV_LINE[$ckey]+set}" ]] && continue

    # --- declare arrays before reading ---
    declare -a DVALS=()
    declare -a CVALS=()

    # split into arrays (skip first column)
    IFS=$'\t' read -r _gene_d_line <<< "${DEPTH_LINE[$dkey]}"
    IFS=$'\t' read -r _gene_c_line <<< "${COV_LINE[$ckey]}"
    DVALS=($(echo "$_gene_d_line" | cut -f2-))
    CVALS=($(echo "$_gene_c_line" | cut -f2-))

    out=("$gbase")

    for ((i=0; i<${#SAMPLES[@]}; i++)); do
        sample="${SAMPLES[$i]}"
        thr="${THR[$sample]:-0.75}"
        d="${DVALS[$i]:-0}"
        c="${CVALS[$i]:-0}"
#if thr NA then set 0.75 but there is no NA
        if [[ "$d" =~ ^[0-9.]+$ && "$c" =~ ^[0-9.]+$ ]]; then
            if (( $(echo "$d >= $thr && $c >= 65" | bc -l) )); then
                out+=(1)
            else
                out+=(0)
            fi
        else
            out+=(0)
        fi
    done

    printf "%s\n" "$(IFS=$'\t'; echo "${out[*]}")" >> "$OUT_FILE"
done

echo "✅ Binary matrix written to: $OUT_FILE"
head -n 5 "$OUT_FILE"

