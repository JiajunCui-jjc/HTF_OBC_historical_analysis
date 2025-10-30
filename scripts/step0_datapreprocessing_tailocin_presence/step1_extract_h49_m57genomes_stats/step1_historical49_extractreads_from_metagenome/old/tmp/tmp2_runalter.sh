#!/bin/bash -l
#$ -l tmem=6G
#$ -l h_vmem=6G
#$ -l h_rt=2:00:0
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -N regen_cov_fix_At_block

source ~/miniconda3/bin/activate phylogeny_snp

# ===========================================================
# 1. Define paths
# ===========================================================
baseAt=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/2025_h49/h10_lopez
answersdir=/SAN/ugi/plant_genom/jiajucui/answer2025_h49/h10_lopez

# ===========================================================
# 2. Loop over all A.thaliana dedup BAMs
# ===========================================================
for bam in ${baseAt}/*_mapped_At_q20.rmdup.sorted.bam; do
    samplename=$(basename "$bam" | sed 's/_mapped_At_q20.rmdup.sorted.bam//')
    ansfile=${answersdir}/answers_for_${samplename}.txt
    covfile=${baseAt}/${samplename}_mapped_At_q20.rmdup.sorted.bysamtoolscoverage.txt

    # Skip if no corresponding answer file
    [[ ! -f "$ansfile" ]] && { echo "⚠️ No answer file for $samplename, skipping."; continue; }

    echo "🔧 Processing $samplename"

    # ===========================================================
    # Step 1. (Re)generate samtools coverage file
    # ===========================================================
    echo "  → Generating coverage file..."
    samtools coverage "$bam" > "$covfile"

    # ===========================================================
    # Step 2. Compute total and mapped bases
    # ===========================================================
    total_len=$(awk '{sum+=$3} END{print sum}' "$covfile")
    mapped_len=$(awk '{sum+=$5} END{print sum}' "$covfile")

    # ===========================================================
    # Step 3. Replace At length block inside the answer file
    # ===========================================================
    tmpfile=$(mktemp)
    awk -v total="$total_len" -v mapped="$mapped_len" '
        BEGIN {replaced=0}
        /^total length of At ref:/ {
            print "total length of At ref:"
            print total
            print "base mapped to At ref:"
            print mapped
            getline  # skip the old numeric line
            replaced=1
            next
        }
        {print}
    ' "$ansfile" > "$tmpfile" && mv "$tmpfile" "$ansfile"

    echo "✅ Updated $samplename"
done

echo "🎯 All coverage files regenerated and At stats updated successfully."
