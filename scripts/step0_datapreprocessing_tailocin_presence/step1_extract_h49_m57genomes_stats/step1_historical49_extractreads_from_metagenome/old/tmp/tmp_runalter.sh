#!/bin/bash -l
#$ -l tmem=6G
#$ -l h_vmem=6G
#$ -l h_rt=3:00:0
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -N insert_At_stats_cleanup

source ~/miniconda3/bin/activate phylogeny_snp

# ===========================================================
# 1. Paths and constants
# ===========================================================
baseAt=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/2025_h49/h10_lopez
answerspath=/SAN/ugi/plant_genom/jiajucui/answer2025_h49/h10_lopez
lenAt=119667750

# ===========================================================
# 2. Loop over all A.thaliana BAMs
# ===========================================================
for bam in ${baseAt}/*_mapped_At_q20.sorted.bam; do
    samplename=$(basename "$bam" | sed 's/_mapped_At_q20.sorted.bam//')
    out_answer=${answerspath}/answers_for_${samplename}.txt

    echo "▶ Processing ${samplename}..."

    # temporary file for inserting At results
    tmpfile=$(mktemp)

    # ===========================================================
    # --- q1.1 + q1.2 stats generation ---
    # ===========================================================
    stats_block=$(mktemp)
    {
        samtools flagstat ${baseAt}/${samplename}_mapped_At_q20.sorted.bam > ${baseAt}/${samplename}_At_flagstats.q20.log
        echo "reads mapped to At ref:"
        grep 'mapped (' ${baseAt}/${samplename}_At_flagstats.q20.log

        echo "q1.1 What is the covered genome proportion of A.thaliana DNA?"
        samtools coverage ${baseAt}/${samplename}_mapped_At_q20.rmdup.sorted.bam > ${baseAt}/${samplename}_mapped_At_q20.rmdup.sorted.bysamtoolscoverage.txt

        echo "total length of At ref:"
        awk '{sum+=$3;} END{print sum;}' ${baseAt}/${samplename}_mapped_At_q20.rmdup.sorted.bysamtoolscoverage.txt

        echo "base mapped to At ref:"
        awk '{sum+=$5;} END{print sum;}' ${baseAt}/${samplename}_mapped_At_q20.rmdup.sorted.bysamtoolscoverage.txt
        echo "q1.2 What is the read depth of A.thaliana DNA?"
	echo "read depth to At ref:"
        awk -v len=$lenAt '{sum+=($3*$7);} END{print sum/len;}' ${baseAt}/${samplename}_mapped_At_q20.rmdup.sorted.bysamtoolscoverage.txt
    } > "$stats_block"

    # ===========================================================
    # --- Insert the At block between q1 and q2 ---
    # ===========================================================
    awk -v block="$(cat "$stats_block")" '
        BEGIN {inserted=0}
        /q2 / && inserted==0 {print block; inserted=1}
        {print}
    ' "$out_answer" > "$tmpfile" && mv "$tmpfile" "$out_answer"

    # ===========================================================
    # --- Cleanup ---
    # ===========================================================
    find "$baseAt" -type f \( \
        -name "${samplename}.collapsed_At.sorted.bam*" -o \
        -name "${samplename}.collapsed.sai" -o \
        -name "${samplename}.combined_At.bam" -o \
        -name "${samplename}.combined_At.sorted.bam*" -o \
        -name "${samplename}.r1r2_At.sorted.bam*" -o \
        -name "${samplename}.r1.sai" -o \
        -name "${samplename}.r2.sai" \
      \) -delete

    rm -f ${baseAt}/${samplename}_At_flagstats.q20.log \
          ${baseAt}/${samplename}_mapped_At_q20.rmdup.sorted.bysamtoolscoverage.txt \
          "$stats_block"

    echo "✅ Inserted At stats + cleaned for ${samplename}"
done

echo "🎯 All samples processed successfully."
