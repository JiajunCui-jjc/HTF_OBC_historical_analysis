#!/bin/bash
#since we show this with mapping stats so we use the otu5only ref to keep them in consistent
# Input files
OTU5_ref="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/HTFreference/OTU5_ref_noHTFhaplotypes/Pseudomonas.OTU5_ref.fasta"
REF="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/HTFreference/tailocin_region.fa"

h43bam_READ_DIR=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h49/Ps_bams_maptoOTU5withouthaplotype_h49_dedup/h43_bams_softlink
format=".mapped_to_Pseudomonas.dd.q20.markeddup.bam"
#.mapped_to_Pseudomonas.dd.q20.markeddup.bam
m57bam="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_m57/all_bams_m57_OTU5refonly_dedup/"
#p4.E6_mapped_Ps_q20.dedup.sorted.bam
formatm57="_mapped_Ps_q20.dedup.sorted.bam"

OUT_DIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step4_tailocin_presence"
THREADS=4
mkdir -p /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence "$OUT_DIR"
cp $REF $OUT_DIR
# -------------------------
# Step 0: Index BAMs
# -------------------------

echo "📂 Indexing H43 BAMs..."
for bam in "$h43bam_READ_DIR"/*"$format"; do
    if [ ! -f "$bam.bai" ]; then
        echo "🔧 Indexing $(basename "$bam")"
        samtools index -@ $THREADS "$bam"
    fi
done

echo "📂 Indexing M57 BAMs..."
for bam in "$m57bam"/*"$formatm57"; do
    if [ ! -f "$bam.bai" ]; then
        echo "🔧 Indexing $(basename "$bam")"
        samtools index -@ $THREADS "$bam"
    fi
done
# -------------------------
# Step 1: Extract region for manual check
# -------------------------
region_chr="22"
region_start=15502
region_end=33558

region_length=$((region_end - region_start + 1))

echo -e "$region_chr\t$region_start\t$region_end" > "$OUT_DIR/region.bed"

# Extract sequence from OTU5 reference for manual validation
seqtk subseq "$OTU5_ref" <(echo -e "$region_chr\t$region_start\t$region_end") > "$OUT_DIR/extracted_tailocin_region.fa"
echo "✅ Extracted tailocin region sequence for manual comparison: $OUT_DIR/extracted_tailocin_region.fa"

# -------------------------
# Step 2: Coverage analysis
# -------------------------

# Output file header
output="$OUT_DIR/tailocin_coverage_summary.tsv"
echo -e "sample\ttailocin_cov_prop\ttailocin_depth" > "$output"

# Function to calculate coverage
process_bam () {
    bamfile=$1
    sample=$2

    # Output depth file
    depth_file="$OUT_DIR/${sample}.depth.txt"

    # Calculate depth in region
    samtools depth -r ${region_chr}:${region_start}-${region_end} "$bamfile" > "$depth_file"

    # Number of positions covered (>0 depth)
    covered_bp=$(awk '$3 > 0 {count++} END {print count+0}' "$depth_file")

    # Total depth across region
    total_depth=$(awk '{sum += $3} END {print sum+0}' "$depth_file")

    # Calculate
    cov_prop=$(awk -v c=$covered_bp -v len=$region_length 'BEGIN{printf "%.4f", c/len}')
    avg_depth=$(awk -v d=$total_depth -v len=$region_length 'BEGIN{printf "%.2f", d/len}')

    echo -e "${sample}\t${cov_prop}\t${avg_depth}" >> "$output"
}

# -------- H43 samples --------
for bam in "$h43bam_READ_DIR"/*"$format"; do
    sample=$(basename "$bam" "$format")
    echo "📘 Processing H43 sample: $sample"
    process_bam "$bam" "$sample"
done

# -------- M57 samples --------
for bam in "$m57bam"/*"$formatm57"; do
    sample=$(basename "$bam" "$formatm57")
    echo "📙 Processing M57 sample: $sample"
    process_bam "$bam" "$sample"
done

echo "✅ Done. Output saved to: $output"

# -------------------------
# Step 4: Calculate average stats
# -------------------------

summary_out="$OUT_DIR/avg_stats_tailocin.txt"

# Split into temp files for h43 and m57
grep -v "^sample" "$output" > "$OUT_DIR/all_samples.tmp"

grep -Ff <(ls "$h43bam_READ_DIR"/*"$format" | xargs -n1 basename | sed "s/$format//") "$OUT_DIR/all_samples.tmp" > "$OUT_DIR/h43_samples.tsv"
grep -Ff <(ls "$m57bam"/*"$formatm57" | xargs -n1 basename | sed "s/$formatm57//") "$OUT_DIR/all_samples.tmp" > "$OUT_DIR/m57_samples.tsv"

# Function to compute mean for a file
compute_avg () {
    file=$1
    label=$2
    n=$(wc -l < "$file")
    if [[ "$n" -eq 0 ]]; then
        echo -e "$label\tNA\tNA"
    else
        avg_cov=$(awk '{sum+=$2} END {printf "%.4f", sum/NR}' "$file")
        avg_depth=$(awk '{sum+=$3} END {printf "%.2f", sum/NR}' "$file")
        echo -e "$label\t$avg_cov\t$avg_depth"
    fi
}

{
    echo -e "group\tavg_tailocin_cov_prop\tavg_tailocin_depth"
    compute_avg "$OUT_DIR/h43_samples.tsv" "H43"
    compute_avg "$OUT_DIR/m57_samples.tsv" "M57"
    compute_avg "$OUT_DIR/all_samples.tmp" "ALL"
} > "$summary_out"

# Clean up
rm "$OUT_DIR/all_samples.tmp"
mkdir -p $OUT_DIR/depth
mv $OUT_DIR/*depth.txt $OUT_DIR/depth
echo "📊 Summary stats saved to: $summary_out"
