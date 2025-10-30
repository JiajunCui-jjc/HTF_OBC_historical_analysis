#!/bin/bash -l

#$ -l tmem=4G  # Allocate 4GB memory per core
#$ -l h_vmem=4G
#$ -l h_rt=8:00:0
#$ -wd /SAN/ugi/plant_genom/jiajucui/
#$ -V
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/

# Activate conda environment
source ~/miniconda3/bin/activate phylogeny_snp

# Define BAM directory
BAM_DIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h49/all_raw_meta_fastq_h49/../Ps_bams_maptoOTU5_with_haplotype_h49_dedup/h49_bams_softlink/"

# Define output directory
OUTPUT_DIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/base_depth/"
mkdir -p "$OUTPUT_DIR"

# Define gene coordinates for epsE_2
GENE_NAME="epsE_2"
CONTIG="30"  # Replace with actual contig name if different
GENE_START=142781
GENE_END=147608

# Output file
OUTPUT_FILE="${OUTPUT_DIR}/epsE_2_base_depth.tsv"
echo -e "Sample\tPosition\tRead_Depth" > "$OUTPUT_FILE"

# Loop through all 46 samples
for i in $(seq 1 46); do
    sample=$(sed -n "${i}p" /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/historical49.txt)

    # Find the BAM file
    BAM_FILE=$(find "$BAM_DIR" -maxdepth 1 -name "${sample}*.bam" | head -n 1)

    # Check if BAM file exists
    if [[ -z "$BAM_FILE" || ! -f "$BAM_FILE" ]]; then
        echo "ERROR: No BAM file found for sample ${sample}" >&2
        continue
    fi

    echo "Processing sample: $sample"

    # Extract per-base depth for the gene region
    samtools depth -r "${CONTIG}:${GENE_START}-${GENE_END}" "$BAM_FILE" | \
    awk -v sample="$sample" '{print sample, $2, $3}' >> "$OUTPUT_FILE"
done

echo "Base mapping depth extraction complete. Results saved in ${OUTPUT_FILE}."

