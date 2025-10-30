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

# Define directories
BAM_DIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_m57/all_bams_m57_OTU5withHTFTFA_dedup/"


OUTPUT_DIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/base_depth/"
mkdir -p "$OUTPUT_DIR"

# Define gene coordinates for epsE_2
GENE_NAME="epsE_2"
CONTIG="30"  # Adjust this if actual contig name is different
GENE_START=142781
GENE_END=147608

# Initialize output file
FINAL_OUTPUT_FILE="${OUTPUT_DIR}/epsE_2_base_depth.tsv"
#echo -e "Sample\tPosition\tRead_Depth" > "$FINAL_OUTPUT_FILE"

# Loop through all 85 samples
for i in $(seq 1 57); do
    sample=$(sed -n "${i}p" /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/modern57.txt)

    # Find the BAM file (handles variable suffixes like _sorted.bam, _mapped.bam, etc.)
    BAM_FILE=$(find "$BAM_DIR" -maxdepth 1 -name "${sample}*.bam" | head -n 1)

    # Check if BAM file exists
    if [[ -z "$BAM_FILE" || ! -f "$BAM_FILE" ]]; then
        echo "ERROR: No BAM file found for sample ${sample}" >&2
        continue
    fi

    echo "Processing BAM file: $BAM_FILE"

    # Ensure BAM file is indexed
    samtools index "$BAM_FILE"

    # Extract per-base depth for epsE_2 gene region and append to final file
    samtools depth -r "${CONTIG}:${GENE_START}-${GENE_END}" "$BAM_FILE" | \
    awk -v sample="$sample" '{print sample, $2, $3}' >> "$FINAL_OUTPUT_FILE"

    echo "Base mapping depth extraction complete for ${sample}!"
done

echo "All samples processed. Final results saved in ${FINAL_OUTPUT_FILE}."

