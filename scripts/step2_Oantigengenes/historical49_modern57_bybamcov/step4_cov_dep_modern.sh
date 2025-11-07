#!/bin/bash -l

#$ -l tmem=4G  # Allocate 4GB memory per core
#$ -l h_vmem=4G
#$ -l h_rt=8:00:0
#$ -wd /SAN/ugi/plant_genom/jiajucui/
#$ -V
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/

# Input: sample index from job array
#i=$SGE_TASK_ID
source ~/miniconda3/bin/activate phylogeny_snp
# Extract sample name from sample list
for i in $(seq 1 57); do
sample=$(sed -n "${i}p" /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/modern57.txt)
#40h_OTU5.txt  all_fasta_m57  all_fastq_h40  all_fastq_m57  archives  HTFreference  LPScluster_ref  modern57.txt
#/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data
# Define BAM directory (corrected path)
BAM_DIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_m57/all_bams_m57_OTU5withHTFTFA_dedup/"
#p12.A11_mapped_Ps_q20.dedup.sorted.bam 

# Find the BAM file (handles variable suffixes like _sorted.bam, _mapped.bam, etc.)
BAM_FILE=$(find "$BAM_DIR" -maxdepth 1 -name "${sample}*.bam" | head -n 1)

# Check if BAM file exists
if [[ -z "$BAM_FILE" || ! -f "$BAM_FILE" ]]; then
    echo "ERROR: No BAM file found for sample ${sample}" >&2
    exit 1
fi

echo "Processing BAM file: $BAM_FILE"

samtools index $BAM_FILE
# Define output directory and file
OUTPUT_DIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/historical49_modern57_bybamcov/"
OUTPUT_FILE="${OUTPUT_DIR}/coverage_h49_m57/${sample}_gene_coverage.tsv"


#GENE_COORDS="${OUTPUT_DIR}/after_step2_coordinates.bed"
GENE_COORDS="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/scripts/step2_Oantigengenes/historical49_modern57_bybamcov/gene_coordinates.bed"

# Header for output file
    echo -e "Sample\tGene\tAvg_Coverage\tCovered_Proportion\tRelative_Read_Depth\tTotal_Ref_Depth" > "$OUTPUT_FILE"

    # **Step 1: Compute total reference depth**
        total_ref_depth=$(samtools depth -aa "$BAM_FILE" | awk -v total=5941411 '{sum += $3} END {print sum / total}')

    # **Step 2: Process each gene**
    while read chrom start end gene strand; do
        [[ "$chrom" == "#"* ]] && continue  # Skip header lines

        # Compute per-base depth using samtools depth
        samtools depth -aa -r "$chrom:$start-$end" "$BAM_FILE" | awk -v sample="$sample" -v gene="$gene" -v total_ref_depth="$total_ref_depth" '
        {
            total_depth += $3;  # Sum depth values
            if ($3 > 0) covered++;  # Count bases with at least 1× coverage
            count++;  # Total bases in gene region
        }
        END {
            if (count > 0) {
                avg_cov = total_depth / count;  # Compute mean depth
                prop_cov = (covered / count) * 100;  # Compute covered proportion (%)
                relative_depth = (total_ref_depth > 0) ? total_depth / (total_ref_depth * count) : "NA";  # Avoid division by zero
                print sample, gene, avg_cov, prop_cov, relative_depth, total_ref_depth;
            } else {
                print sample, gene, 0, 0, "NA", total_ref_depth;  # If no coverage, output "NA" for relative depth
            }
        }' >> "$OUTPUT_FILE"

    done < "$GENE_COORDS"

    echo "Coverage analysis complete for ${sample}! Results saved in ${OUTPUT_FILE}."
done
