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
#source ~/miniconda3/bin/activate phylogeny_snp

# Extract sequences using samtools
#samtools faidx ../sl_p25.C2_annotated_ref/p25.C2.contigs.second_polished.pilon.annotated.with_Tail_Fibers_haps.fasta -r gene_coordinates.bed > extracted_genes.fasta

# Input files
REF="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/HTFreference/sl_p25.C2_annotated_ref/p25.C2.contigs.second_polished.pilon.annotated.with_Tail_Fibers_haps.fasta"
BED="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/scripts/step2_Oantigengenes/historical49_modern57_bybamcov/gene_coordinates.bed"
OUTPUT="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/historical49_modern57_bybamcov/output.fasta"

# Clear output file if exists
> "$OUTPUT"

# Read BED file line by line
while IFS=$'\t' read -r chrom start end gene strand; do
    if [[ "$strand" == "-" ]]; then
        # Extract reverse complement for negative strand
        #seq=$(samtools faidx --reverse-complement "$REF" "$chrom:$start-$end" | tail -n +2)
        seq=$(samtools faidx "$REF" "$chrom:$start-$end" | tail -n +2)
    else
        # Extract sequence for positive strand
        seq=$(samtools faidx "$REF" "$chrom:$start-$end" | tail -n +2)
    fi

    # Write to output FASTA with proper headers
    echo ">$gene|$chrom:$start-$end" >> "$OUTPUT"
    echo "$seq" >> "$OUTPUT"
done < "$BED"

echo "Sequences saved to $OUTPUT"

#the strands have been double check by samtools faidx manually
