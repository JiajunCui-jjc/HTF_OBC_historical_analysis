#!/bin/bash

# Load environment
source /home/jiajucui/miniconda3/bin/activate phylogeny_snp

# Define directories

# Define directories
RESULT_DIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/checkespE2_57m"
SAMPLE_LIST="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fasta_m57/modern57.txt"
ASSEMBLY_DIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fasta_m57"
MINIMAP_OUTPUT_DIR="${RESULT_DIR}/minimap_results"
QUERY="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/scripts/step2_Oantigengenes/checkespE2_contigs/espE2_p25c2.fasta"
TOOLS="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/shfortailocin/tools"

#the partially matched only
#RETRIEVE_PAF="${RESULT_DIR}/minimap_results/sub_30_retrievepaf.txt"
#the full 7 full + 30 partial
RETRIEVE_PAF="${RESULT_DIR}/minimap_results/sumpaf.txt"


#OUTPUT_FASTA="${RESULT_DIR}/sub_30_retrieved_epsE2_sequences.fasta"
#OUTPUT_FASTA2="${RESULT_DIR}/sub_30_retrieved_epsE2_sequences_withp25.c2.fasta"
mkdir -p ${RESULT_DIR}/intermediate_seqs/
OUTPUT_FASTA="${RESULT_DIR}/intermediate_seqs/sub_7full_30_retrieved_epsE2_sequences.fasta"
#OUTPUT_FASTA2="${RESULT_DIR}/intermediate_seqs/sub_7full_30_retrieved_epsE2_sequences_withp25.c2.fasta"


# Ensure output file is empty
> "$OUTPUT_FASTA"

# Process each line of the PAF file
while read -r line; do
    # Extract relevant fields
    sample=$(echo "$line" | awk '{print $1}')        # Sample name
    # Skip p26.B7
    if [[ "$sample" == "p26.B7" ]]; then
        continue
    fi
    strand=$(echo "$line" | awk '{print $6}')        # Strand orientation
    contig=$(echo "$line" | awk '{print $7}')        # Reference contig
    Q_start=$(echo "$line" | awk '{print $9}')       # Query start position
    Q_end=$(echo "$line" | awk '{print $10}')         # Query end position
    A_len=$(echo "$line" | awk '{print $5}')         # Alignment length

    full_length=4827  # Expected extraction length
    fasta_gz="${ASSEMBLY_DIR}/${sample}.fasta.gz"
    fasta="${ASSEMBLY_DIR}/${sample}.fasta"
    echo $strand
    # Check if the FASTA file exists
    if [[ ! -f "$fasta_gz" ]]; then
        echo "Warning: FASTA file not found for $sample"
        continue
    fi

    # Decompress FASTA if needed
    if [[ ! -f "$fasta" ]]; then
        gunzip -c "$fasta_gz" > "$fasta"
    fi

    # Compute extraction coordinates based on strand
    if [[ "$strand" == "+" ]]; then
        extract_start=$((Q_start + 1))
        extract_end=$((Q_end + (full_length - A_len)))
    else
        extract_start=$((Q_start - (full_length - A_len) + 1))
        extract_end=$Q_end
    fi
    # Extract sequence using samtools faidx and reverse complement if needed
    if [[ "$strand" == "+" ]]; then
        samtools faidx "$fasta" "$contig:$extract_start-$extract_end" >> "$OUTPUT_FASTA"
    else
        samtools faidx "$fasta" "$contig:$extract_start-$extract_end" | seqtk seq -r - >> "$OUTPUT_FASTA"
    fi

    # Rename sequence header to sample name
    sed -i "s/>$contig/>$sample/" "$OUTPUT_FASTA"

    echo "Extracted: $sample ($extract_start-$extract_end) from $contig"

done < "$RETRIEVE_PAF"
# Manual extraction for p26.B7 at 74:43178-47770
P26_B7_FASTA="${ASSEMBLY_DIR}/p26.B7.fasta.gz"
P26_B7_CONTIG="74"
P26_B7_START=43179
P26_B7_END=47951

if [[ -f "$P26_B7_FASTA" ]]; then
    echo "Extracting p26.B7 from $P26_B7_CONTIG:$P26_B7_START-$P26_B7_END"
    
    gunzip -c "$P26_B7_FASTA" > "${ASSEMBLY_DIR}/p26.B7.fasta"
    samtools faidx ${ASSEMBLY_DIR}/p26.B7.fasta "$P26_B7_CONTIG:$P26_B7_START-$P26_B7_END" >> "$OUTPUT_FASTA"

    # Rename sequence header
    sed -i "s/>$P26_B7_CONTIG/>p26.B7/" "$OUTPUT_FASTA"
else
    echo "Warning: FASTA file not found for p26.B7"
fi

echo "Extraction completed. Sequences saved in: $OUTPUT_FASTA"



#then add 3bp to the end of the two samples p20.D4 and p6.A10




#cat /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/scripts/step2_Oantigengenes/checkespE2_contigs/espE2_p25c2.fasta $OUTPUT_FASTA > $OUTPUT_FASTA2
