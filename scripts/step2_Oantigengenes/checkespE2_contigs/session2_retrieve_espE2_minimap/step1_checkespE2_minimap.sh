#!/bin/bash

# Load environment
source /home/jiajucui/miniconda3/bin/activate phylogeny_snp

# Define directories
RESULT_DIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/checkespE2_57m"
SAMPLE_LIST="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fasta_m57/modern57.txt"
ASSEMBLY_DIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fasta_m57"
MINIMAP_OUTPUT_DIR="${RESULT_DIR}/minimap_results"
QUERY="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/scripts/step2_Oantigengenes/checkespE2_contigs/espE2_p25c2.fasta"
TOOLS="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/shfortailocin/tools"

# Create necessary directories
#rm -r "$MINIMAP_OUTPUT_DIR"
mkdir -p "$MINIMAP_OUTPUT_DIR"
mkdir -p "$RESULT_DIR/tmp"

# Load Minimap2 path
MINIMAP_CMD="$TOOLS/minimap2/minimap2"

# Process each sample
while read -r sample; do
    echo "Processing Sample: $sample"

    # Construct file paths
    fasta="${ASSEMBLY_DIR}/${sample}.fasta.gz"
    paf_file="${MINIMAP_OUTPUT_DIR}/${sample}_alignment.paf"

    # Check if the FASTA file exists before running Minimap2
    if [[ -f "$fasta" ]]; then
        $MINIMAP_CMD -cx asm5 "$fasta" "$QUERY" > "$paf_file"
    else
        echo "Warning: File not found - $fasta"
    fi

done < "$SAMPLE_LIST"




cd $MINIMAP_OUTPUT_DIR
awk '{print FILENAME "\t" $0}' *.paf | sed 's/_alignment.paf//g' | less -S > sumpaf.txt

#subgroup the samples...
# Input file (you can change this to your actual filename)
input="sumpaf.txt"
ls -lh | grep ' 0 ' | cut -d ' ' -f 12 | cut -d '_' -f 1 > sub_20_absent_nomatch.txt
# Output files
full="sub_7_fulllengthmatched.txt"
partial="sub_30_retrievepaf.txt"

# Subset based on column 11 (matched length)
awk '$11 == 4827' "$input" > "$full"
awk '$11 != 4827' "$input" > "$partial"
