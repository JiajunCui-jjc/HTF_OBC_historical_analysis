#!/bin/bash

# Load environment
source /home/jiajucui/miniconda3/bin/activate phylogeny_snp

# Define directories

RESULT_DIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/checkespE2_57m/intermediate_seqs/OantigenP_but_espE2notpresent_p22D1_p25C11_p12F2_p24B5_p5C1usetagH1"
mkdir -p $RESULT_DIR
SAMPLE_LIST="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/scripts/step2_Oantigengenes/checkespE2_contigs/session2_retrieve_espE2_minimap/step3_OantigenP_but_espE2notpresent_p22D1_p25C11_p12F2_p24B5_p5C1usetagH1/fivesamples.txt"
ASSEMBLY_DIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fasta_m57"

MINIMAP_OUTPUT_DIR="${RESULT_DIR}/minimap_results"
QUERY="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/scripts/step2_Oantigengenes/checkespE2_contigs/session2_retrieve_espE2_minimap/step3_OantigenP_but_espE2notpresent_p22D1_p25C11_p12F2_p24B5_p5C1usetagH1/tagH1.fasta"
TOOLS="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/shfortailocin/tools"

# Create necessary directories
rm -r "$MINIMAP_OUTPUT_DIR"
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

#then check the sumpaf to see the coordinate:

#for example before we have:
#tagH_1|utg000001l_p25.C2:1015755-1017122
#>epsE_2|utg000001l_p25.C2:1017115-1021941
#so the epsE_2 start should be start_epsE2=(end_tagH1-7)

#p12.F2  tagH_1|utg000001l_p25.C2:1015755-1017122        1368    0       1097    +       1       1043802 1004158 1005255 1077    1097    60      NM:i:20 ms:i:697        AS:i:697        nn:i:0  tp:A:P  cm:i:73 s1:i:843        s2:i:0  de:f:0.0182     rl:i:0  cg:Z:1097M
#p22.D1  tagH_1|utg000001l_p25.C2:1015755-1017122        1368    0       1097    +       275     69660   41503   42600   1077    1097    60      NM:i:20 ms:i:697        AS:i:697        nn:i:0  tp:A:P  cm:i:73 s1:i:843        s2:i:0  de:f:0.0182     rl:i:0  cg:Z:1097M
#p25.C11 tagH_1|utg000001l_p25.C2:1015755-1017122        1368    0       1109    -       932     13858   2709    3818    1088    1109    60      NM:i:21 ms:i:689        AS:i:689        nn:i:0  tp:A:P  cm:i:68 s1:i:805        s2:i:0  de:f:0.0189     rl:i:0  cg:Z:1109M
#p5.C1   tagH_1|utg000001l_p25.C2:1015755-1017122        1368    0       1097    -       66      1004901 38767   39864   1077    1097    60      NM:i:20 ms:i:697        AS:i:697        nn:i:0  tp:A:P  cm:i:73 s1:i:843        s2:i:0  de:f:0.0182     rl:i:0  cg:Z:1097M
#no p24.B5 tagH1
#so the epsE_2 should be around: 
#p12.F2: strands: +
#				start: (1005255+(1368-1097)-7)=1005519
#				end: 1005519+4827=1010346
#				the coordinate: 1:1005519-1010346
#p22.D1: strands: +
#				start: (42600+(1368-1097)-7)=42864
#				end: 42864+4827=47691
#				the coordinate: 275:42864-47691
#p25.C11: strands: -
#			   start: 2709-(1368-1109)+7+1=2458
#				 end: 2458-4827... not complete but have a look though, 
#				 the coordinate: 932:0-2458
#p5.C1: strands: -
#			   start: 38767-(1368-1097)+7+1=38504
#				 end: 38504-4827=33677
#				 the coordinate: 66:33677-38504

gunzip -c ${ASSEMBLY_DIR}/p25.C11.fasta.gz > ${ASSEMBLY_DIR}/p25.C11.fasta
gunzip -c ${ASSEMBLY_DIR}/p22.D1.fasta.gz > ${ASSEMBLY_DIR}/p22.D1.fasta
gunzip -c ${ASSEMBLY_DIR}/p12.F2.fasta.gz > ${ASSEMBLY_DIR}/p12.F2.fasta
gunzip -c ${ASSEMBLY_DIR}/p5.C1.fasta.gz > ${ASSEMBLY_DIR}/p5.C1.fasta

samtools faidx ${ASSEMBLY_DIR}/p12.F2.fasta  1:1005519-1010346  > $RESULT_DIR/p12.F2.espE2.fasta

samtools faidx ${ASSEMBLY_DIR}/p22.D1.fasta  275:42864-47691  > $RESULT_DIR/p22.D1.espE2.fasta

samtools faidx ${ASSEMBLY_DIR}/p25.C11.fasta  932:0-2458 | seqtk seq -r - > $RESULT_DIR/p25.C11.espE2half.fasta

samtools faidx ${ASSEMBLY_DIR}/p5.C1.fasta  66:33677-38504 | seqtk seq -r - > $RESULT_DIR/p5.C1.espE2.fasta

