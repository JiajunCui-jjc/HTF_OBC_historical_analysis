#!/bin/bash
source /home/jiajucui/miniconda3/bin/activate phylogeny_snp

# === Paths and settings ===
REF="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/HTFreference/tailocin_region.fa"
READ_DIR="the2fastq"
THREADS=4

# === Prepare fastq directory and link files ===
mkdir -p "$READ_DIR"

#ln -s /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h40/all_bams_h40/step1_bams_mappedtoPv_withHTFTFA/Pv_only_fastq_haventrmdup/PL0240.inclHTF.fastq.gz "$READ_DIR"
#ln -s /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h40/all_bams_h40/step1_bams_mappedtoPv_withHTFTFA/Pv_only_fastq_haventrmdup/64.GBR_1933b_S36.inclHTF.fastq.gz "$READ_DIR"

# === Index the reference if not done ===
if [ ! -f "${REF}.mmi" ]; then
    echo "📦 Indexing $REF..."
    minimap2 -d "${REF}.mmi" "$REF"
fi

# === Main loop: Assemble then map ===
for fq in "$READ_DIR"/*.fastq.gz; do
    sample=$(basename "$fq" .fastq.gz)
    outdir="./tmp_assemble_from_raw_Psfastq_thenminmap_${sample}"
    mkdir -p "$outdir"

    echo "🧬 Assembling $sample with SPAdes..."
    spades.py -s "$fq" -o "$outdir/spades_output" -t $THREADS --careful -k 21,33 > "$outdir/spades.log" 2>&1

    contigs="$outdir/spades_output/contigs.fasta"
    if [ ! -f "$contigs" ]; then
        echo "❌ Assembly failed or contigs not generated for $sample"
        continue
    fi

    echo "🔍 Mapping assembled contigs of $sample to tailocin region..."

    # BAM output
    minimap2 -t $THREADS -ax asm5 "${REF}.mmi" "$contigs" | \
        samtools view -@ $THREADS -b - | \
        samtools sort -@ $THREADS -o "$outdir/${sample}.bam" -

    samtools index "$outdir/${sample}.bam"
    samtools flagstat "$outdir/${sample}.bam" > "$outdir/${sample}.flagstat.txt"

    # PAF output
    minimap2 -t $THREADS -c "${REF}.mmi" "$contigs" > "$outdir/${sample}.paf"

    echo "✅ Done: $sample → $outdir"
done

echo "🎯 All assemblies and mappings completed."
