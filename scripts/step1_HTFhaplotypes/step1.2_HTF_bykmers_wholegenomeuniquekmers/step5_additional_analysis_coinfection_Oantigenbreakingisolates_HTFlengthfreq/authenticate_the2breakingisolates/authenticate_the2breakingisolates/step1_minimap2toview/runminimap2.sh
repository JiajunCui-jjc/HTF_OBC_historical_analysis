#!/bin/bash

REF="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/HTFreference/tailocin_region.fa"
READ_DIR="the2fastq"
OUT_DIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step4_additionalanalysis/step1_2breakingisolates/bam_output"
THREADS=4
mkdir -p $READ_DIR 
ln -s /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h40/all_bams_h40/step1_bams_mappedtoPv_withHTFTFA/Pv_only_fastq_haventrmdup/PL0240.inclHTF.fastq.gz $READ_DIR
ln -s /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h40/all_bams_h40/step1_bams_mappedtoPv_withHTFTFA/Pv_only_fastq_haventrmdup/64.GBR_1933b_S36.inclHTF.fastq.gz $READ_DIR
mkdir -p  /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step4_additionalanalysis/  /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step4_additionalanalysis/step1_2breakingisolates "$OUT_DIR"

# Index reference
if [ ! -f "${REF}.mmi" ]; then
    echo "📦 Indexing $REF..."
    minimap2 -d "${REF}.mmi" "$REF"
fi

for fq in "$READ_DIR"/*.fastq.gz; do
    sample=$(basename "$fq" .fastq.gz)
    echo "🧬 Mapping $sample to tailocin region..."

    # BAM output
    minimap2 -t $THREADS -ax sr "${REF}.mmi" "$fq" | \
        samtools view -@ $THREADS -b - | \
        samtools sort -@ $THREADS -o "$OUT_DIR/${sample}.bam" -

    samtools index "$OUT_DIR/${sample}.bam"
    samtools flagstat "$OUT_DIR/${sample}.bam" > "$OUT_DIR/${sample}.flagstat.txt"

    # PAF output (summary mapping in tab format)
    minimap2 -t $THREADS -c "${REF}.mmi" "$fq" > "$OUT_DIR/${sample}.paf"
done

echo "✅ Finished. BAMs, flagstats, and PAFs are in $OUT_DIR/"
