#!/bin/bash -l
#$ -l tmem=20G
#$ -l h_vmem=20G
#$ -l h_rt=10:00:0
#$ -N h40keykmerstep3
#$ -V
#$ -wd /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers
#$ -e /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers
#$ -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers

source ~/miniconda3/bin/activate phylogeny_snp
# --- Configuration ---fff
#first run rmdup bam and fastq
original="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h40/all_bams_h40/step1_bams_mappedtoPv_withHTFTFA"
bamdir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h40/all_bams_h40/step2_rmdup"
outdir="$bamdir/Pv_only_fastq_rmdup"

mkdir -p "$bamdir" "$outdir"

echo "=== Step 1: Mark duplicates and generate flagstats ==="
for bamfile in "$original"/*.bam; do
    samplename=$(basename "$bamfile" .mapped_to_Pseudomonas.dd.q20.bam)
    echo "Processing $samplename..."
    
    # Mark duplicates
    samtools markdup -@ 4 \
        "$original/${samplename}.mapped_to_Pseudomonas.dd.q20.bam" \
        "$bamdir/${samplename}.mapped_to_Pseudomonas.dd.q20.markeddup.bam"
    
    # Generate flagstat log
    samtools flagstat \
        "$bamdir/${samplename}.mapped_to_Pseudomonas.dd.q20.markeddup.bam" \
        > "$bamdir/${samplename}_Ps_q20_flagstats.log"
done

echo "=== Step 2: Convert markeddup BAM to FASTQ ==="
for markeddup in "$bamdir"/*.markeddup.bam; do
    sample=$(basename "$markeddup" .mapped_to_Pseudomonas.dd.q20.markeddup.bam)
    echo "Converting $sample to FASTQ..."
    
    echo "🔹 $sample: filtering and generating FASTQ..."
    samtools view -h -F 1024 "$bam" | samtools sort -n -@ 4 -o "$outdir/${sample}_filtered.bam"
    samtools fastq -n "$outdir/${sample}_filtered.bam" | gzip > "$outdir/${sample}.rmdup.inclHTF.fastq.gz"
    rm "$outdir/${sample}_filtered.bam"
done


#second, run calling jf
wd=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers
cd $wd

KMER=31
JFSIZE=300000000

h40_Pv_fastq="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h40/all_bams_h40/step2_rmdup/Pv_only_fastq_rmdup"
m57_raw_fastq=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_m57
#this is the fastq that coverted from bam files mapped to OTU5 ref with HTF and TFA... for h40
#for m57 we use the pv raw fastq

mkdir -p  step2_isolate_kmers step2_isolate_kmers/isolate_jf step2_isolate_kmers/isolate_jf/dump step3_query_kmers 
dump_out=$wd/step2_isolate_kmers/isolate_jf/dump



#h40
# Loop through all pair1 files
for fq in $h40_Pv_fastq/*.rmdup.inclHTF.fastq.gz; do
    sample=$(basename "$fq" .inclHTF.fastq.gz)
    iso_name=$sample
    echo "Processing $sample..."

    jellyfish count -m $KMER -s $JFSIZE -C \
        -o $wd/isolate_jf/${sample}.jf \
        <(zcat "$fq")
    jellyfish dump -c "$wd/step2_isolate_kmers/isolate_jf/${iso_name}.jf" > "$dump_out/${iso_name}_dump.txt"
done

