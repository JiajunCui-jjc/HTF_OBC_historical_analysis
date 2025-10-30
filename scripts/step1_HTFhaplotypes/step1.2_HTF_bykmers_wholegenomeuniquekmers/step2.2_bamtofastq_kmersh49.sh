#!/bin/bash -l
#$ -l tmem=20G
#$ -l h_vmem=20G
#$ -l h_rt=10:00:0
#$ -N h49keykmerstep3
#$ -V
#$ -wd /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers
#$ -e /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/logs
#$ -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/logs

source ~/miniconda3/bin/activate phylogeny_snp
# --- Configuration ---fff
#first run rmdup bam and fastq
bamdir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h49/Ps_bams_maptoOTU5_with_haplotype_h49_dedup/h49_bams_softlink/"
outdir="$bamdir/OTU5withHTFTFA_fastq_rmdup"

mkdir -p "$bamdir" "$outdir"

#.mapped_to_Pseudomonas.dd.q20.markeddup.bam 
echo "=== Step 2: Convert markeddup BAM to FASTQ ==="
for bam in "$bamdir"/*.markeddup.bam; do
    sample=$(basename "$bam" .mapped_to_Pseudomonas.dd.q20.markeddup.bam)
    echo "Converting $sample to FASTQ..."
    
    echo "🔹 $sample: filtering and generating FASTQ..."
    samtools view -h -F 1024 "$bam" | samtools sort -n -@ 4 -o "$outdir/${sample}_filtered.bam"
    samtools fastq -n "$outdir/${sample}_filtered.bam" | gzip > "$outdir/${sample}.rmdup.inclHTF.fastq.gz"
    rm "$outdir/${sample}_filtered.bam"
done


