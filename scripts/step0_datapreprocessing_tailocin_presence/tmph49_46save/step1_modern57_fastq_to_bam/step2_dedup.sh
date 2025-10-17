#!/bin/bash -l
#$ -l tmem=8G              # memory per job
#$ -l h_vmem=8G
#$ -l h_rt=2:00:0          # runtime
#$ -cwd
#$ -j y
#$ -V
#$ -t 1-57

# ==========================================================
# SGE array job: deduplicate each BAM using samtools rmdup
# ==========================================================

set -euo pipefail


indir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_m57/all_bams_m57"
before="${indir}/beforededup"
#first mkdir beforededup
#mv *.bam *.bai beforededup/
mkdir -p ${before}/tmp/
# --- Determine which BAM to process ---
bam=$(ls ${before}/*.bam | sed -n ${SGE_TASK_ID}p)
sample=$(basename "$bam" .bam)
dedupbam="${indir}/${sample}.dedup.bam"

echo "[$(date)] Processing $sample ..."
samtools sort -o "${before}/tmp/${sample}.sorted.bam" "$bam"
samtools rmdup "${before}/tmp/${sample}.sorted.bam" "$dedupbam"
samtools index "$dedupbam"
rm "${before}/tmp/${sample}.sorted.bam"
echo "[$(date)] Done: $dedupbam"
#after all check
#for f in beforededup/*.bam; do
#    base=$(basename "$f" .bam)
#    if [[ ! -f ${base}.dedup.bam ]]; then
#        echo "Missing dedup for: $base"
#    fi
#done
