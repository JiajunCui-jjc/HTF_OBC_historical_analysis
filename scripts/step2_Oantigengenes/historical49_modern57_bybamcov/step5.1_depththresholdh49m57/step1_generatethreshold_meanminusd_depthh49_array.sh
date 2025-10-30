#!/bin/bash -l
#$ -l tmem=6G
#$ -l h_vmem=6G
#$ -l h_rt=1:00:0
#$ -cwd
#$ -V
#$ -j y
#$ -N h49_depthstats
#$ -t 1-49
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/

set -euo pipefail

# =============================================================
# Compute genome-wide depth statistics per sample
# Run as an SGE array: one sample per task (1–49)
# =============================================================

RESULTS_DIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/historical49_modern57_bybamcov/coverage_h49_m57"
BAM_DIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h49/all_raw_meta_fastq_h49/../Ps_bams_maptoOTU5_with_haplotype_h49_dedup/h49_bams_softlink"
SAMPLE_LIST="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/historical49.txt"
GENOME_LENGTH=5941411  # reference genome length

OUT_DIR="$RESULTS_DIR"
TMP_DIR="${OUT_DIR}/tmp"
mkdir -p "$TMP_DIR"

OUT_FILE="${OUT_DIR}/depth_thresholds_mean_minus_025sd.tsv"
TMP_FILE="${TMP_DIR}/tmp_${SGE_TASK_ID}.tsv"

# --- pick sample for this array task ---
sample=$(sed -n "${SGE_TASK_ID}p" "$SAMPLE_LIST")
[[ -z "$sample" ]] && { echo "No sample for task ID $SGE_TASK_ID"; exit 1; }

# --- find BAM ---
BAM_FILE=$(find "$BAM_DIR" -maxdepth 1 -name "${sample}*.bam" | head -n1)
if [[ -z "$BAM_FILE" || ! -f "$BAM_FILE" ]]; then
    echo "No BAM for $sample"
    exit 0
fi

echo "[$(date)] Processing $sample ($BAM_FILE)"

# --- compute mean, SD, threshold ---
samtools depth -aa "$BAM_FILE" | \
awk -v total="$GENOME_LENGTH" -v s="$sample" '
{
    d[NR]=$3
    sum+=$3
}
END {
    n=NR
    if (n>0) {
        mean=sum/total
        for (i=1;i<=n;i++){diff=d[i]-mean; var+=diff*diff}
        sd=sqrt(var/n)
        thr=mean-(0.25*sd)
        if(thr<0) thr=0
        rel_thr=(mean>0)? thr/mean : 0
        printf "%s\t%.6f\t%.6f\t%.6f\t%.6f\n", s, mean, sd, thr, rel_thr
    }
}' > "$TMP_FILE"

echo "[$(date)] Done $sample → $TMP_FILE"

# =============================================================
# After all tasks finish, merge the results manually:
# =============================================================
#   cd "/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/historical49_modern57_bybamcov/coverage_h49_m57"
#   echo -e "Sample\tMean_total_depth\tSD_total_depth\tMean_minus_SD\trelative_threshold" > depth_thresholds_mean_minus_sd.tsv
#   cat tmp/tmp_*.tsv >> depth_thresholds_mean_minus_sd.tsv
#   rm -r tmp
# =============================================================

