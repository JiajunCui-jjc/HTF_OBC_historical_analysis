#!/bin/bash -l
#$ -l tmem=8G
#$ -l h_vmem=8G
#$ -l h_rt=12:00:0
#$ -cwd
#$ -N m57maptoOTU5
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/
#$ -j y
#$ -V
#$ -t 1-57   # adjust to your total number of samples

# ==========================================================
# Map paired-end reads → dedup → sort → index (final clean version)
# ==========================================================


# ---- Paths ----
WD="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_m57"
OUTDIR="${WD}/all_bams_m57_OTU5refonly_dedup"
#REF=/SAN/ugi/plant_genom/jiajucui/1_initial_data/reference_genome_Ps_with_tailocin_haplotypes/Pseudomonas.plate25.C2.pilon.contigs_renamed.with_Tail_Fiber_Haps.fasta
REF="/SAN/ugi/plant_genom/jiajucui/1_initial_data/reference_genome_Ps/Pseudomonas.OTU5_ref.fasta"

#use this ref to make sure when calling snp we are using the same ref as historical samples
mkdir -p "${OUTDIR}/tmp"

# ---- Determine sample ----
SAMPLE=$(ls ${WD}/*.pair1.truncated.gz | sed -n ${SGE_TASK_ID}p | xargs -n1 basename | sed 's/.pair1.truncated.gz//')
echo "[$(date)] Processing ${SAMPLE}"
#run adapterremoval with nsqs like in historical preprocessing, since its modern we use truncated (most are long reads)
# ---- Input files ----
R1="${WD}/${SAMPLE}.pair1.truncated.gz"
R2="${WD}/${SAMPLE}.pair2.truncated.gz"

# ---- Output files ----
BAM_FINAL="${OUTDIR}/${SAMPLE}_mapped_Ps_q20.dedup.sorted.bam"

# ---- 1. BWA alignment ----
bwa aln -t 2 -l 1024 -f "${OUTDIR}/tmp/${SAMPLE}.r1.sai" "${REF}" "${R1}"
bwa aln -t 2 -l 1024 -f "${OUTDIR}/tmp/${SAMPLE}.r2.sai" "${REF}" "${R2}"

bwa sampe -r "@RG\tID:${SAMPLE}\tSM:${SAMPLE}" \
  -f "${OUTDIR}/tmp/${SAMPLE}.sam" \
  "${REF}" \
  "${OUTDIR}/tmp/${SAMPLE}.r1.sai" "${OUTDIR}/tmp/${SAMPLE}.r2.sai" \
  "${R1}" "${R2}"

# ---- 2. Convert SAM → BAM (q20 filter) ----
samtools view -@ 2 -F 4 -q 20 -Sbh \
  -o "${OUTDIR}/tmp/${SAMPLE}.filtered.bam" \
  "${OUTDIR}/tmp/${SAMPLE}.sam"

samtools sort -@ 2 -o "${OUTDIR}/tmp/${SAMPLE}.sorted_pre_dedup.bam" "${OUTDIR}/tmp/${SAMPLE}.filtered.bam"
samtools rmdup "${OUTDIR}/tmp/${SAMPLE}.sorted_pre_dedup.bam" "${OUTDIR}/tmp/${SAMPLE}.dedup.bam"

# ---- 4. Final sort and index ----
samtools sort -@ 2 -o "${BAM_FINAL}" "${OUTDIR}/tmp/${SAMPLE}.dedup.bam"
samtools index "${BAM_FINAL}"

# ---- 5. Cleanup intermediates ----
rm -f "${OUTDIR}/tmp/${SAMPLE}.sam" \
      "${OUTDIR}/tmp/${SAMPLE}.r1.sai" \
      "${OUTDIR}/tmp/${SAMPLE}.r2.sai" \
      "${OUTDIR}/tmp/${SAMPLE}.filtered.bam" \

echo "[$(date)] ✅ Done: ${SAMPLE} → ${BAM_FINAL}"
