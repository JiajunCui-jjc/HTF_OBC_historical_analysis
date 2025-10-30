#!/bin/bash -l
#$ -l tmem=15G
#$ -l h_vmem=15G
#$ -l h_rt=2:00:0
#$ -cwd
#$ -V
#$ -N lopez10_espE2_assembly
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/
#$ -t 1-10

echo "Task id is $SGE_TASK_ID"

# ==========================================================
# 1. Environment and directories
# ==========================================================
source ~/miniconda3/bin/activate phylogeny_snp

bam_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/2025_h49/h10_lopez"
outbase="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/fish_historical/using39m_espE2fasta"
espE2_ref="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/checkespE2_57m/final_afterexpasy/completeopenframe.39fasta.txt"
threads=4

mkdir -p "${outbase}"/{fastq_from_bam/lopez10,mapped_bam_espE2,fastq_espE2only,tmp,assemblies,logs}

# ==========================================================
# 2. Filter BAMs: exclude HB* and PL*
# ==========================================================
bam_list=($(ls "$bam_dir"/*_Ps.sorted.bam | grep -vE '/HB|/PL' | sort -u | sed 's/\.\(collapsed\|r1r2\)_Ps\.sorted\.bam//' | uniq))
sample_base="${bam_list[$((SGE_TASK_ID-1))]}"
sample=$(basename "$sample_base")
echo "🔍 Processing sample base: $sample"

# ==========================================================
# 3. Collect the two BAMs (collapsed and r1r2 if exist)
# ==========================================================
bam_collapsed="${sample_base}.collapsed_Ps.sorted.bam"
bam_paired="${sample_base}.r1r2_Ps.sorted.bam"

# ----------------------------------------------------------
# Step A: Extract non-At reads (all reads in Ps BAMs)
# ----------------------------------------------------------
# collapsed
if [[ -f "$bam_collapsed" ]]; then
    fq_collapsed="${outbase}/fastq_from_bam/lopez10/${sample}_collapsed_nonAt.fq.gz"
    echo "[Step 1A] Extracting collapsed (single-end) non-At reads..."
    samtools fastq -@ $threads -n "$bam_collapsed" | gzip -c > "$fq_collapsed"
else
    fq_collapsed=""
fi

# paired
if [[ -f "$bam_paired" ]]; then
    fq_r1="${outbase}/fastq_from_bam/lopez10/${sample}_r1r2_nonAt.R1.fq.gz"
    fq_r2="${outbase}/fastq_from_bam/lopez10/${sample}_r1r2_nonAt.R2.fq.gz"
    echo "[Step 1B] Extracting paired non-At reads..."
    samtools fastq -@ $threads -n \
        -1 "$fq_r1" -2 "$fq_r2" \
        -0 /dev/null -s /dev/null \
        "$bam_paired"
else
    fq_r1=""; fq_r2=""
fi

# ==========================================================
# 4. Map to espE2 reference
# ==========================================================
# collapsed
if [[ -n "$fq_collapsed" ]]; then
    sai_c="${outbase}/tmp/${sample}_collapsed.sai"
    bam_c="${outbase}/mapped_bam_espE2/${sample}_collapsed.espE2.bam"
    bwa aln -t $threads "$espE2_ref" "$fq_collapsed" > "$sai_c"
    bwa samse "$espE2_ref" "$sai_c" "$fq_collapsed" | samtools sort -@ $threads -o "$bam_c"
    samtools index "$bam_c"
    rm -f "$sai_c"
fi

# paired
if [[ -n "$fq_r1" && -n "$fq_r2" ]]; then
    sai_r1="${outbase}/tmp/${sample}_R1.sai"
    sai_r2="${outbase}/tmp/${sample}_R2.sai"
    bam_p="${outbase}/mapped_bam_espE2/${sample}_r1r2.espE2.bam"
    bwa aln -t $threads -l 1024 "$espE2_ref" "$fq_r1" > "$sai_r1"
    bwa aln -t $threads -l 1024 "$espE2_ref" "$fq_r2" > "$sai_r2"
    bwa sampe "$espE2_ref" "$sai_r1" "$sai_r2" "$fq_r1" "$fq_r2" \
        | samtools sort -@ $threads -o "$bam_p"
    samtools index "$bam_p"
    rm -f "$sai_r1" "$sai_r2"
fi

# ==========================================================
# 5. Extract espE2-mapped reads from both & merge
# ==========================================================
mkdir -p ${outbase}/fastq_espE2only/lopez10/
fq_c_mapped="${outbase}/fastq_espE2only/lopez10/${sample}_collapsed_espE2mapped.fastq"
fq_p_mapped="${outbase}/fastq_espE2only/lopez10/${sample}_r1r2_espE2mapped.fastq"
fq_merged="$outbase/fastq_espE2only/${sample}_espE2mapped.fastq"

if [[ -f "${outbase}/mapped_bam_espE2/${sample}_collapsed.espE2.bam" ]]; then
    samtools fastq -@ $threads -F 4 "${outbase}/mapped_bam_espE2/${sample}_collapsed.espE2.bam" > "$fq_c_mapped"
fi

if [[ -f "${outbase}/mapped_bam_espE2/${sample}_r1r2.espE2.bam" ]]; then
    samtools fastq -@ $threads -F 4 "${outbase}/mapped_bam_espE2/${sample}_r1r2.espE2.bam" > "$fq_p_mapped"
fi

# merge (if both exist)
cat "$fq_c_mapped" "$fq_p_mapped" 2>/dev/null > "$fq_merged"
echo "✅ Merged FASTQ: $fq_merged"

