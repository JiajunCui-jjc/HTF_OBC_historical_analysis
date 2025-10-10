#!/bin/bash -l
#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=12:00:0
#$ -wd /SAN/ugi/plant_genom/jiajucui/
#$ -V
#$ -N markdup_AtPs_h46
#$ -t 1-46
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/

# ============================================================
# Mark duplicates (no -r) + coverage/depth for At & Ps
# ============================================================

source ~/miniconda3/bin/activate phylogeny_snp

i=$SGE_TASK_ID
list="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/scripts/step0_datapreprocessing_tailocin_presence/step1_extract_h40genomes_stats/step3_dedupstats_At_Ps/h46.txt"
sample=$(sed -n ${i}p ${list})     # clean IDs already match metadata

echo "===== [$(date)] Sample: ${sample} ====="

# --- Reference lengths ---
lenAt=119667750
lenPs=5941411

# --- Directories ---
indir_At="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/2024_413"
indir_Ps="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/tailocin_46/markdup"
outdir_tmp="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/2025_h46_markdup_stats"
outdir_merge="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/scripts/step0_datapreprocessing_tailocin_presence/step1_extract_h40genomes_stats/step3_dedupstats_At_Ps"
mkdir -p ${outdir_tmp} ${outdir_merge}

# --- Output files ---
At_table="${outdir_merge}/At_covdepth_h46.txt"
Ps_table="${outdir_merge}/Ps_covdepth_h46.txt"
merged_file="${outdir_merge}/allh46_markdup_covdepth_summary.txt"

# --- Input and output BAMs ---

bam_Ps="${indir_Ps}/${sample}.mapped_to_Pseudomonas.dd.q20.markeddup.bam"
# --- Identify correct BAM (paired vs single) ---
bam_At1="${indir_At}/${sample}_mapped_At_q20.sorted.bam"       # paired-end
bam_At2="${indir_At}/${sample}_mapped_At.q20.sorted.bam"       # single-end (Lopez)

if [[ -f "${bam_At1}" ]]; then
    bam_At="${bam_At1}"
    mode="paired"
elif [[ -f "${bam_At2}" ]]; then
    bam_At="${bam_At2}"
    mode="single"
else
    echo "[${sample}] ❌ Missing both BAM variants, skipping..."
    bam_At=""
    mode="none"
fi

bam_At_marked="${outdir_tmp}/${sample}_mapped_At_q20.sorted.markeddup.bam"

# --- Deduplication ---
if [[ "${mode}" == "paired" ]]; then
    echo "[${sample}] Running samtools markdup (paired-end)..."
    samtools markdup -r -@ 4 "${bam_At}" "${bam_At_marked}"
elif [[ "${mode}" == "single" ]]; then
    echo "[${sample}] Running samtools rmdup (single-end)..."
    samtools rmdup "${bam_At}" "${bam_At_marked}" 2> "${outdir_tmp}/${sample}_At_rmdup.log"
else
    echo "[${sample}] ❌ Skipped: no valid BAM found."
fi
# ============================================================
# Step 2. function to compute coverage/depth
# ============================================================
calc_cov_depth() {
    bamfile=$1
    reflen=$2
    exclude=$3
    if [[ -f "${bamfile}" ]]; then
        samtools coverage "${bamfile}" 2>/dev/null | \
        awk -v ex="${exclude}" '!/^#/ && $1 !~ ex' | \
        awk -v len=$reflen '
            {sumcov+=$5; sumlen+=$3; sumd+=$3*$7}
            END{
                if(sumlen>0){cov=sumcov/sumlen; dep=sumd/len; print cov"\t"dep}
                else{print "NA\tNA"}
            }'
    else
        echo -e "NA\tNA"
    fi
}

# ============================================================
# Step 3. Compute coverage/depth (exclude tailocin contigs for Ps)
# ============================================================
read At_cov At_depth <<<$(calc_cov_depth "${bam_At_marked}" ${lenAt} "^$")
read Ps_cov Ps_depth <<<$(calc_cov_depth "${bam_Ps}" ${lenPs} "^(TFA_|HTF_)")

# ============================================================
# Step 4. Append to species-specific tables + master summary
# ============================================================
{
    flock -x 200
    # --- A. thaliana ---
    if ! grep -q "samplename" "${At_table}" 2>/dev/null; then
        echo -e "samplename\tAt_cov_markdup\tAt_depth_markdup" > "${At_table}"
    fi
    echo -e "${sample}\t${At_cov}\t${At_depth}" >> "${At_table}"

    # --- Pseudomonas ---
    if ! grep -q "samplename" "${Ps_table}" 2>/dev/null; then
        echo -e "samplename\tPs_cov_markdup\tPs_depth_markdup" > "${Ps_table}"
    fi
    echo -e "${sample}\t${Ps_cov}\t${Ps_depth}" >> "${Ps_table}"

    # --- Combined summary ---
    if ! grep -q "samplename" "${merged_file}" 2>/dev/null; then
        echo -e "samplename\tAt_cov_markdup\tAt_depth_markdup\tPs_cov_markdup\tPs_depth_markdup" > "${merged_file}"
    fi
    echo -e "${sample}\t${At_cov}\t${At_depth}\t${Ps_cov}\t${Ps_depth}" >> "${merged_file}"
} 200>"${merged_file}.lock"

echo "[${sample}] ✅ Done (At_cov=${At_cov}, Ps_cov=${Ps_cov})"
