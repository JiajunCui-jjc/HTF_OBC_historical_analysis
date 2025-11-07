#!/bin/bash -l
#$ -l tmem=16G
#$ -l h_vmem=16G
#$ -l h_rt=12:00:0
#$ -cwd
#$ -V
#$ -N h49_callvcf
#$ -t 1-49
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/
#$ -hold_jid m57maptoOTU5,m30maptoOTU5,OTU5refs2025_h39

# ============================================================
# SNP calling for 3 historical Pseudomonas viridiflava genomes
# Reference includes tailocin haplotypes (HTF + TFA)
# ============================================================
export TMPDIR=/SAN/ugi/plant_genom/jiajucui/tmp_bcftools
mkdir -p "$TMPDIR"

source ~/miniconda3/bin/activate phylogeny_snp

# --- reference without HTF + TFA regions ---
REF=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/HTFreference/OTU5_ref_noHTFhaplotypes/ref.fasta


# --- input sample list ---
LIST="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/historical49.txt"

# --- BAM directory ---
BAM_DIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data//all_fastq_h49/Ps_bams_maptoOTU5withouthaplotype_h49_dedup/h49_bams_softlink/"

# --- output directory ---
OUT_TOP="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step3_phylogeny_tree"
VCF_DIR="${OUT_TOP}/vcfs_h49"

mkdir -p "${VCF_DIR}"

# ============================================================

sample=$(sed -n ${SGE_TASK_ID}p ${LIST})
bam="${BAM_DIR}/${sample}.mapped_to_Pseudomonas.dd.q20.markeddup.bam"
#samtools index $bam
if [[ ! -f "$bam" ]]; then
    echo "[ERROR] BAM not found: $bam"
    exit 1
fi

echo "[INFO] Processing historical sample ${sample}"

tmpvcf="${VCF_DIR}/${sample}.tmp.vcf.gz"
finalvcf="${VCF_DIR}/${sample}.ploidy1_filtered.sorted.vcf.gz"

if [[ -f "${finalvcf}" ]]; then
    echo "[SKIP] ${sample} already processed."
    exit 0
fi

# --- SNP calling ---
/SAN/ugi/plant_genom/software/bcftools-1.11/bcftools mpileup \
    -Ou -f "${REF}" "${bam}" | \
/SAN/ugi/plant_genom/software/bcftools-1.11/bcftools call \
    -c --skip-variants indels --ploidy 1 -Oz | \
/SAN/ugi/plant_genom/software/bcftools-1.11/bcftools view \
    -i '%QUAL>=20' -Oz -o "${tmpvcf}"

/SAN/ugi/plant_genom/software/bcftools-1.11/bcftools sort \
    "${tmpvcf}" -Oz -o "${finalvcf}"
/SAN/ugi/plant_genom/software/bcftools-1.11/bcftools index -t "${finalvcf}"

rm -f "${tmpvcf}"

echo "[DONE] ${sample}"
