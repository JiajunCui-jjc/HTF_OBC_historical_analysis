#!/bin/bash -l
#$ -l tmem=6G
#$ -l h_vmem=6G
#$ -l h_rt=12:00:0
#$ -cwd
#$ -V
#$ -N m55_callvcf
#$ -t 1-55
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/
#$ -hold_jid m57maptoOTU5,m30maptoOTU5,OTU5refs2025_h39

# ============================================================
# SNP calling for 55 modern Pseudomonas viridiflava genomes
# Reference includes tailocin haplotypes (HTF + TFA)
# ============================================================

source ~/miniconda3/bin/activate phylogeny_snp
REF="/SAN/ugi/plant_genom/jiajucui/1_initial_data/reference_genome_Ps/Pseudomonas.OTU5_ref.fasta"
#REF="/SAN/ugi/plant_genom/jiajucui/1_initial_data/reference_genome_Ps_with_tailocin_haplotypes/Pseudomonas.plate25.C2.pilon.contigs_renamed.with_Tail_Fiber_Haps.fasta"
LIST="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/modern55.txt"
BAM_DIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_m57/all_bams_m57"
OUT_TOP="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step3_phylogeny_tree"
VCF_DIR="${OUT_TOP}/vcfs_m55"

mkdir -p "${VCF_DIR}"

sample=$(sed -n ${SGE_TASK_ID}p ${LIST})
bam="${BAM_DIR}/${sample}_mapped_Ps_q20.sorted.bam"

if [[ ! -f "$bam" ]]; then
    echo "[ERROR] BAM not found: $bam"
    exit 1
fi

echo "[INFO] Processing ${sample}"

tmpvcf="${VCF_DIR}/${sample}.tmp.vcf.gz"
finalvcf="${VCF_DIR}/${sample}.ploidy1_filtered.sorted.vcf.gz"

if [[ -f "${finalvcf}" ]]; then
    echo "[SKIP] ${sample} already processed."
    exit 0
fi

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
