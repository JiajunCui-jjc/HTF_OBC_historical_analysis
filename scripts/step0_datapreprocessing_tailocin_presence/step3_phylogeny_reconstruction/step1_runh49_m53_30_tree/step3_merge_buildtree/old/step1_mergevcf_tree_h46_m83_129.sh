#!/bin/bash -l
#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=24:00:0
#$ -cwd
#$ -V
#$ -N hm129_merge_tree
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/
#$ -hold_jid h46_callvcf,m30_callvcf,m53_callvcf,m55_callvcf

# ============================================================
# Merge 129 VCFs (h46 + m53 + m30), filter SNPs, and build tree
# ============================================================

source ~/miniconda3/bin/activate r2t

# ---- Working and output directories ----
OUT_TOP="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step3_phylogeny_tree"
RESULTS_DIR="${OUT_TOP}/results_h46_m83_129"
TREE_DIR="${RESULTS_DIR}/tree"

mkdir -p "${RESULTS_DIR}" "${TREE_DIR}"

# ---- Input VCF directories ----
VCF_H46="${OUT_TOP}/vcfs_h46"
VCF_M53="${OUT_TOP}/vcfs_m53"
VCF_M30="${OUT_TOP}/vcfs_m30"
VCF_COMBINED="${RESULTS_DIR}/vcfs_129_h46_m83"

mkdir -p "${VCF_COMBINED}"

# ============================================================
# Step 1. Link all VCFs into one combined folder
# ============================================================
echo "[INFO] Linking all 129 VCFs..."
cd "${VCF_COMBINED}"

ln -sf ${VCF_H46}/*.ploidy1_filtered.sorted.vcf.gz* .
ln -sf ${VCF_M53}/*.ploidy1_filtered.sorted.vcf.gz* .
ln -sf ${VCF_M30}/*.ploidy1_filtered.sorted.vcf.gz* .

echo "[INFO] Total VCFs linked:"
ls *.vcf.gz | wc -l
cd "${RESULTS_DIR}"

#BED="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/scripts/step0_datapreprocessing_tailocin_presence/step3_phylogeny_reconstruction/step1_runh46_m53_30_tree/step3_merge_buildtree/exclude_HTF_TFA_auto.bed"

# ============================================================
# Step 2. Merge all VCFs
# ============================================================
echo "[INFO] Merging all VCFs..."
/SAN/ugi/plant_genom/software/bcftools-1.11/bcftools merge \
    ${VCF_COMBINED}/*.ploidy1_filtered.sorted.vcf.gz \
    -Oz -o merged_h46_m83_129samples.vcf.gz

/SAN/ugi/plant_genom/software/bcftools-1.11/bcftools index -t merged_h46_m83_129samples.vcf.gz


# ============================================================
# Step 3. SNP filtering and counts
# ============================================================
echo "[INFO] Filtering and counting SNPs..."

nofall=$(bcftools view -g ^miss -H merged_h46_m83_129samples.vcf.gz | wc -l)

bcftools view -M2 merged_h46_m83_129samples.vcf.gz | bgzip > nomulti_h46_m83_129samples.vcf.gz
bcftools index -t nomulti_h46_m83_129samples.vcf.gz
nofnomulti=$(bcftools view -g ^miss -H nomulti_h46_m83_129samples.vcf.gz | wc -l)

bcftools view -m2 -M2 merged_h46_m83_129samples.vcf.gz | gzip > bialleliconly_h46_m83_129samples.vcf.gz
nofbiallelic=$(bcftools view -g ^miss -H bialleliconly_h46_m83_129samples.vcf.gz | wc -l)

echo "allmergedSNPs,nomultiSNPs,biallelicSNPs" > SNPstats.txt
echo "${nofall},${nofnomulti},${nofbiallelic}" >> SNPstats.txt
echo "[INFO] SNP stats written to ${RESULTS_DIR}/SNPstats.txt"

# ============================================================
# Step 4. Convert to PHYLIP and build ML tree
# ============================================================
echo "[INFO] Converting to PHYLIP..."
python /SAN/ugi/plant_genom/jiajucui/phylogeny/phylogeny_snp/vcf2phylip/vcf2phylip.py \
    -i bialleliconly_h46_m83_129samples.vcf.gz -m 129

mv bialleliconly_h46_m83_129samples.min129.phy "${TREE_DIR}/"

cd "${TREE_DIR}"
echo "[INFO] Running IQ-TREE..."
iqtree -T AUTO \
    -s bialleliconly_h46_m83_129samples.min129.phy \
    -m MFP -B 1000 -alrt 1000 | tee iqtree_129.log

echo "[DONE] Tree building complete. Results in:"
echo " - ${RESULTS_DIR}/"
echo " - ${TREE_DIR}/"
