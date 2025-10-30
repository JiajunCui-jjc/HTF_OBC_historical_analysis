#!/bin/bash -l
#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=24:00:0
#$ -cwd
#$ -V
#$ -N h40_m53_merge_tree
#$ -hold_jid h46_callvcf,m53_callvcf
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/
#$ -hold_jid h46_callvcf,m30_callvcf,m53_callvcf,m55_callvcf

# ============================================================
# Merge h40 (ATUE5-only historical) + m53 (modern ATUE5)
# Filter SNPs and build tree
# ============================================================

source ~/miniconda3/bin/activate r2t

# ---- Working directories ----
OUT_TOP="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step3_phylogeny_tree"
RESULTS_DIR="${OUT_TOP}/results_h40_m53"
TREE_DIR="${RESULTS_DIR}/tree"
mkdir -p "${RESULTS_DIR}" "${TREE_DIR}"

# ---- Input VCF directories ----
VCF_H46="${OUT_TOP}/vcfs_h46"
VCF_M53="${OUT_TOP}/vcfs_m53"
VCF_COMBINED="${RESULTS_DIR}/vcfs_h40_m53"
mkdir -p "${VCF_COMBINED}"

# ---- Historical samples to exclude ----
EXCLUDE_SAMPLES=("HB0828" "HB0863" "PL0066" "PL0108" "PL0203" "PL0258")
#BED="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/scripts/step0_datapreprocessing_tailocin_presence/step3_phylogeny_reconstruction/step1_runh46_m53_30_tree/step3_merge_buildtree/exclude_HTF_TFA_auto.bed"

# ============================================================
# Step 1. Link all historical (except 6 non-ATUE5) + modern 53 VCFs
# ============================================================
echo "[INFO] Linking h40 + m53 VCFs..."
cd "${VCF_COMBINED}"

# Link historical 40 (filter out non-ATUE5)
for f in ${VCF_H46}/*.ploidy1_filtered.sorted.vcf.gz; do
    sample=$(basename "$f" .ploidy1_filtered.sorted.vcf.gz)
    skip=false
    for ex in "${EXCLUDE_SAMPLES[@]}"; do
        [[ "$sample" == "$ex" ]] && skip=true && break
    done
    if [ "$skip" = false ]; then
        ln -sf "$f" .
        ln -sf "$f.tbi" .
    fi
done

# Link 53 modern ATUE5
ln -sf ${VCF_M53}/*.ploidy1_filtered.sorted.vcf.gz* .

echo "[INFO] Total linked VCFs:"
ls *.vcf.gz | wc -l
cd "${RESULTS_DIR}"

# ============================================================
# Step 2. Merge all selected VCFs
# ============================================================
echo "[INFO] Merging h40 + m53 VCFs..."
/SAN/ugi/plant_genom/software/bcftools-1.11/bcftools merge \
    ${VCF_COMBINED}/*.ploidy1_filtered.sorted.vcf.gz \
    -Oz -o merged_h40_m53_93samples.vcf.gz

/SAN/ugi/plant_genom/software/bcftools-1.11/bcftools index -t merged_h40_m53_93samples.vcf.gz

# ============================================================
# Step 3. SNP filtering and counts
# ============================================================
echo "[INFO] Filtering and counting SNPs..."

nofall=$(bcftools view -g ^miss -H merged_h40_m53_93samples.vcf.gz | wc -l)

bcftools view -M2 merged_h40_m53_93samples.vcf.gz | bgzip > nomulti_h40_m53_93samples.vcf.gz
bcftools index -t nomulti_h40_m53_93samples.vcf.gz
nofnomulti=$(bcftools view -g ^miss -H nomulti_h40_m53_93samples.vcf.gz | wc -l)

bcftools view -m2 -M2 merged_h40_m53_93samples.vcf.gz | gzip > bialleliconly_h40_m53_93samples.vcf.gz
nofbiallelic=$(bcftools view -g ^miss -H bialleliconly_h40_m53_93samples.vcf.gz | wc -l)

echo "allmergedSNPs,nomultiSNPs,biallelicSNPs" > SNPstats.txt
echo "${nofall},${nofnomulti},${nofbiallelic}" >> SNPstats.txt

echo "[INFO] SNP stats written to ${RESULTS_DIR}/SNPstats.txt"
# ============================================================
# Step 4. Convert to PHYLIP and run IQ-TREE
# ============================================================
echo "[INFO] Converting to PHYLIP..."
python /SAN/ugi/plant_genom/jiajucui/phylogeny/phylogeny_snp/vcf2phylip/vcf2phylip.py \
    -i bialleliconly_h40_m53_93samples.vcf.gz -m 93

mv bialleliconly_h40_m53_93samples.min93.phy "${TREE_DIR}/"

cd "${TREE_DIR}"
echo "[INFO] Running IQ-TREE..."
iqtree -T AUTO \
    -s bialleliconly_h40_m53_93samples.min93.phy \
    -m MFP -B 1000 -alrt 1000 | tee iqtree_93.log

echo "[DONE] Tree built successfully. Results in:"
echo " - ${RESULTS_DIR}/"
echo " - ${TREE_DIR}/"
