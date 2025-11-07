#!/bin/bash -l
#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=24:00:0
#$ -cwd
#$ -V
#$ -N m53_merge_bivcf
#$ -hold_jid m53_callvcf
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/

# ============================================================
# Merge m53 (modern ATUE5)
# Filter SNPs and build tree
# ============================================================

source ~/miniconda3/bin/activate r2t

# ---- Working directories ----
OUT_TOP="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step3_phylogeny_tree"
RESULTS_DIR="${OUT_TOP}/results_m53_LDdecay"
mkdir -p "${RESULTS_DIR}"

# ---- Input VCF directories ----
VCF_M53="${OUT_TOP}/vcfs_m53"


# ============================================================
# Step 1. Link all historical (except 6 non-ATUE5) + modern 53 VCFs
# ============================================================
echo "[INFO] Linking h43 + m53 VCFs..."
cd "${VCF_M53}"



echo "[INFO] Total linked VCFs:"
ls *.vcf.gz | wc -l
cd "${RESULTS_DIR}"

# ============================================================
# Step 2. Merge all selected VCFs
# ============================================================
echo "[INFO] Merging m53 VCFs..."
/SAN/ugi/plant_genom/software/bcftools-1.11/bcftools merge \
    ${VCF_M53}/*.ploidy1_filtered.sorted.vcf.gz \
    -Oz -o merged_m53samples.vcf.gz

/SAN/ugi/plant_genom/software/bcftools-1.11/bcftools index -t merged_m53samples.vcf.gz

# ============================================================
# Step 3. SNP filtering and counts
# ============================================================
echo "[INFO] Filtering and counting SNPs..."

nofall=$(bcftools view -g ^miss -H merged_m53samples.vcf.gz | wc -l)

bcftools view -M2 merged_m53samples.vcf.gz | bgzip > nomulti_m53samples.vcf.gz
bcftools index -t nomulti_m53samples.vcf.gz
nofnomulti=$(bcftools view -g ^miss -H nomulti_m53samples.vcf.gz | wc -l)

bcftools view -m2 -M2 merged_m53samples.vcf.gz | gzip > bialleliconly_m53samples.vcf.gz
nofbiallelic=$(bcftools view -g ^miss -H bialleliconly_m53samples.vcf.gz | wc -l)

echo "Keeping full-info (no missing genotypes)..."
VCF_FULL="${RESULTS_DIR}/biallelic_fullinfo_m53samples.vcf.gz"

bcftools view -g ^miss bialleliconly_m53samples.vcf.gz -Oz -o "${VCF_FULL}"
bcftools index -f "${VCF_FULL}"


#vcftools LD need pure GT col, but not GT:PL by ploidy 1
# ============================================================
# Step 4. Fix genotype (GT) format for LD calculation
# ============================================================


echo "[INFO] Fixing GT field in ${VCF_FULL} for vcftools LD..."

# Output file path
VCF_LDREADY="${RESULTS_DIR}/biallelic_fullinfo_m53samples_GTfixed.vcf.gz"

# --- 1️⃣  Keep only GT field (remove PL, DP, etc.) ---
/SAN/ugi/plant_genom/software/bcftools-1.11/bcftools annotate \
    -x INFO,^FORMAT/GT "${VCF_FULL}" -Oz -o "${RESULTS_DIR}/tmp_GTonly.vcf.gz"
/SAN/ugi/plant_genom/software/bcftools-1.11/bcftools index -f "${RESULTS_DIR}/tmp_GTonly.vcf.gz"

# --- 2️⃣  Normalize malformed or missing GT entries ---
# (convert “0:0,.” → “0”, fill missing with 0)
export BCFTOOLS_PLUGINS=/SAN/ugi/plant_genom/software/bcftools-1.11/plugins

/SAN/ugi/plant_genom/software/bcftools-1.11/bcftools +setGT \
    "${RESULTS_DIR}/tmp_GTonly.vcf.gz" \
    -- -t q -n 0 -i 'FMT/GT="."' | \
/SAN/ugi/plant_genom/software/bcftools-1.11/bcftools view \
    -m2 -M2 -v snps -Oz -o "${VCF_LDREADY}"
/SAN/ugi/plant_genom/software/bcftools-1.11/bcftools index -f "${VCF_LDREADY}"

# --- 3️⃣  Clean temporary files ---
rm -f "${RESULTS_DIR}/tmp_GTonly.vcf.gz"*

# --- 4️⃣  Check the result ---
echo "[CHECK] Example genotypes (should be plain 0 or 1):"
bcftools view -H "${VCF_LDREADY}" | head -3 | cut -f 1-15

echo "[DONE] LD-ready VCF saved to:"
echo "       ${VCF_LDREADY}"




echo "allmergedSNPs,nomultiSNPs,biallelicSNPs" > SNPstats.txt
echo "${nofall},${nofnomulti},${nofbiallelic}" >> SNPstats.txt

echo "[INFO] SNP stats written to ${RESULTS_DIR}/SNPstats.txt"



