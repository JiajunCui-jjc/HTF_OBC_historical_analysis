#!/bin/bash
#$ -l tmem=50G
#$ -l h_vmem=50G
#$ -l h_rt=2:00:0
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N spadeh49
#$ -t 1-49
#$ -hold_jid lopez10_espE2_assembly
#$ -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/fish_historical/logs
#$ -e /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/fish_historical/logs


# ===============================================================
#  SPAdes + minimap2 pipeline (array self-submitting version)
# ===============================================================

# --- Activate env ---
source /home/jiajucui/miniconda3/bin/activate phylogeny_snp

# --- Base dirs ---
output_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/fish_historical/using39m_espE2fasta"
reference_genome="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/checkespE2_57m/final_afterexpasy/completeopenframe.39fasta.txt"
trimmed_fastq_dir="${output_dir}/fastq_espE2only"
assembly_dir="${output_dir}/assemblies"
mapping_dir="${output_dir}/mappings"
tmp_dir="${output_dir}/tmp"
tools="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/shfortailocin/tools"

mkdir -p "$assembly_dir" "$mapping_dir" "$tmp_dir"

SAMPLE_LIST="${output_dir}/sample_list.txt"

# ===============================================================
#  Stage 1 — if not in array, build list and resubmit
# ===============================================================
#    echo "🧩 Generating sample list from ${trimmed_fastq_dir}"
#    find "$trimmed_fastq_dir" -type f -name "*_espE2mapped.fastq" \
#        | xargs -n1 basename \
#        | sed 's/_espE2mapped.fastq//' \
#        | sort > "$SAMPLE_LIST"
#
#    N=$(wc -l < "$SAMPLE_LIST")
#    echo "✅ Found $N samples. Submitting array job..."


# ===============================================================
#  Stage 2 — inside the array: run SPAdes + minimap2
# ===============================================================
SAMPLE=$(sed -n "${SGE_TASK_ID}p" "$SAMPLE_LIST")
if [[ -z "$SAMPLE" ]]; then
    echo "❌ Invalid SGE_TASK_ID=${SGE_TASK_ID}."
    exit 1
fi

fq="${trimmed_fastq_dir}/${SAMPLE}_espE2mapped.fastq"
assembly_out="${assembly_dir}/${SAMPLE}"
contig_file="${assembly_out}/contigs.fasta"
paf_file="${mapping_dir}/${SAMPLE}_mapped.paf"
nocontigs_file="${mapping_dir}/nocontigs.txt"
mkdir -p "$assembly_out"
touch "$nocontigs_file"

echo "=== [Task ${SGE_TASK_ID}] Processing: ${SAMPLE} ==="
echo "FASTQ: ${fq}"

# ---- SPAdes ----
if [[ ! -f "$contig_file" ]]; then
    echo "[SPAdes] assembling..."
    spades.py -s "$fq" -o "$assembly_out" -k 17,21
else
    echo "[SPAdes] contigs already exist – skipping."
fi

# ---- Check assembly ----
if [[ ! -f "$contig_file" ]]; then
    echo "$SAMPLE" >> "$nocontigs_file"
    echo "❌ No contigs for $SAMPLE"
    exit 0
fi

# ---- minimap2 ----
echo "[minimap2] mapping..."
"$tools/minimap2/minimap2" -cx asm5 "$reference_genome" "$contig_file" > "$paf_file"

if [[ ! -s "$paf_file" ]]; then
    echo "$SAMPLE nomappedcontig" >> "$nocontigs_file"
    echo "❌ No mapping for $SAMPLE"
    exit 0
fi

# ---- Count contigs ----
contig_count=$(grep -c "^>" "$contig_file" || true)
echo -e "${SAMPLE}\t${contig_count}" >> "${mapping_dir}/contig_counts.txt"

echo "✅ Finished $SAMPLE (contigs: $contig_count)"

