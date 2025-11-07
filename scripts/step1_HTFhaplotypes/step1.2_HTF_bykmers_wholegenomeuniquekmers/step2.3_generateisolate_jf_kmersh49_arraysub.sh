#!/bin/bash -l
#$ -l tmem=20G
#$ -l h_vmem=20G
#$ -l h_rt=10:00:0
#$ -N h49keykmerstep4
#$ -V
#$ -j y
#$ -cwd
#$ -hold_jid h49keykmerstep3
#$ -e /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/logs
#$ -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/logs
#$ -t 1-49

# === Activate conda ===
source ~/miniconda3/bin/activate phylogeny_snp

# === Config ===
bamdir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h49/Ps_bams_maptoOTU5_with_haplotype_h49_dedup/h49_bams_softlink"
outdir="$bamdir/OTU5withHTFTFA_fastq_rmdup"
wd="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers"

KMER=31
JFSIZE=300000000

jf_dir="$wd/step2_isolate_kmers/isolate_jf"
dump_dir="$jf_dir/dump"
mkdir -p "$jf_dir" "$dump_dir" logs

# === Build list ===
listfile="$wd/h49_fastq_list.txt"
ls -1 $outdir/*.rmdup.inclHTF.fastq.gz > "$listfile"

fq=$(sed -n "${SGE_TASK_ID}p" "$listfile")
if [ ! -s "$fq" ]; then
    echo "❌ File missing: $fq"
    exit 1
fi

# === Extract clean sample name ===
sample=$(basename "$fq" .rmdup.inclHTF.fastq.gz)
echo "[$(date)] Processing sample: $sample"

jf_out="$jf_dir/${sample}.jf"
dump_out="$dump_dir/${sample}_dump.txt"

# === Run jellyfish count safely ===
zcat "$fq" | jellyfish count -m "$KMER" -s "$JFSIZE" -C -o "$jf_out" /dev/fd/0
exit_code=$?
if [ $exit_code -ne 0 ]; then
    echo "❌ jellyfish count failed for $sample"
    exit 1
fi

# === Dump results ===
jellyfish dump -c "$jf_out" > "$dump_out"
if [ $? -ne 0 ]; then
    echo "❌ jellyfish dump failed for $sample"
    exit 1
fi

echo "✅ Done: $sample"

