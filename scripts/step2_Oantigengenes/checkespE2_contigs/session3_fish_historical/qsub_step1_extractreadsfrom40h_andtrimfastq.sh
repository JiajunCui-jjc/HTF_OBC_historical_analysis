#!/bin/bash -l

#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=24:00:0
#$ -wd /SAN/ugi/plant_genom/jiajucui/
#$ -V
#$ -N s2025_h49_extractfastq
#$ -t 1-49
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/

echo "Task id is $SGE_TASK_ID"
#variables
i=$SGE_TASK_ID
source ~/miniconda3/bin/activate phylogeny_snp



# ======= User Config =======
bam_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/2025_h49"
outbase="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/fish_historical/using39m_espE2fasta"
espE2_ref="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/espE2_rescue/checkespE2_57m/final_afterexpasy/completeopenframe.39fasta.txt"
threads=4
mkdir -p $outbase
# ======= Output dirs =======
mkdir -p "$outbase"/{fastq_from_bam,mapped_bam_espE2,fastq_espE2only,tmp}

# ======= BAM List + Select Sample =======
bam_list=($(ls "$bam_dir"/*_after_removal_mappedAt.bam | sort))
bam="${bam_list[$((SGE_TASK_ID-1))]}"
sample=$(basename "$bam" _after_removal_mappedAt.bam)

echo "🔍 Processing: $sample"

# ======= Step 1: BAM → FASTQ =======
fastq_raw="$outbase/fastq_from_bam/${sample}.fastq"
samtools fastq "$bam" > "$fastq_raw"

# ======= Step 2: bwa aln + samse + sort =======
sai_tmp="$outbase/tmp/${sample}.sai"
bam_out="$outbase/mapped_bam_espE2/${sample}.espE2.bam"
bwa aln -t $threads "$espE2_ref" "$fastq_raw" > "$sai_tmp"
bwa samse "$espE2_ref" "$sai_tmp" "$fastq_raw" | samtools sort -@ $threads -o "$bam_out"
samtools index "$bam_out"

# ======= Step 3: Extract mapped reads only =======
fastq_out="$outbase/fastq_espE2only/${sample}_espE2mapped.fastq"
samtools fastq -F 4 "$bam_out" > "$fastq_out"

echo "✅ Extracted reads to: $fastq_out"

# ======= Final Cleanup =======
rm -f "$sai_tmp" "$fastq_raw"  


