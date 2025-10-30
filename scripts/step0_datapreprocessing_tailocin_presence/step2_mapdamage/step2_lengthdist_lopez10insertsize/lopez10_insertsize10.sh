#!/bin/bash -l

#$ -l tmem=5G
#$ -l h_vmem=5G
#$ -l h_rt=24:00:0
#$ -wd /SAN/ugi/plant_genom/jiajucui/
#$ -V
#$ -N slopez10
#$ -t 1-10
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/
echo "Task id is $SGE_TASK_ID"

source ~/miniconda3/bin/activate phylogeny_snp

# ================================
# Inputs
# ================================
i=${SGE_TASK_ID}
sample_list=/SAN/ugi/plant_genom/jiajucui/1_initial_data/new_sequences/h10.txt
samplename=$(sed -n ${i}p "${sample_list}")

refAt=/SAN/ugi/plant_genom/jiajucui/1_initial_data/reference_genome_At/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
refPs=/SAN/ugi/plant_genom/jiajucui/1_initial_data/reference_genome_Ps/Pseudomonas.OTU5_ref.fasta

r1o=/SAN/ugi/plant_genom/jiajucui/1_initial_data/new_sequences/${samplename}.R1.fastq.gz
r2o=/SAN/ugi/plant_genom/jiajucui/1_initial_data/new_sequences/${samplename}.R2.fastq.gz

# ================================
# Outputs
# ================================

base_out=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step1_h49metadata/step2_lengthdist_git_andnsqs10insertsize
outdir_prop=${base_out}/tmreadprop
outdir_ins=${base_out}/insert_sizes

#trimdir=${base_out}/2_trimmed_merged
trimdir=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h49/all_raw_meta_fastq_h49/ 
#outdir_qc=${base_out}/3_quality_control
outdir_At=${base_out}/4_mapping_to_A_thaliana
outdir_Ps=${base_out}/4_mapping_to_pseudomonas
tmp=${base_out}/tmp
#mkdir -p $base_out "$trimdir" "$outdir_qc" "$outdir_At" "$outdir_Ps" "$outdir_prop" "$outdir_ins"
mkdir -p $base_out "$outdir_At" "$outdir_Ps" "$outdir_prop" "$outdir_ins" $tmp
echo "=== ${samplename} ==="

# ================================
# STEP 1: Trim/Merge + FastQC
# ================================
#AdapterRemoval --file1 ${r1o} --file2 ${r2o} \
#  --basename ${trimdir}/${samplename} \
#  --trimns --trimqualities --collapse --gzip --threads 2
#/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h49/all_raw_meta_fastq_h49/ 
#fastqc -t 2 -o ${outdir_qc} \
#  ${trimdir}/${samplename}.collapsed.gz \
#  ${trimdir}/${samplename}.pair1.truncated.gz \
#  ${trimdir}/${samplename}.pair2.truncated.gz

# Collapsed proportion (tmreads)
trimedandmerged=$(zless ${trimdir}/${samplename}.collapsed.gz | grep -c '^+$' || echo 0)
allreadsfq=$(zless ${r1o} | grep -c '^+$' || echo 0)
res4=$(echo "scale=5; ($allreadsfq>0)?${trimedandmerged}/${allreadsfq}:0" | bc)
if [[ ! -s ${outdir_prop}/tmreads.txt ]]; then
  echo -e "Sample\tCollapsedReads\tR1Reads\tCollapsedProp" > ${outdir_prop}/tmreads.txt
fi
echo -e "${samplename}\t${trimedandmerged}\t${allreadsfq}\t0${res4}" >> ${outdir_prop}/tmreads.txt

# ================================
# STEP 2: Map to A. thaliana (collapsed SE and paired PE)
# ================================
# Collapsed → At
bwa aln -t 2 -l 1024 $refAt ${trimdir}/${samplename}.collapsed.gz > ${outdir_At}/${samplename}.collapsed.sai
bwa samse $refAt ${outdir_At}/${samplename}.collapsed.sai ${trimdir}/${samplename}.collapsed.gz \
  | samtools view -b - \
  | samtools sort -o ${outdir_At}/${samplename}.collapsed_At.sorted.bam
samtools view -@ 2 -F 4 -q 20 -bh -o ${outdir_At}/${samplename}.collapsed_At.q20.bam ${outdir_At}/${samplename}.collapsed_At.sorted.bam

samtools sort -o $tmp/${samplename}.collapsed_At.q20.sorted.bam ${outdir_At}/${samplename}.collapsed_At.q20.bam
 
samtools rmdup $tmp/${samplename}.collapsed_At.q20.sorted.bam $tmp/${samplename}.collapsed_At.q20.sorted.rmdup.bam 

samtools sort -o ${outdir_At}/${samplename}.collapsed_At.q20.rmdup.sorted.bam ${tmp}/${samplename}.collapsed_At.q20.sorted.rmdup.bam 

samtools index ${outdir_At}/${samplename}.collapsed_At.q20.rmdup.sorted.bam

# Paired → At
bwa aln -t 2 -l 1024 $refAt ${trimdir}/${samplename}.pair1.truncated.gz > ${outdir_At}/${samplename}.r1.sai
bwa aln -t 2 -l 1024 $refAt ${trimdir}/${samplename}.pair2.truncated.gz > ${outdir_At}/${samplename}.r2.sai
bwa sampe $refAt ${outdir_At}/${samplename}.r1.sai ${outdir_At}/${samplename}.r2.sai \
  ${trimdir}/${samplename}.pair1.truncated.gz ${trimdir}/${samplename}.pair2.truncated.gz \
  | samtools view -b - \
  | samtools sort -o ${outdir_At}/${samplename}.r1r2_At.sorted.bam
samtools view -@ 2 -F 4 -q 20 -bh -o ${outdir_At}/${samplename}.r1r2_At.q20.bam ${outdir_At}/${samplename}.r1r2_At.sorted.bam

# Sort for rmdup
samtools sort -@ 2 \
  -o $tmp/${samplename}.r1r2_At.q20.sorted.bam \
  ${outdir_At}/${samplename}.r1r2_At.q20.bam

# Remove duplicates (paired mode)
samtools rmdup \
  $tmp/${samplename}.r1r2_At.q20.sorted.bam \
  $tmp/${samplename}.r1r2_At.q20.sorted.rmdup.bam

# Final sort and index
samtools sort -@ 2 \
  -o ${outdir_At}/${samplename}.r1r2_At.q20.rmdup.sorted.bam \
  $tmp/${samplename}.r1r2_At.q20.sorted.rmdup.bam

samtools index ${outdir_At}/${samplename}.r1r2_At.q20.rmdup.sorted.bam

# Clean At intermediates
rm -f \
  ${outdir_At}/${samplename}.collapsed.sai \
  ${outdir_At}/${samplename}.r1.sai \
  ${outdir_At}/${samplename}.r2.sai \
  ${outdir_At}/${samplename}.collapsed_At.q20.bam \
  ${outdir_At}/${samplename}.r1r2_At.q20.bam \
  $tmp/${samplename}.collapsed_At.q20.sorted.bam \
  $tmp/${samplename}.collapsed_At.q20.sorted.rmdup.bam \
  $tmp/${samplename}.r1r2_At.q20.sorted.bam \
  $tmp/${samplename}.r1r2_At.q20.sorted.rmdup.bam
# ================================
# STEP 3: Extract unmapped from At and map to Pseudomonas
# ================================

#paried unmapped and collapsed unmapped to ps seperately... then combine...
#logic is first map to At after collapse or truncated (from two libraries), then map to ps from (collapsed after rm At, and truncated after rm At).
# Paired unmapped → Ps (re-map from BAM)
# Extract unmapped paired reads from At BAM
samtools view -bf 12 ${outdir_At}/${samplename}.r1r2_At.sorted.bam > ${outdir_Ps}/${samplename}.r1r2_unmapped.bam
#	Use -bf 12 (both reads unmapped) instead of -bf 4 (only R1 unmapped but R2 could be mapped).
# Convert to FASTQ mates
samtools fastq -@ 4 -n \
  -1 ${outdir_Ps}/${samplename}.r1r2_unmapped.R1.fq.gz \
  -2 ${outdir_Ps}/${samplename}.r1r2_unmapped.R2.fq.gz \
  -0 /dev/null -s /dev/null \
  ${outdir_Ps}/${samplename}.r1r2_unmapped.bam

# Align each mate separately
bwa aln -t 2 -l 1024 $refPs ${outdir_Ps}/${samplename}.r1r2_unmapped.R1.fq.gz > ${outdir_Ps}/${samplename}.R1.sai
bwa aln -t 2 -l 1024 $refPs ${outdir_Ps}/${samplename}.r1r2_unmapped.R2.fq.gz > ${outdir_Ps}/${samplename}.R2.sai

# Proper sampe: requires R1.sai R2.sai + R1.fq R2.fq
bwa sampe $refPs \
  ${outdir_Ps}/${samplename}.R1.sai \
  ${outdir_Ps}/${samplename}.R2.sai \
  ${outdir_Ps}/${samplename}.r1r2_unmapped.R1.fq.gz \
  ${outdir_Ps}/${samplename}.r1r2_unmapped.R2.fq.gz \
| samtools view -@ 4 -F 4 -q 20 -bh - \
| samtools sort -o ${outdir_Ps}/${samplename}.r1r2_Ps.q20.sorted.bam


# --- Deduplicate (paired-end mode)
samtools rmdup \
  ${outdir_Ps}/${samplename}.r1r2_Ps.q20.sorted.bam \
  ${outdir_Ps}/${samplename}.r1r2_Ps.q20.rmdup.bam

# --- Final sort and index
samtools sort -@ 4 \
  -o ${outdir_Ps}/${samplename}.r1r2_Ps.q20.rmdup.sorted.bam \
  ${outdir_Ps}/${samplename}.r1r2_Ps.q20.rmdup.bam

samtools index ${outdir_Ps}/${samplename}.r1r2_Ps.q20.rmdup.sorted.bam

# Cleanup intermediates to save space
rm ${outdir_Ps}/${samplename}.r1r2_unmapped.bam \
   ${outdir_Ps}/${samplename}.r1r2_unmapped.R1.fq.gz \
   ${outdir_Ps}/${samplename}.r1r2_unmapped.R2.fq.gz \
   ${outdir_Ps}/${samplename}.R1.sai \
   ${outdir_Ps}/${samplename}.R2.sai \
   ${outdir_Ps}/${samplename}.r1r2_Ps.q20.sorted.bam \
   ${outdir_Ps}/${samplename}.r1r2_Ps.q20.rmdup.bam


# Collapsed unmapped → Ps (re-map from BAM)
samtools view -bf 4 ${outdir_At}/${samplename}.collapsed_At.sorted.bam > ${outdir_Ps}/${samplename}.collapsed_unmapped.bam
bwa aln -t 2 -l 1024 $refPs -b ${outdir_Ps}/${samplename}.collapsed_unmapped.bam > ${outdir_Ps}/${samplename}.collapsed_unmapped.sai
bwa samse $refPs ${outdir_Ps}/${samplename}.collapsed_unmapped.sai ${outdir_Ps}/${samplename}.collapsed_unmapped.bam \
  | samtools view -b - \
  | samtools sort -o ${outdir_Ps}/${samplename}.collapsed_Ps.sorted.bam
samtools view -@ 2 -F 4 -q 20 -bh -o ${outdir_Ps}/${samplename}.collapsed_Ps.q20.bam ${outdir_Ps}/${samplename}.collapsed_Ps.sorted.bam


# --- 2. Sort for rmdup
samtools sort -@ 2 \
  -o $tmp/${samplename}.collapsed_Ps.q20.sorted.bam \
  ${outdir_Ps}/${samplename}.collapsed_Ps.q20.bam

# --- 3. Remove duplicates (single-end mode)
samtools rmdup -s \
  $tmp/${samplename}.collapsed_Ps.q20.sorted.bam \
  $tmp/${samplename}.collapsed_Ps.q20.sorted.rmdup.bam

# --- 4. Final sort and index
samtools sort -@ 2 \
  -o ${outdir_Ps}/${samplename}.collapsed_Ps.q20.rmdup.sorted.bam \
  $tmp/${samplename}.collapsed_Ps.q20.sorted.rmdup.bam

samtools index ${outdir_Ps}/${samplename}.collapsed_Ps.q20.rmdup.sorted.bam
# --- 5. Cleanup intermediates
rm -f \
  ${outdir_Ps}/${samplename}.collapsed_Ps.q20.bam \
  ${outdir_Ps}/${samplename}.collapsed_Ps.sorted.bam \
  $tmp/${samplename}.collapsed_Ps.q20.sorted.bam \
  $tmp/${samplename}.collapsed_Ps.q20.sorted.rmdup.bam

# Clean Ps intermediates
rm ${outdir_Ps}/${samplename}.collapsed_unmapped.sai

# ================================
# STEP 4: Molecular size distributions (q20 only)
#  - Collapsed: read length from BAM SEQ
#  - Paired: TLEN (absolute insert size) from q20 paired BAM
# ================================
# Headers once

#{outdir_At}/${samplename}.collapsed_At.q20.rmdup.sorted.bam
#{outdir_At}/${samplename}.r1r2_At.q20.rmdup.sorted.bam
#{outdir_Ps}/${samplename}.r1r2_Ps.q20.rmdup.sorted.bam
#{outdir_Ps}/${samplename}.collapsed_Ps.q20.rmdup.sorted.bam


if [[ ! -s ${outdir_ins}/ALL_lengths.tsv ]]; then
  echo -e "Length\tCount\tSample\tType\tSpecies" > ${outdir_ins}/ALL_lengths.tsv
fi

# ---- A. thaliana collapsed ----
samtools view ${outdir_At}/${samplename}.collapsed_At.q20.rmdup.sorted.bam | cut -f 10 \
  | awk '{print length($0)}' | sort -n | uniq -c \
  | awk -v s=${samplename} '{print $2"\t"$1"\t"s"\tCollapsed\tArabidopsis"}' \
  > ${outdir_ins}/${samplename}_At_collapsed_length.txt

# ---- A. thaliana paired (TLEN) ----
samtools view -f 66 ${outdir_At}/${samplename}.r1r2_At.q20.rmdup.sorted.bam | cut -f 9 \
  | awk '{if($1<0) $1=-$1} {c[$1]++} END{for (k in c) print k"\t"c[k]}' \
  | sort -k1,1n \
  | awk -v s=${samplename} '{print $1"\t"$2"\t"s"\tPaired\tArabidopsis"}' \
  > ${outdir_ins}/${samplename}_At_paired_length.txt

# ---- Pseudomonas collapsed ----
samtools view ${outdir_Ps}/${samplename}.collapsed_Ps.q20.rmdup.sorted.bam | cut -f 10 \
  | awk '{print length($0)}' | sort -n | uniq -c \
  | awk -v s=${samplename} '{print $2"\t"$1"\t"s"\tCollapsed\tPseudomonas"}' \
  > ${outdir_ins}/${samplename}_Ps_collapsed_length.txt

# ---- Pseudomonas paired (TLEN) ----
samtools view -f 66 ${outdir_Ps}/${samplename}.r1r2_Ps.q20.rmdup.sorted.bam | cut -f 9 \
  | awk '{if($1<0) $1=-$1} {c[$1]++} END{for (k in c) print k"\t"c[k]}' \
  | sort -k1,1n \
  | awk -v s=${samplename} '{print $1"\t"$2"\t"s"\tPaired\tPseudomonas"}' \
  > ${outdir_ins}/${samplename}_Ps_paired_length.txt

# Per-sample all lengths
cat ${outdir_ins}/${samplename}_*_length.txt > ${outdir_ins}/${samplename}_all_lengths.txt
# Append to global (safe even if repeated runs; dedup later if needed)
cat ${outdir_ins}/${samplename}_*_length.txt >> ${outdir_ins}/ALL_lengths.tsv

# ================================
# STEP 5: Mean insert sizes (paired q20)
# ================================
if [[ ! -s ${outdir_ins}/insert_size.txt ]]; then
  echo -e "Sample\tSpecies\tMeanInsert" > ${outdir_ins}/insert_size.txt
fi

meanAt=$(samtools view -f 66 ${outdir_At}/${samplename}.r1r2_At.q20.rmdup.sorted.bam | cut -f 9 \
         | awk '{if($1<0) $1=-$1; sum+=$1} END{if(NR>0) print sum/NR; else print 0}')
meanPs=$(samtools view -f 66 ${outdir_Ps}/${samplename}.r1r2_Ps.q20.rmdup.sorted.bam | cut -f 9 \
         | awk '{if($1<0) $1=-$1; sum+=$1} END{if(NR>0) print sum/NR; else print 0}')

echo -e "${samplename}\tArabidopsis\t${meanAt}" >> ${outdir_ins}/insert_size.txt
echo -e "${samplename}\tPseudomonas\t${meanPs}" >> ${outdir_ins}/insert_size.txt

# ================================
# STEP 6: Join tmreads with mean insert sizes (per run append)
# ================================
# Build a two-column wide table MeanInsert_At and MeanInsert_Ps
awk 'BEGIN{FS=OFS="\t"} NR>1{m[$1 FS $2]=$3} END{for(k in m) print k, m[k]}' ${outdir_ins}/insert_size.txt \
 | awk -F'\t' 'BEGIN{OFS="\t"} $2=="Arabidopsis"{a[$1]=$3} $2=="Pseudomonas"{p[$1]=$3} END{print "Sample","MeanInsert_At","MeanInsert_Ps"; for(s in a) print s, a[s], (p[s]==""?0:p[s])}' \
 > ${outdir_ins}/mean_insert_wide.tmp

join -t $'\t' -1 1 -2 1 <(sort -k1,1 ${outdir_prop}/tmreads.txt) <(sort -k1,1 ${outdir_ins}/mean_insert_wide.tmp) \
  > ${outdir_ins}/tmreads_insert_merged.txt
rm -f ${outdir_ins}/mean_insert_wide.tmp

# ================================
# STEP 7: Remove unnecessary intermediates
#  (keeps both *.sorted.bam and *.q20.bam + .bai as requested)
# ================================
# At: we already removed *.sai; no SAM kept.
# Ps: remove temporary unmapped BAMs and *.sai already done above.

echo "Done: ${samplename}"
