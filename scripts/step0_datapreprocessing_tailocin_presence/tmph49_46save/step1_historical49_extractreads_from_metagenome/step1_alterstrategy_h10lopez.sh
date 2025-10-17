#!/bin/bash -l

#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=24:00:0
#$ -wd /SAN/ugi/plant_genom/jiajucui/
#$ -V
#$ -N lopez10
#$ -t 1-10
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/

echo "Task id is $SGE_TASK_ID"
#variables
i=$SGE_TASK_ID
source ~/miniconda3/bin/activate phylogeny_snp
# ===========================================================
# 1. Variables and paths
# ===========================================================
i=${SGE_TASK_ID}
pathtmp=2025_h49/h10_lopez
basePs=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}
baseAt=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}
answerspath=/SAN/ugi/plant_genom/jiajucui/answer2025_h49/h10_lopez

mkdir -p "$basePs" "$baseAt" "$answerspath"

sample_list=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/historical49.txt
samplename=$(grep -vE '^(HB|PL)' "${sample_list}" | sed -n ${i}p)

out_answer=${answerspath}/answers_for_${samplename}.txt
echo "${samplename}:" > "$out_answer"

refAt=/SAN/ugi/plant_genom/jiajucui/1_initial_data/reference_genome_At/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
refPs=/SAN/ugi/plant_genom/jiajucui/1_initial_data/reference_genome_Ps_with_tailocin_haplotypes/Pseudomonas.plate25.C2.pilon.contigs_renamed.with_Tail_Fiber_Haps.fasta

#r1o="$wd/data/all_fastq_h49/all_raw_meta_fastq_h49/${samplename}.R1.fastq.gz"
#r2o="$wd/data/all_fastq_h49/all_raw_meta_fastq_h49/${samplename}.R2.fastq.gz"


#collapse:
#trim and merge
#AdapterRemoval --file1 ${r1o} --file2 ${r2o} --basename /SAN/ugi/plant_genom/jiajucui/2_trimmed_merged/${pathtmp}/${samplename} --trimns --trimqualities --collapse --gzip --threads 2 &&
#echo 'trim and merge done' &&
#fastqc
#fastqc -t 2 -o /SAN/ugi/plant_genom/jiajucui/3_quality_control/${pathtmp}/ /SAN/ugi/plant_genom/jiajucui/2_trimmed_merged/${pathtmp}/${samplename}.collapsed.gz /SAN/ugi/plant_genom/jiajucui/2_trimmed_merged/${pathtmp}/${samplename}.pair1.truncated.gz /SAN/ugi/plant_genom/jiajucui/2_trimmed_merged/${pathtmp}/${samplename}.pair2.truncated.gz

raw=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h49/all_raw_meta_fastq_h49/
collapsed=$raw/${samplename}.collapsed.gz
r1t=$raw/${samplename}.pair1.truncated.gz
r2t=$raw/${samplename}.pair2.truncated.gz

# ===========================================================
# 2. Map to A.thaliana (collapsed + paired)
# ===========================================================
bwa aln -t 4 -l 1024 $refAt $collapsed > ${baseAt}/${samplename}.collapsed.sai
bwa samse $refAt ${baseAt}/${samplename}.collapsed.sai $collapsed \
  | samtools view -b - | samtools sort -o ${baseAt}/${samplename}.collapsed_At.sorted.bam
samtools index ${baseAt}/${samplename}.collapsed_At.sorted.bam

bwa aln -t 4 -l 1024 $refAt $r1t > ${baseAt}/${samplename}.r1.sai
bwa aln -t 4 -l 1024 $refAt $r2t > ${baseAt}/${samplename}.r2.sai
bwa sampe $refAt ${baseAt}/${samplename}.r1.sai ${baseAt}/${samplename}.r2.sai $r1t $r2t \
  | samtools view -b - | samtools sort -o ${baseAt}/${samplename}.r1r2_At.sorted.bam
samtools index ${baseAt}/${samplename}.r1r2_At.sorted.bam

# merge both
samtools merge -n ${baseAt}/${samplename}.combined_At.bam \
  ${baseAt}/${samplename}.collapsed_At.sorted.bam \
  ${baseAt}/${samplename}.r1r2_At.sorted.bam
samtools sort -o ${baseAt}/${samplename}.combined_At.sorted.bam ${baseAt}/${samplename}.combined_At.bam
samtools index ${baseAt}/${samplename}.combined_At.sorted.bam

samtools flagstat ${baseAt}/${samplename}.combined_At.sorted.bam > ${baseAt}/${samplename}_Atsam_flagstats.log
echo "q1 What is the percentage of A.thaliana DNA?" >> "$out_answer"
echo "reads in total (allrawreads):" >> "$out_answer"
grep 'in total' ${baseAt}/${samplename}_Atsam_flagstats.log >> "$out_answer"

# Dedup for A.thaliana (sort → rmdup)
samtools view -@ 4 -F 4 -q 20 -bh -o ${baseAt}/${samplename}_mapped_At_q20.bam ${baseAt}/${samplename}.combined_At.sorted.bam
samtools sort -o ${baseAt}/${samplename}_mapped_At_q20.sorted.bam ${baseAt}/${samplename}_mapped_At_q20.bam
samtools rmdup ${baseAt}/${samplename}_mapped_At_q20.sorted.bam ${baseAt}/${samplename}_mapped_At_q20.rmdup.bam
samtools sort -o ${baseAt}/${samplename}_mapped_At_q20.rmdup.sorted.bam ${baseAt}/${samplename}_mapped_At_q20.rmdup.bam
samtools index ${baseAt}/${samplename}_mapped_At_q20.rmdup.sorted.bam
rm ${baseAt}/${samplename}_mapped_At_q20.bam ${baseAt}/${samplename}_mapped_At_q20.rmdup.bam

# ===========================================================
# 3. Extract unmapped → Pseudomonas (collapsed + paired)
# ===========================================================
# collapsed
samtools view -bf 4 ${baseAt}/${samplename}.collapsed_At.sorted.bam > ${basePs}/${samplename}.collapsed_unmapped.bam
bwa aln -t 4 -l 1024 $refPs -b ${basePs}/${samplename}.collapsed_unmapped.bam > ${basePs}/${samplename}.collapsed_Ps.sai
bwa samse $refPs ${basePs}/${samplename}.collapsed_Ps.sai ${basePs}/${samplename}.collapsed_unmapped.bam \
  | samtools view -b - | samtools sort -o ${basePs}/${samplename}.collapsed_Ps.sorted.bam
samtools view -@ 4 -F 4 -q 20 -bh -o ${basePs}/${samplename}.collapsed_Ps.q20.bam ${basePs}/${samplename}.collapsed_Ps.sorted.bam
samtools index ${basePs}/${samplename}.collapsed_Ps.q20.bam

# paired
samtools view -bf 12 ${baseAt}/${samplename}.r1r2_At.sorted.bam > ${basePs}/${samplename}.r1r2_unmapped.bam
samtools fastq -@ 4 -n \
  -1 ${basePs}/${samplename}.r1r2_unmapped.R1.fq.gz \
  -2 ${basePs}/${samplename}.r1r2_unmapped.R2.fq.gz \
  -0 /dev/null -s /dev/null \
  ${basePs}/${samplename}.r1r2_unmapped.bam
bwa aln -t 4 -l 1024 $refPs ${basePs}/${samplename}.r1r2_unmapped.R1.fq.gz > ${basePs}/${samplename}.R1.sai
bwa aln -t 4 -l 1024 $refPs ${basePs}/${samplename}.r1r2_unmapped.R2.fq.gz > ${basePs}/${samplename}.R2.sai
bwa sampe $refPs \
  ${basePs}/${samplename}.R1.sai ${basePs}/${samplename}.R2.sai \
  ${basePs}/${samplename}.r1r2_unmapped.R1.fq.gz ${basePs}/${samplename}.r1r2_unmapped.R2.fq.gz \
  | samtools view -b - | samtools sort -o ${basePs}/${samplename}.r1r2_Ps.sorted.bam
samtools view -@ 4 -F 4 -q 20 -bh -o ${basePs}/${samplename}.r1r2_Ps.q20.bam ${basePs}/${samplename}.r1r2_Ps.sorted.bam
samtools index ${basePs}/${samplename}.r1r2_Ps.q20.bam

# merge (after indexing both)
samtools merge -n ${basePs}/${samplename}.mapped_to_Pseudomonas.dd.q20.bam \
  ${basePs}/${samplename}.collapsed_Ps.q20.bam \
  ${basePs}/${samplename}.r1r2_Ps.q20.bam
samtools index ${basePs}/${samplename}.mapped_to_Pseudomonas.dd.q20.bam

# ===========================================================
# 4. Pseudomonas stats (exclude, rmdup, cov, NM)
# ===========================================================
in_bam="${basePs}/${samplename}.mapped_to_Pseudomonas.dd.q20.bam"
dedup_bam="${basePs}/${samplename}.mapped_to_Pseudomonas.dd.q20.rmdup.bam"
nm_dir="${basePs}/NM225"
mkdir -p "$nm_dir"

# --- pre-dedup read count ---
bam_excl_predup="${basePs}/${samplename}.mapped_to_Pseudomonas.dd.q20.excluding_HTF_TFA.bam"
samtools view -h "$in_bam" \
  | awk 'BEGIN{OFS="\t"} /^@/ {print; next} {gsub(/\r/,"",$3); ref=toupper($3)} !(ref ~ /^HTF/ || ref ~ /^TFA/)' \
  | samtools view -b -o "$bam_excl_predup"
samtools index "$bam_excl_predup"
samtools flagstat "$bam_excl_predup" > "${basePs}/${samplename}_Ps_q20_excluding_HTF_TFA_flagstats.log"
echo "q2 What is the percentage of Pseudomonas DNA after removing mapped reads to A.thaliana?" >> "$out_answer"
echo "reads mapped to Pseudomonas (-q20) after removing:" >> "$out_answer"
grep 'mapped (' "${basePs}/${samplename}_Ps_q20_excluding_HTF_TFA_flagstats.log" >> "$out_answer"

# --- dedup (sort → rmdup) ---
samtools sort -o "${in_bam%.bam}.sorted.bam" "$in_bam"
samtools rmdup "${in_bam%.bam}.sorted.bam" "$dedup_bam"
samtools sort -o "${dedup_bam%.bam}.sorted.bam" "$dedup_bam"
mv "${dedup_bam%.bam}.sorted.bam" "$dedup_bam"
samtools index "$dedup_bam"

# --- exclude HTF/TFA after dedup for coverage/depth ---
bam_excl_postdup="${basePs}/${samplename}.mapped_to_Pseudomonas.dd.q20.rmdup.excluding_HTF_TFA.bam"
covfile="${basePs}/${samplename}.mapped_to_Pseudomonas.dd.q20.rmdup.excluding_HTF_TFA.coverage.txt"

samtools view -h "$dedup_bam" \
  | awk 'BEGIN{OFS="\t"} /^@/ {print; next} {gsub(/\r/,"",$3); ref=toupper($3)} !(ref ~ /^HTF/ || ref ~ /^TFA/)' \
  | samtools view -b -o "$bam_excl_postdup"
samtools index "$bam_excl_postdup"

# build a list of contigs to keep (exclude HTF/TFA)
contig_list="${basePs}/non_HTF_TFA_contigs.txt"
samtools idxstats "$dedup_bam" \
  | awk -F '\t' '!($1 ~ /^HTF/ || $1 ~ /^TFA/ || $1=="*"){print $1"\t0\t"$2}' \
  > "$contig_list"
echo 'q4 What is the percentage of the Pseudomonas genome that is covered by at least 1 read?' >> "$out_answer"
samtools depth -aa -b "$contig_list" "$bam_excl_postdup" > "$covfile"

echo 'pseudomonas bases in total:' >> "$out_answer"
wc -l "$covfile" | cut -d ' ' -f 1 >> "$out_answer"
echo 'bases with at least 1 read:' >> "$out_answer"
awk '($3 >= 1) {count++} END {print count}' "$covfile" >> "$out_answer"

# --- NM distribution ---
samtools view "$bam_excl_postdup" \
  | awk -F '\t' '{for(i=12;i<=NF;i++){if($i~/^NM:i:/){split($i,a,":");print a[3];break}}}' \
  | sort | uniq -c | sort -k2,2n | awk '{print $2"\t"$1}' \
  > "${nm_dir}/${samplename}_NM.txt"

# --- Average depth ---
echo 'q6 Average depth of Ps genome?' >> "$out_answer"
len=5941411
awk -v len=$len '{sum+=$3} END {print sum/len}' "$covfile" >> "$out_answer"

# --- Cleanup ---
rm -f "$bam_excl_predup" "${bam_excl_predup}.bai" "$bam_excl_postdup" "${bam_excl_postdup}.bai" "$covfile"
rm -f ${basePs}/${samplename}.collapsed_Ps.sai ${basePs}/${samplename}.R1.sai ${basePs}/${samplename}.R2.sai
rm -f ${basePs}/${samplename}.collapsed_unmapped.bam ${basePs}/${samplename}.r1r2_unmapped.bam
rm -f ${basePs}/${samplename}.r1r2_unmapped.R1.fq.gz ${basePs}/${samplename}.r1r2_unmapped.R2.fq.gz
rm -f ${basePs}/${samplename}_Ps_q20_excluding_HTF_TFA_flagstats.log
rm -f ${basePs}/${samplename}_removalAt_mapped_to_ps.bam
rm -f ${basePs}/${samplename}.unsort_mapped_to_Pseudomonas.dd.q20.bam
rm -f ${basePs}/${samplename}.mapped_to_Pseudomonas.dd.q20.bam

#then ln -s
# ===========================================================
# 5. Symlink dedup BAMs (Lopez10 subset under h49 structure)
# ===========================================================
echo "Linking deduplicated BAMs for Lopez10 subset into h49 structure..."

# Target directories (Lopez10 subfolders)
At_linkdir=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h49/At_bams_h49_dedup/h10_lopez
Ps_linkdir=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h49/Ps_bams_h49_dedup/h10_lopez

mkdir -p "$At_linkdir" "$Ps_linkdir"

# Symlink A.thaliana dedup BAM
ln -sf ${baseAt}/${samplename}_mapped_At_q20.rmdup.sorted.bam \
       ${At_linkdir}/${samplename}_mapped_At_q20.rmdup.sorted.bam
ln -sf ${baseAt}/${samplename}_mapped_At_q20.rmdup.sorted.bam.bai \
       ${At_linkdir}/${samplename}_mapped_At_q20.rmdup.sorted.bam.bai

# Symlink Pseudomonas dedup BAM
ln -sf ${basePs}/${samplename}.mapped_to_Pseudomonas.dd.q20.rmdup.bam \
       ${Ps_linkdir}/${samplename}.mapped_to_Pseudomonas.dd.q20.rmdup.bam
ln -sf ${basePs}/${samplename}.mapped_to_Pseudomonas.dd.q20.rmdup.bam.bai \
       ${Ps_linkdir}/${samplename}.mapped_to_Pseudomonas.dd.q20.rmdup.bam.bai

echo "✅ Symlinks created under h49 structure for ${samplename}"
echo "✅ Finished successfully for ${samplename}"
