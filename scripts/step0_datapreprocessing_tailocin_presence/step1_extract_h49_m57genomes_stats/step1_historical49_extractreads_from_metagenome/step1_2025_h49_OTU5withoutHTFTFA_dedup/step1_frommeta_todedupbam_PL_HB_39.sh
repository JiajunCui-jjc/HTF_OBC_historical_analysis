#!/bin/bash -l

#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=24:00:0
#$ -wd /SAN/ugi/plant_genom/jiajucui/
#$ -V
#$ -N OTU5refs2025_h39
#$ -t 1-39
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/

echo "Task id is $SGE_TASK_ID"
#variables
i=$SGE_TASK_ID
source ~/miniconda3/bin/activate phylogeny_snp
#variables
#need to modify: reference refPS,
#sample name and the path to store the answer
#thecollapsed file sotrage
#the path of intermediate files: ${pathtmp}
pathtmp=2025_maptoOTU5_nohaplotype_h49
answerspath=answer2025_maptoOTU5_nohaplotype_h49
mkdir -p /SAN/ugi/plant_genom/jiajucui/${answerspath}/
#mkdir -p /SAN/ugi/plant_genom/jiajucui/3_quality_control/${pathtmp}/
samplename=$(grep -E '^(HB|PL)' /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/historical49.txt | sed -n ${i}p)
echo "${samplename}:" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt
refAt=/SAN/ugi/plant_genom/jiajucui/1_initial_data/reference_genome_At/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
#refPs=/SAN/ugi/plant_genom/jiajucui/1_initial_data/reference_genome_Ps_with_tailocin_haplotypes/Pseudomonas.plate25.C2.pilon.contigs_renamed.with_Tail_Fiber_Haps.fasta
refPs="/SAN/ugi/plant_genom/jiajucui/1_initial_data/reference_genome_Ps/Pseudomonas.OTU5_ref.fasta"

#for PL samples: we start with raw r1 and r2
#r1o="$wd/data/all_fastq_h49/all_raw_meta_fastq_h49/${samplename}.R1.fastq.gz"
#r2o="$wd/data/all_fastq_h49/all_raw_meta_fastq_h49/${samplename}.R2.fastq.gz"


#collapse:
#trim and merge
#AdapterRemoval --file1 ${r1o} --file2 ${r2o} --basename /SAN/ugi/plant_genom/jiajucui/2_trimmed_merged/${pathtmp}/${samplename} --trimns --trimqualities --collapse --gzip --threads 2 &&
#echo 'trim and merge done' &&
#fastqc
#fastqc -t 2 -o /SAN/ugi/plant_genom/jiajucui/3_quality_control/${pathtmp}/ /SAN/ugi/plant_genom/jiajucui/2_trimmed_merged/${pathtmp}/${samplename}.collapsed.gz /SAN/ugi/plant_genom/jiajucui/2_trimmed_merged/${pathtmp}/${samplename}.pair1.truncated.gz /SAN/ugi/plant_genom/jiajucui/2_trimmed_merged/${pathtmp}/${samplename}.pair2.truncated.gz



collapsed=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h49/all_raw_meta_fastq_h49/${samplename}.fastq.gz



mkdir -p /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/
out_dir=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}

# ===========================================================
# q1 What is the percentage of A.thaliana DNA?
# ===========================================================

bwa aln -t 8 -l 1024 -f ${out_dir}/${samplename}.collapsed.sai ${refAt} ${collapsed}

bwa samse -r @RG\\tID:${samplename}\\tSM:${samplename} \
  -f ${out_dir}/${samplename}.collapsed_At.sam ${refAt} ${out_dir}/${samplename}.collapsed.sai ${collapsed}

samtools flagstat ${out_dir}/${samplename}.collapsed_At.sam > ${out_dir}/${samplename}_Atsam_flagstats.log

echo "q1 What is the percentage of A.thaliana DNA?" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt
echo "reads in total (allrawreads):" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt
grep 'in total' ${out_dir}/${samplename}_Atsam_flagstats.log >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt

# convert mapped reads to BAM with mapQ≥20
samtools view -@ 4 -F 4 -q 20 -Sbh -o ${out_dir}/${samplename}_mapped_At_q20.bam ${out_dir}/${samplename}.collapsed_At.sam
rm ${out_dir}/${samplename}.collapsed.sai

samtools flagstat ${out_dir}/${samplename}_mapped_At_q20.bam > ${out_dir}/${samplename}_At_q20_flagstats.log
echo "reads mapped to At (-q20):" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt
grep 'mapped (' ${out_dir}/${samplename}_At_q20_flagstats.log >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt


# ===========================================================
# q1.1 Covered genome proportion (after markdup removal)
# ===========================================================
echo "q1.1 What is the covered genome proportion of A.thaliana DNA?" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt

# Remove duplicates and sort
samtools sort -o ${out_dir}/${samplename}_mapped_At_q20.sorted.bam ${out_dir}/${samplename}_mapped_At_q20.bam
samtools markdup -r -@ 4 ${out_dir}/${samplename}_mapped_At_q20.sorted.bam ${out_dir}/${samplename}_mapped_At_q20.sorted.markdup.bam
samtools sort -o ${out_dir}/${samplename}_mapped_At_q20.markdup.sorted.bam ${out_dir}/${samplename}_mapped_At_q20.sorted.markdup.bam
samtools index ${out_dir}/${samplename}_mapped_At_q20.markdup.sorted.bam
rm ${out_dir}/${samplename}_mapped_At_q20.markdup.bam ${out_dir}/${samplename}_mapped_At_q20.bam

samtools coverage ${out_dir}/${samplename}_mapped_At_q20.markdup.sorted.bam > ${out_dir}/${samplename}_mapped_At_q20.markdup.sorted.bysamtoolscoverage.txt

echo "total length of At ref:" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt
awk '{sum+=$3;} END{print sum;}' ${out_dir}/${samplename}_mapped_At_q20.markdup.sorted.bysamtoolscoverage.txt >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt

echo "base mapped to At ref:" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt
awk '{sum+=$5;} END{print sum;}' ${out_dir}/${samplename}_mapped_At_q20.markdup.sorted.bysamtoolscoverage.txt >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt


# ===========================================================
# q1.2 Average depth
# ===========================================================
echo "q1.2 What is the read depth of A.thaliana DNA?" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt
echo "read depth to At ref:" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt
awk '{sum+=($3*$7);} END{print sum/119667750;}' ${out_dir}/${samplename}_mapped_At_q20.markdup.sorted.bysamtoolscoverage.txt >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt


# ===========================================================
# Cleanup intermediate files
# ===========================================================
rm ${out_dir}/${samplename}_mapped_At_q20.markdup.sorted.bysamtoolscoverage.txt


# ===========================================================
# q2: extract unmapped reads to pass to Pseudomonas
# ===========================================================
echo "q2 What is the percentage of Pseudomonas DNA after removing mapped reads to A.thaliana?" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt

# Extract unmapped reads (-bf 4)
samtools view -bf 4 ${out_dir}/${samplename}.collapsed_At.sam > ${out_dir}/${samplename}_after_removal_mappedAt.bam

rm ${out_dir}/${samplename}.collapsed_At.sam

echo "All calculation done. Final BAMs retained:"
echo " - ${samplename}_mapped_At_q20.markdup.sorted.bam"
echo " - ${samplename}_after_removal_mappedAt.bam"






##use bwa to realign the unmapped bam file
bwa aln -t 2 -l 1024 ${refPs} -b /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_after_removal_mappedAt.bam > /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_realign.sai

bwa samse -r @RG\\tID:${samplename}\\tSM:${samplename} -f /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_mapped_to_ps.sam ${refPs} /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_realign.sai /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_after_removal_mappedAt.bam
#Keep the mapped reads and create a compressed BAM file using samtools
samtools view -@ 2 -F 4 -Sbh -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_mapped_to_ps.bam /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_mapped_to_ps.sam
#Sort the BAM file by chromosome and position using samtools
samtools sort -@ 2 -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_mapped_to_ps.sorted.bam /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_mapped_to_ps.bam

samtools view -q 20 /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_mapped_to_ps.sorted.bam -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}.unsort_mapped_to_Pseudomonas.dd.q20.bam
samtools sort -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}.mapped_to_Pseudomonas.dd.q20.bam /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}.unsort_mapped_to_Pseudomonas.dd.q20.bam

#rm sai sam bam and unsort unq20 bams:
rm /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_realign.sai
rm /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_mapped_to_ps.sam
rm /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_mapped_to_ps.sorted.bam




#for Ps
# ==============================
# Paths and variables
# ==============================
base_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}"
in_bam="${base_dir}/${samplename}.mapped_to_Pseudomonas.dd.q20.bam"
markdup_bam="${base_dir}/${samplename}.mapped_to_Pseudomonas.dd.q20.markeddup.bam"
out_answer="/SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt"
nm_dir="${base_dir}/NM225"
mkdir -p "$nm_dir"

# ==============================
# 1. Mapping percentage (pre-dedup)
# ==============================
samtools flagstat "$in_bam" > "${base_dir}/${samplename}_Ps_q20_flagstats.log"
echo "q2 What is the percentage of Pseudomonas DNA after removing mapped reads to A.thaliana?" >> "$out_answer"
grep 'mapped (' "${base_dir}/${samplename}_Ps_q20_flagstats.log" >> "$out_answer"

# ==============================
# 2. Deduplication (sort → markdup → sort + index)
# ==============================
samtools sort -o "${in_bam%.bam}.sorted.bam" "$in_bam"
samtools markdup -r "${in_bam%.bam}.sorted.bam" "${markdup_bam%.bam}.unsorted.bam"
samtools sort -o "$markdup_bam" "${markdup_bam%.bam}.unsorted.bam"
samtools index "$markdup_bam"
rm -f "${markdup_bam%.bam}.unsorted.bam" "${in_bam%.bam}.sorted.bam"

# ==============================
# 3. Coverage and depth (after dedup)
# ==============================
covfile="${base_dir}/${samplename}.mapped_to_Pseudomonas.dd.q20.markeddup.coverage.txt"
samtools depth -aa "$markdup_bam" > "$covfile"

echo "q4 What is the percentage of the Pseudomonas genome that is covered by at least 1 read?" >> "$out_answer"
echo "pseudomonas bases in total:" >> "$out_answer"
wc -l "$covfile" | cut -d ' ' -f 1 >> "$out_answer"
echo "bases with at least 1 read:" >> "$out_answer"
awk '($3 >= 1) {count++} END {print count}' "$covfile" >> "$out_answer"

# ==============================
# 4. NM (edit distance) distribution
# ==============================
samtools view "$markdup_bam" \
  | awk -F '\t' '{for(i=12;i<=NF;i++){if($i ~ /^NM:i:/){split($i,a,":"); print a[3]; break}}}' \
  | sort | uniq -c | sort -k2,2n | awk '{print $2"\t"$1}' \
  > "${nm_dir}/${samplename}_NM.txt"

# ==============================
# 5. Average depth
# ==============================
echo 'q6 Average depth of Ps genome?' >> "$out_answer"
len=5941411
awk -v len=$len '{sum+=$3} END {print sum/len}' "$covfile" >> "$out_answer"

# ==============================
# 6. Cleanup
# ==============================
rm -f "$covfile"
rm -f ${base_dir}/${samplename}_Ps_q20_flagstats.log
rm -f ${base_dir}/${samplename}.unsort_mapped_to_Pseudomonas.dd.q20.bam \
      ${base_dir}/${samplename}_removalAt_mapped_to_ps.sam \
      ${base_dir}/${samplename}_removalAt_mapped_to_ps.sorted.bam \
      ${base_dir}/${samplename}_removalAt_realign.sai
rm -f ${base_dir}/${samplename}_removalAt_mapped_to_ps.bam
echo "All calculation done. Only dedup q20 BAM retained for ${samplename}."

# ==============================
# 7. Symlink dedup BAMs (Lopez10 subset under h49 structure)
# ==============================
At_linkdir=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h49/At_bams_maptoOTU5withouthaplotype_h49_dedup/
Ps_linkdir=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h49/Ps_bams_maptoOTU5withouthaplotype_h49_dedup/
mkdir -p "$At_linkdir" "$Ps_linkdir"

ln -sf /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_mapped_At_q20.markdup.sorted.bam* "$At_linkdir/"
ln -sf ${base_dir}/${samplename}.mapped_to_Pseudomonas.dd.q20.markeddup.bam* "$Ps_linkdir/"

rm -f /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_mapped_At_q20.sorted.markdup.bam

echo "✅ Symlinks created under h49 structure (mapped to OTU5 without haplotype) for ${samplename}"
