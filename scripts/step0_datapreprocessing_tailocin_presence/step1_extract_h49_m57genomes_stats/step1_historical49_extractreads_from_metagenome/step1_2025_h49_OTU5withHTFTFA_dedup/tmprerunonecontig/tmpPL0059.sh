#!/bin/bash -l

#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=24:00:0
#$ -wd /SAN/ugi/plant_genom/jiajucui/
#$ -V
#$ -N s2025_h39
#$ -t 6
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
pathtmp=2025_h49
answerspath=answer2025_h49_notused_withHTFTFA
mkdir -p /SAN/ugi/plant_genom/jiajucui/${answerspath}/
#mkdir -p /SAN/ugi/plant_genom/jiajucui/3_quality_control/${pathtmp}/
samplename=$(grep -E '^(HB|PL)' /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/historical49.txt | sed -n ${i}p)
echo "${samplename}:" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt
refAt=/SAN/ugi/plant_genom/jiajucui/1_initial_data/reference_genome_At/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
refPs='/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/HTFreference/OTU5_ref/ref.withhyplotype.fasta'

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

#bwa aln -t 8 -l 1024 -f ${out_dir}/${samplename}.collapsed.sai ${refAt} ${collapsed}

#bwa samse -r @RG\\tID:${samplename}\\tSM:${samplename} \
#  -f ${out_dir}/${samplename}.collapsed_At.sam ${refAt} ${out_dir}/${samplename}.collapsed.sai ${collapsed}

#samtools flagstat ${out_dir}/${samplename}.collapsed_At.sam > ${out_dir}/${samplename}_Atsam_flagstats.log

#echo "q1 What is the percentage of A.thaliana DNA?" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt
#echo "reads in total (allrawreads):" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt
#grep 'in total' ${out_dir}/${samplename}_Atsam_flagstats.log >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt

# convert mapped reads to BAM with mapQ≥20
#samtools view -@ 4 -F 4 -q 20 -Sbh -o ${out_dir}/${samplename}_mapped_At_q20.bam ${out_dir}/${samplename}.collapsed_At.sam
#rm ${out_dir}/${samplename}.collapsed.sai

#samtools flagstat ${out_dir}/${samplename}_mapped_At_q20.bam > ${out_dir}/${samplename}_At_q20_flagstats.log
#echo "reads mapped to At (-q20):" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt
#grep 'mapped (' ${out_dir}/${samplename}_At_q20_flagstats.log >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt


# ===========================================================
# q1.1 Covered genome proportion (after markdup removal)
# ===========================================================
#echo "q1.1 What is the covered genome proportion of A.thaliana DNA?" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt

# Remove duplicates and sort
#samtools sort -o ${out_dir}/${samplename}_mapped_At_q20.sorted.bam ${out_dir}/${samplename}_mapped_At_q20.bam
#samtools markdup -r -@ 4 ${out_dir}/${samplename}_mapped_At_q20.sorted.bam ${out_dir}/${samplename}_mapped_At_q20.sorted.markdup.bam
#samtools sort -o ${out_dir}/${samplename}_mapped_At_q20.markdup.sorted.bam ${out_dir}/${samplename}_mapped_At_q20.sorted.markdup.bam
#samtools index ${out_dir}/${samplename}_mapped_At_q20.markdup.sorted.bam
#rm ${out_dir}/${samplename}_mapped_At_q20.markdup.bam ${out_dir}/${samplename}_mapped_At_q20.bam

#samtools coverage ${out_dir}/${samplename}_mapped_At_q20.markdup.sorted.bam > ${out_dir}/${samplename}_mapped_At_q20.markdup.sorted.bysamtoolscoverage.txt

#echo "total length of At ref:" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt
#awk '{sum+=$3;} END{print sum;}' ${out_dir}/${samplename}_mapped_At_q20.markdup.sorted.bysamtoolscoverage.txt >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt

#echo "base mapped to At ref:" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt
#awk '{sum+=$5;} END{print sum;}' ${out_dir}/${samplename}_mapped_At_q20.markdup.sorted.bysamtoolscoverage.txt >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt


# ===========================================================
# q1.2 Average depth
# ===========================================================
#echo "q1.2 What is the read depth of A.thaliana DNA?" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt
#echo "read depth to At ref:" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt
#awk '{sum+=($3*$7);} END{print sum/119667750;}' ${out_dir}/${samplename}_mapped_At_q20.markdup.sorted.bysamtoolscoverage.txt >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt


# ===========================================================
# Cleanup intermediate files
# ===========================================================
#rm ${out_dir}/${samplename}_mapped_At_q20.markdup.sorted.bysamtoolscoverage.txt


# ===========================================================
# q2: extract unmapped reads to pass to Pseudomonas
# ===========================================================
echo "q2 What is the percentage of Pseudomonas DNA after removing mapped reads to A.thaliana?" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt

# Extract unmapped reads (-bf 4)
#samtools view -bf 4 ${out_dir}/${samplename}.collapsed_At.sam > ${out_dir}/${samplename}_after_removal_mappedAt.bam

#rm ${out_dir}/${samplename}.collapsed_At.sam

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
# 1. Exclude HTF and TFA (pre-dedup) for mapped read count
# ==============================
bam_excl_predup="${base_dir}/${samplename}.mapped_to_Pseudomonas.dd.q20.excluding_HTF_TFA.bam"
samtools view -h "$in_bam" \
  | awk 'BEGIN{OFS="\t"} /^@/ {print; next} {gsub(/\r/,"",$3); ref=toupper($3)} !(ref ~ /^HTF/ || ref ~ /^TFA/)' \
  | samtools view -b -o "$bam_excl_predup"
samtools index "$bam_excl_predup"

samtools flagstat "$bam_excl_predup" > "${base_dir}/${samplename}_Ps_q20_excluding_HTF_TFA_flagstats.log"

echo "q2 What is the percentage of Pseudomonas DNA after removing mapped reads to A.thaliana?" >> "$out_answer"
echo "reads mapped to Pseudomonas (-q20) after removing:" >> "$out_answer"
grep 'mapped (' "${base_dir}/${samplename}_Ps_q20_excluding_HTF_TFA_flagstats.log" >> "$out_answer"
#checked by 
#samtools view -h PL0137.mapped_to_Pseudomonas.dd.q20.markeddup.bam \
#| awk 'BEGIN{OFS="\t"} /^@/ {print; next} {gsub(/\r/,"",$3); ref=toupper($3)} !(ref ~ /^HTF/ || ref ~ /^TFA/)' \
#| grep 'HTF'
#print only SQ header, no reads left



# ==============================
# 2. Deduplication (keep HTF/TFA)
# ==============================
# Input must be coordinate-sorted
samtools sort -o "${in_bam%.bam}.sorted.bam" "$in_bam"

# Remove duplicates
samtools markdup -r "${in_bam%.bam}.sorted.bam" "${markdup_bam%.bam}.unsorted.bam"

# Resort and index (important!)
samtools sort -o "$markdup_bam" "${markdup_bam%.bam}.unsorted.bam"

samtools index "$markdup_bam"

rm "${markdup_bam%.bam}.unsorted.bam"

# ==============================
# 3. Exclude HTF/TFA after dedup for coverage/depth
# ==============================
bam_excl_postdup="${base_dir}/${samplename}.mapped_to_Pseudomonas.dd.q20.markeddup.excluding_HTF_TFA.bam"
covfile="${base_dir}/${samplename}.mapped_to_Pseudomonas.dd.q20.markeddup.excluding_HTF_TFA.coverage.txt"
samtools view -h "$markdup_bam" \
  | awk 'BEGIN{OFS="\t"} /^@/ {print; next} {gsub(/\r/,"",$3); ref=toupper($3)} !(ref ~ /^HTF/ || ref ~ /^TFA/)' \
  | samtools view -b -o "$bam_excl_postdup"
samtools index "$bam_excl_postdup"

# build a list of contigs to keep (exclude HTF/TFA)
contig_list="${base_dir}/non_HTF_TFA_contigs.txt"
samtools idxstats "$markdup_bam" \
  | awk -F '\t' '!($1 ~ /^HTF/ || $1 ~ /^TFA/ || $1=="*"){print $1"\t0\t"$2}' \
  > "$contig_list"
echo 'q4 What is the percentage of the Pseudomonas genome that is covered by at least 1 read?' >> "$out_answer"

# compute depth only on those contigs
samtools depth -aa -b "$contig_list" "$bam_excl_postdup" > "$covfile"

echo 'pseudomonas bases in total:' >> "$out_answer"
wc -l "$covfile" | cut -d ' ' -f 1 >> "$out_answer"
echo 'bases with at least 1 read:' >> "$out_answer"
awk '($3 >= 1) {count++} END {print count}' "$covfile" >> "$out_answer"


# ==============================
# 4. NM (edit distance) distribution
# ==============================
samtools view "$bam_excl_postdup" \
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
rm -f "$bam_excl_predup" "${bam_excl_predup}.bai"
rm -f "$bam_excl_postdup" "${bam_excl_postdup}.bai"
rm -f "$covfile"

echo "All calculation done. Only dedup q20 BAM (including HTF and TFA) retained for ${samplename}."


#rm sai sam bam and unsort unq20 bams:
#keep for other refs
rm /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_realign.sai
rm /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_mapped_to_ps.sam
rm /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_mapped_to_ps.sorted.bam

rm -f ${base_dir}/${samplename}_Ps_q20_excluding_HTF_TFA_flagstats.log
rm -f ${base_dir}/${samplename}_removalAt_mapped_to_ps.bam
rm -f ${base_dir}/${samplename}.unsort_mapped_to_Pseudomonas.dd.q20.bam
rm -f ${base_dir}/${samplename}.mapped_to_Pseudomonas.dd.q20.bam


echo "Removed intermediate files for ${samplename}"
echo "All done."
#mkdir -p /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h49/At_bams_h49_dedup/
mkdir -p /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h49/Ps_bams_maptoOTU5_with_haplotype_h49_dedup/
#ln -s /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_mapped_At_q20.markdup.sorted.bam \
#/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h49/At_bams_h49_dedup/

ln -s /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}.mapped_to_Pseudomonas.dd.q20.markeddup.bam \
/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h49/Ps_bams_maptoOTU5_with_haplotype_h49_dedup/

rm /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_mapped_At_q20.sorted.markdup.bam

