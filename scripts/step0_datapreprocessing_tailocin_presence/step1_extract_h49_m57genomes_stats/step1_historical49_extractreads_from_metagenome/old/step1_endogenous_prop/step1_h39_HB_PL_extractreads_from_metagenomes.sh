#!/bin/bash -l

#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=24:00:0
#$ -wd /SAN/ugi/plant_genom/jiajucui/
#$ -V
#$ -N s2025_${pathtmp}
#$ -t 1-40
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
wd='/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/'
pathtmp=2025_${pathtmp}
answerspath=answer2025_${pathtmp}
mkdir -p /SAN/ugi/plant_genom/jiajucui/${answerspath}/
mkdir -p  /SAN/ugi/plant_genom/jiajucui/2_trimmed_merged/${pathtmp}/
mkdir -p /SAN/ugi/plant_genom/jiajucui/3_quality_control/${pathtmp}/
samplename=$(grep -E '^(HB|PL)' /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/historical49.txt | sed -n ${i}p)
echo "${samplename}:" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt
refAt=/SAN/ugi/plant_genom/jiajucui/1_initial_data/reference_genome_At/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
refPs=/SAN/ugi/plant_genom/jiajucui/1_initial_data/reference_genome_Ps/Pseudomonas.OTU5_ref.fasta
#for PL samples: we start with raw r1 and r2
r1o="$wd/data/all_fastq_h49/all_raw_meta_fastq_h49/${samplename}.R1.fastq.gz"
r2o="$wd/data/all_fastq_h49/all_raw_meta_fastq_h49/${samplename}.R2.fastq.gz"
 

#collapse:
#trim and merge
AdapterRemoval --file1 ${r1o} --file2 ${r2o} --basename /SAN/ugi/plant_genom/jiajucui/2_trimmed_merged/${pathtmp}/${samplename} --trimns --trimqualities --collapse --gzip --threads 2 &&
echo 'trim and merge done' &&
#fastqc
fastqc -t 2 -o /SAN/ugi/plant_genom/jiajucui/3_quality_control/${pathtmp}/ /SAN/ugi/plant_genom/jiajucui/2_trimmed_merged/${pathtmp}/${samplename}.collapsed.gz /SAN/ugi/plant_genom/jiajucui/2_trimmed_merged/${pathtmp}/${samplename}.pair1.truncated.gz /SAN/ugi/plant_genom/jiajucui/2_trimmed_merged/${pathtmp}/${samplename}.pair2.truncated.gz



collapsed=/SAN/ugi/plant_genom/jiajucui/2_trimmed_merged/${pathtmp}/${samplename}.collapsed.gz
mkdir -p /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/
bwa aln -t 8 -l 1024 -f /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}.collapsed.sai ${refAt} ${collapsed}
#Convert  reads into a standard alignment format (SAM)
bwa samse -r @RG\\tID:${samplename}\\tSM:${samplename} -f /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}.collapsed_At.sam ${refAt} /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}.collapsed.sai ${collapsed}
samtools flagstat /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}.collapsed_At.sam > /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_Atsam_flagstats.log 
echo "q1 What is the percentage of A.thaliana DNA?" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt
echo "reads in total (allrawreads):" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt
grep 'in total' /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_Atsam_flagstats.log  >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt


#convert mapped reads into bam with map quality score greater than 80 
samtools view -@ 2 -F 4  -q 20 -Sbh -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_mapped_At_q20.bam /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}.collapsed_At.sam
#rm sai, sam will be removed later
rm /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}.collapsed.sai
#for q1 What is the percentage of A.thaliana DNA?

samtools flagstat /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_mapped_At_q20.bam > /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_At_q20_flagstats.log 
echo "reads mapped to At (-q20):" >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt
grep 'mapped (' /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_At_q20_flagstats.log >> /SAN/ugi/plant_genom/jiajucui/${answerspath}/answers_for_${samplename}.txt

#then rm unsortedbam 
rm /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_mapped_At_q20.bam

#q2 What is the percentage of Pseudomonas DNA? after removing mapped reads to A.thaliana

####################
##separate pseudomonas from remains DNA (removal)
samtools view -bf 4 /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}.collapsed_At.sam > /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_after_removal_mappedAt.bam &&

#then rm Atsam
rm /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}.collapsed_At.sam

##use bwa to realign the unmapped bam file
bwa aln -t 2 -l 1024 ${refPs} -b /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_after_removal_mappedAt.bam > /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_realign.sai

bwa samse -r @RG\\tID:${samplename}\\tSM:${samplename} -f /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_mapped_to_ps.sam ${refPs} /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_realign.sai /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_after_removal_mappedAt.bam
#Keep the mapped reads and create a compressed BAM file using samtools
samtools view -@ 2 -F 4 -Sbh -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_mapped_to_ps.bam /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_mapped_to_ps.sam
#Sort the BAM file by chromosome and position using samtools
samtools sort -@ 2 -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_mapped_to_ps.sorted.bam /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_mapped_to_ps.bam

samtools view -q 20 /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_mapped_to_ps.sorted.bam -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}.mapped_to_Pseudomonas.dd.q20.bam


#rm sai sam bam and unsort unq20 bams:
rm /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/${pathtmp}/${samplename}_after_removal_mappedAt.bam
rm /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_realign.sai
rm /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_mapped_to_ps.sam
rm /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_mapped_to_ps.sorted.bam

samtools flagstat /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}.mapped_to_Pseudomonas.dd.q20.bam > /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_Ps_q20_flagstats.log
#output
echo "q2 What is the percentage of Pseudomonas DNA after removing mapped reads to A.thaliana?" >> /SAN/ugi/plant_genom/jiajucui/answers${pathtmp}/answers_for_${samplename}.txt
echo "reads mapped to Pseudomonas (-q20) after removing:" >> /SAN/ugi/plant_genom/jiajucui/answers${pathtmp}/answers_for_${samplename}.txt
grep 'mapped (' /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_Ps_q20_flagstats.log >> /SAN/ugi/plant_genom/jiajucui/answers${pathtmp}/answers_for_${samplename}.txt

echo 'all calculation done'

#rm sai sam bam and unsort unq20 bams:
#keep for other refs
rm /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_realign.sai
rm /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_mapped_to_ps.sam
rm /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/${pathtmp}/${samplename}_removalAt_mapped_to_ps.sorted.bam



