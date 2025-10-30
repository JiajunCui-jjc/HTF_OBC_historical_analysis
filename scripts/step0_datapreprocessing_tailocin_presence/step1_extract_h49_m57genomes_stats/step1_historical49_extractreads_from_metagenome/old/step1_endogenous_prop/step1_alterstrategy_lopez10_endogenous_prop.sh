#!/bin/bash -l

#$ -l tmem=12G
#$ -l h_vmem=12G
#$ -l h_rt=24:00:0
#$ -wd /SAN/ugi/plant_genom/jiajucui/
#$ -V
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/
#for the 10 samples from lopez et al 2022 that have some reads longer than 75bp, we would lose info by collapsing the reads (overlap reads). So we have an alternative strategy to extract as much reads as possible:
# first collapse reads to capture the short r1 r2 reads that could overlap (the molecular insert size is shorter than 150bp), second combine r1 and r2 after aln separetely and then merge as r1r2bam to capture the long 75bp reads of r1 and r2,
# then we combine the two sources to finally has the bam files
#variables
i=${1}
samplename=$(cat /SAN/ugi/plant_genom/jiajucui/1_initial_data/new_sequences/new_sampleindex.txt | sed -n $i'p' | awk '{print $2}')
echo "${samplename}:" >> /SAN/ugi/plant_genom/jiajucui/lopez10answers/answers_for_${samplename}.txt
refAt=/SAN/ugi/plant_genom/jiajucui/1_initial_data/reference_genome_At/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
refPs=/SAN/ugi/plant_genom/jiajucui/1_initial_data/reference_genome_Ps/Pseudomonas.OTU5_ref.fasta
r1o=/SAN/ugi/plant_genom/jiajucui/1_initial_data/new_sequences/${samplename}.R1.fastq.gz
r2o=/SAN/ugi/plant_genom/jiajucui/1_initial_data/new_sequences/${samplename}.R2.fastq.gz

source ~/miniconda3/bin/activate dataprepython3.10
#collapse_At.sam is from collapsed sai
#r1r2_At.sam is merged from r1.sai and r2.sai
#combined_At.sam is combined from the above two sam files


#collapse:
#trim and merge
AdapterRemoval --file1 ${r1o} --file2 ${r2o} --basename /SAN/ugi/plant_genom/jiajucui/2_trimmed_merged/lopez10/${samplename} --trimns --trimqualities --collapse --gzip --threads 2 &&
echo 'trim and merge done' &&
#fastqc
fastqc -t 2 -o /SAN/ugi/plant_genom/jiajucui/3_quality_control/secondlopez10qc/ /SAN/ugi/plant_genom/jiajucui/2_trimmed_merged/lopez10/${samplename}.collapsed.gz /SAN/ugi/plant_genom/jiajucui/2_trimmed_merged/lopez10/${samplename}.pair1.truncated.gz /SAN/ugi/plant_genom/jiajucui/2_trimmed_merged/lopez10/${samplename}.pair2.truncated.gz

r1=/SAN/ugi/plant_genom/jiajucui/2_trimmed_merged/lopez10/${samplename}.pair1.truncated.gz
r2=/SAN/ugi/plant_genom/jiajucui/2_trimmed_merged/lopez10/${samplename}.pair2.truncated.gz

#generate sai
bwa aln -t 2 -l 1024 -f /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}.collapsed.sai ${refAt} /SAN/ugi/plant_genom/jiajucui/2_trimmed_merged/lopez10/${samplename}.collapsed.gz &&
echo 'bwa aln done' &&
#generate the sam files from collapsed sai
bwa samse -f /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}.collapsed_At.sam ${refAt} /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}.collapsed.sai /SAN/ugi/plant_genom/jiajucui/2_trimmed_merged/lopez10/${samplename}.collapsed.gz &&
echo 'collapsed bwa samse done' &&
#create a compressed BAM file using samtools
samtools view -b -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}.collapsed_At.bam /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}.collapsed_At.sam
#Sort the BAM file by chromosome and position using samtools
samtools sort -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}.collapsed_At.sorted.bam /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}.collapsed_At.bam
#for restore
rm /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}.collapsed_At.sam &&
echo 'rm collpased.sam' &&

#r1 and r2 first aln separetely and then merge as r1r2bam
#trim files were given by adapterremoval previously names as pair1.truncated.gz and 2. When try single-end mode, adapterremoval give 1.CAN_1994a_S53.truncated.gz for r1 with size 298mb. with paired-end mode, give pair1.truncated.gz wtih size 238mb. I guess we should use paired mode as we need to merge them later.
#ref:https://adapterremoval.readthedocs.io/en/stable/examples.html#trimming-paired-end-reads
#This command generates the files output_paired.pair1.truncated and output_paired.pair2.truncated, which contain trimmed pairs of reads which were not collapsed, output_paired.singleton.truncated containing reads where one mate was discarded, output_paired.collapsed containing merged reads, and output_paired.collapsed.truncated containing merged reads that have been trimmed due to the --trimns or --trimqualities options. 

r1=/SAN/ugi/plant_genom/jiajucui/2_trimmed_merged/lopez10/${samplename}.pair1.truncated.gz
r2=/SAN/ugi/plant_genom/jiajucui/2_trimmed_merged/lopez10/${samplename}.pair2.truncated.gz

bwa aln -t 2 -l 1024 -f /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}.r1.truncated.sai ${refAt} ${r1} &&
echo 'bwa aln r1 done' &&
bwa aln -t 2 -l 1024 -f /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}.r2.truncated.sai ${refAt} ${r2} &&
echo 'bwa aln r2 done' &&
# Convert  reads into a standard alignment format (SAM)
bwa sampe -f /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}.r1r2_At.sam ${refAt} /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}.r1.truncated.sai /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}.r2.truncated.sai ${r1} ${r2} &&
echo 'bwa sampe done' &&
#create a compressed BAM file using samtools
samtools view -b -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}.r1r2_At.bam /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}.r1r2_At.sam
#sort the BAM file by chromosome and position using samtools
samtools sort -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}.r1r2_At.sorted.bam /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}.r1r2_At.bam &&
echo 'samtobam and sort done for two' &&

#for restore
rm /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}.r1r2_At.sam &&
echo 'rm r1r2.sam' &&


#generate the index

samtools index /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}.r1r2_At.sorted.bam

samtools index /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}.collapsed_At.sorted.bam

echo 'part1 done'

#use samtools merge to combine sam from collapsed sai and sam merged after algin to both fastq reads
#-n means merge by read name
samtools merge -n /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}.combined_At.bam /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}.r1r2_At.sorted.bam /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}.collapsed_At.sorted.bam
#need to be sorted again
samtools sort -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}.combined_At.sorted.bam /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}.combined_At.bam &&
echo 'samtools cat combine done' &&

samtools flagstat /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}.combined_At.sorted.bam > /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}_Atsam_flagstats.log 
#convert mapped reads into bam with map quality score greater than 80 
#-b means output format is bam, -h means write the header
samtools view -@ 2 -F 4 -q 20 -bh -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}_mapped_At.q20.bam /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}.combined_At.sorted.bam


#for q1 What is the percentage of A.thaliana DNA?
echo "q1 What is the percentage of A.thaliana DNA?" >> /SAN/ugi/plant_genom/jiajucui/lopez10answers/answers_for_${samplename}.txt
echo "reads in total (allrawreads):" >> /SAN/ugi/plant_genom/jiajucui/lopez10answers/answers_for_${samplename}.txt
grep 'in total' /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}_Atsam_flagstats.log  >> /SAN/ugi/plant_genom/jiajucui/lopez10answers/answers_for_${samplename}.txt

samtools flagstat /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}_mapped_At.q20.bam > /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}_At_flagstats.q20.log 
echo "reads mapped to At ref:" >> /SAN/ugi/plant_genom/jiajucui/lopez10answers/answers_for_${samplename}.txt
grep 'mapped (' /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}_At_flagstats.q20.log >> /SAN/ugi/plant_genom/jiajucui/lopez10answers/answers_for_${samplename}.txt



####################
#q2 What is the percentage of Pseudomonas DNA? after removing mapped reads to A.thaliana
###method 1
#-b means the output format is bam, -f means only keep the unmapped reads
samtools view -bf 4 /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}.combined_At.sorted.bam > /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}_after_removal_mappedAt.bam

#sort it 

samtools sort -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}_after_removal_mappedAt.sorted.bam /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}_after_removal_mappedAt.bam
##use bwa to realign the unmapped bam file
bwa aln -t 2 -l 1024 ${refPs} -b /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}_after_removal_mappedAt.sorted.bam > /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/lopez10/${samplename}_removalAt_realign.sorted.sai

##sai to sam
bwa samse -f /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/lopez10/${samplename}_removalAt_mapped_to_ps.sam ${refPs} /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/lopez10/${samplename}_removalAt_realign.sorted.sai /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/lopez10/${samplename}_after_removal_mappedAt.sorted.bam

#Keep the mapped reads and create a compressed BAM file using samtools -S means the input format is sam and q20
samtools view -@ 2 -F 4 -q 20 -Sbh -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/lopez10/${samplename}_removalAt_mapped_to_ps.q20.bam /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/lopez10/${samplename}_removalAt_mapped_to_ps.sam

#Sort the BAM file by chromosome and position using samtools
samtools sort -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/lopez10/${samplename}_removalAt_mapped_to_ps.q20.sorted.bam /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/lopez10/${samplename}_removalAt_mapped_to_ps.q20.bam

#output
echo "q2 What is the percentage of Pseudomonas DNA after removing mapped reads to A.thaliana?" >> /SAN/ugi/plant_genom/jiajucui/lopez10answers/answers_for_${samplename}.txt

echo "reads mapped to Pseudomonas after removing:" >> /SAN/ugi/plant_genom/jiajucui/lopez10answers/answers_for_${samplename}.txt
grep 'mapped (' /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/lopez10/${samplename}_Ps_flagstats.q20.log >> /SAN/ugi/plant_genom/jiajucui/lopez10answers/answers_for_${samplename}.txt

echo 'all calculation done'

