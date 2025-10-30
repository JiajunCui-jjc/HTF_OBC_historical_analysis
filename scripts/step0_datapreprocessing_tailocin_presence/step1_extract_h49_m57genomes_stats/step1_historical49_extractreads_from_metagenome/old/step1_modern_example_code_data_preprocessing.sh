#!/bin/bash -l
#$ -l tmem=3G
#$ -l h_vmem=3G
#$ -l h_rt=24:00:0
#$ -wd /SAN/ugi/plant_genom/jiajucui/
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/
#$ -V
#modern
#/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_m57/p23.B2.pair1.truncated.gz p1.G2.pair1.truncated.gzf
#variables
#except the first mDC3000
source ~/miniconda3/bin/activate phylogeny_snp
samplename=p1.G2
#and p23.B2     
wd=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_m57/
mkdir -p $wd/all_bams_m57/intermediate_tmp2
echo "${samplename}:" | tee -a ${wd}/tobamtmp2.log
refPs=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/HTFreference/OTU5_ref/Pseudomonas.OTU5_ref.fasta
rawtrimed=$wd
r1=$rawtrimed/${samplename}.pair1.truncated.gz
r2=$rawtrimed/${samplename}.pair2.truncated.gz
bwa aln -t 2 -l 1024 -f ${wd}/all_bams_m57/intermediate_tmp2/${samplename}.r1.collapsed.sai ${refPs} ${r1} &&
echo 'bwa aln r1 done' &&
bwa aln -t 2 -l 1024 -f ${wd}/all_bams_m57/intermediate_tmp2/${samplename}.r2.collapsed.sai ${refPs} ${r2} &&
echo 'bwa aln r2 done' &&
#Convert  reads into a standard alignment format (SAM)
bwa sampe -r @RG\\tID:${samplename}\\tSM:${samplename} -f ${wd}/all_bams_m57/intermediate_tmp2/${samplename}.collapsed_Ps.sam ${refPs} ${wd}/all_bams_m57/intermediate_tmp2/${samplename}.r1.collapsed.sai ${wd}/all_bams_m57/intermediate_tmp2/${samplename}.r2.collapsed.sai ${r1} ${r2} &&
echo 'bwa sampe done' &&

#convert mapped reads into bam with map quality score greater than 80 
samtools view -@ 2 -F 4  -q 20 -Sbh -o ${wd}/all_bams_m57/intermediate_tmp2/${samplename}_mapped_Ps_q20unsort.bam ${wd}/all_bams_m57/intermediate_tmp2/${samplename}.collapsed_Ps.sam

samtools sort -@ 2 -o ${wd}/all_bams_m57/intermediate_tmp2/${samplename}_mapped_Ps_q20.bam ${wd}/all_bams_m57/intermediate_tmp2/${samplename}_mapped_Ps_q20unsort.bam


