#!/bin/bash -l

#$ -l tmem=15G
#$ -l h_vmem=15G
#$ -l h_rt=10:00:0
#$ -wd /SAN/ugi/plant_genom/jiajucui/
#$ -V
#$ -N h49_Psmapdamage 
#$ -t 1-49
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/
#$ -wd /SAN/ugi/plant_genom/jiajucui/

echo "Task id is $SGE_TASK_ID"
#conda create -n  mapDamageR4
#conda activate mapdamageR4 and conda install -c conda-forge r-base=4.1 (4.2 or more is not supported) ref:https://github.com/rstudio/gt/issues/1004
#install mapDamage
#or not install the dep 5 R libs ref:https://ginolhac.github.io/mapDamage/#a1
#then use conda to install https://anaconda.org/bioconda/mapdamage2 
#conda install -c bioconda mapdamage2

source /home/jiajucui/miniconda3/bin/activate phylogeny_snp
#variables
i=$SGE_TASK_ID
samplename=$(cat /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/historical49.txt | sed -n $i'p')
#/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step2_mapdamageplot
#(base) [jiajucui@pchuckle step2_mapdamageplot]$ ls
#stats1stbaseall49  toAt  tops
wd=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step2_mapdamageplot
refPs=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/HTFreference/OTU5_ref_noHTFhaplotypes/ref.fasta

refAt=/SAN/ugi/plant_genom/jiajucui/1_initial_data/reference_genome_At/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
mkdir -p  ${wd}/tops/
cd ${wd}/tops/
#echo "samplename: ${samplename}"  | tee -a ${wd}/md.log
#mapDamage -i /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/2023_37HB/${samplename}*.bam -r ${refAt}

#mapdamage.reader.BAMError: Incomplete readgroup found: PL0042 is missing 'LB'. Either fix BAM or use --merge-libraries
#to add the 3to5 freeq file:
#https://github.com/ginolhac/mapDamage/blob/86fbef74fdedfdea71d0c67356c72c037761fe3e/mapdamage/Rscripts/mapDamage.R
bam=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h49/Ps_bams_maptoOTU5withouthaplotype_h49_dedup/h49_bams_softlink/${samplename}.mapped_to_Pseudomonas.dd.q20.markeddup.bam
mapDamage -i $bam -r ${refPs} --merge-libraries

mapDamage -d ${samplename}* -y 0.05 --plot-only

echo 'md done'

