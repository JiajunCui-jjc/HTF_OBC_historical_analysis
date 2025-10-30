#!/bin/bash -l
#$ -l tmem=12G                # memory per core
#$ -l h_vmem=12G
#$ -l h_rt=12:00:00           # runtime
#$ -cwd
#$ -j y
#$ -V
#$ -S /bin/bash
#$ -N PLHB39_reads
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/
# Define the base directories
source /home/jiajucui/miniconda3/bin/activate phylogeny_snp

# Define the base directories
bam_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h49/Ps_bams_maptoOTU5_with_haplotype_h49_dedup/h10_lopez"
fastq_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h49/all_raw_meta_fastq_h49/"
output_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.1_HTF_bylocalassembly/historical49paf"
tmp_dir="${output_dir}/tmp"
readids_dir="${output_dir}/readids"
trimmed_fastq_dir="${output_dir}/trimmed_tailocin_fastq"

# Create output directories
mkdir -p $output_dir
mkdir -p "$tmp_dir"
mkdir -p "$readids_dir"
mkdir -p "$trimmed_fastq_dir"
# Define regions and haplotypes with corrected lengths and positions
#TFA_p23.B8:1-513
#TFA_p26.D6:1-513
#TFA_p21.F9:1-498
#TFA_p25.C2:1-513
#TFA_p5.D5:1-513
#TFA_p25.A12:1-546
#TFA_p7.G11:1-546
#HTF_p23.B8:1-1803
#HTF_p26.D6:1-1803
#HTF_p21.F9:1-1245
#HTF_p25.C2:1-1803
#HTF_p5.D5:1-1803
#HTF_p25.A12:1-1383
#HTF_p7.G11:1-1830
#already double checked samtools view bam region, every region works

#but since SPAdes --isolate requires uniform coverage, the length of regions should be as similar as possible so we use tailocin genes:
formatted_regions="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/HTFreference/formatted_regions.txt"
regions=()
while IFS= read -r line; do
  regions+=("$line")
done < "$formatted_regions"



for bam in "$bam_dir"/*.bam; do
  # Extract sample name from BAM file name
   samplename=$(basename "$bam" .mapped_to_Pseudomonas.dd.q20.rmdup.bam)
 
  # Initialize read IDs file for the sample
  read_ids="${tmp_dir}/${samplename}_read_ids.txt"
  > "$read_ids"
  samtools index $bam 
  # Check references in the BAM file
  samtools idxstats "$bam" | cut -f 1 > "${tmp_dir}/${samplename}_references.txt"
  
  # Extract read IDs for each region and haplotype
  for region in "${regions[@]}"; do
    ref=$(echo "$region" | cut -d: -f1)
    if grep -q "$ref" "${tmp_dir}/${samplename}_references.txt"; then
      samtools view "$bam" "$region" | awk '{print $1}' >> "$read_ids"
    else
      echo "Reference $ref not found in $bam"
    fi
  done
  
  sort "$read_ids" | uniq > "${readids_dir}/${samplename}_all_read_ids.txt"
  
  # Subset the reads from the raw FASTQ files
  # Find the raw FASTQ file using pattern matching
# Determine the correct raw FASTQ file name
  raw_fastq="${fastq_dir}/${samplename}.collapsed.gz"

  # Find the raw FASTQ file using the determined pattern
  subset_fastq="${trimmed_fastq_dir}/${samplename}_subset.fastq.gz"
  
  
  # Run seqtk subseq and handle errors
  /SAN/ugi/plant_genom/jiajucui/tools/seqtk/seqtk subseq "$raw_fastq" "${readids_dir}/${samplename}_all_read_ids.txt" | gzip > "$subset_fastq"
  echo "finish $samplename" 
  
  # Remove intermediate files
  rm -f "$read_ids" "${tmp_dir}/${samplename}_references.txt"
done

echo "Processing complete. Check the directories for results."
