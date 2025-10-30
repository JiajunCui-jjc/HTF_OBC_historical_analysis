#!/bin/bash
#$ -l tmem=100G
#$ -l h_vmem=100G
#$ -l h_rt=2:30:0
#$ -S /bin/bash
#$ -N h49_localassembly
#$ -hold_jid PLHB39_reads,lopez10_reads
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/

# Activate the conda environment
source /home/jiajucui/miniconda3/bin/activate phylogeny_snp

# Define the base directories
bam_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h49/Ps_bams_maptoOTU5_with_haplotype_h49_dedup/h49_bams_softlink"
output_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.1_HTF_bylocalassembly/historical49paf"
formatted_regions="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/HTFreference/formatted_regions.txt"
mkdir -p $output_dir
tmp_dir="${output_dir}/tmp"
readids_dir="${output_dir}/readids"

trimmed_fastq_dir="${output_dir}/trimmed_tailocin_fastq/"
assembly_dir="${output_dir}/assemblies"
mapping_dir="${output_dir}/mappings"
nocontigs_file="${mapping_dir}/nocontigs.txt"
contigmapping_file="${mapping_dir}/contigmapping.txt"

rm -r $assembly_dir $mapping_dir 

mkdir -p $assembly_dir  $mapping_dir 

tailocin_fasta="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/HTFreference/tailocin_regiondoublecheck.fa"
reference_genome="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/HTFreference/OTU5_ref/Pseudomonas.plate25.C2.pilon.contigs_renamed.with_Tail_Fiber_Haps.fasta"
tailocin_region="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/HTFreference/regions.txt"

#in the formated_region, all the TFA and HTF are reversed, make sure all of them are strand + like tailocin first
#I double checked with omega and found the *p25.C2 are the reference indeed but in reverse strand, so do the reverse and dont need to include the extracted ref chunk for TFA and HTF. only 7 and 7
# Read the tailocin region and process it
>"$tailocin_fasta"
while read -r region name; do
  if [ "$name" == "tailocin" ]; then
    samtools faidx "$reference_genome" "$region" | sed "1s/.*/>$name/" >> "$tailocin_fasta"
  else
    samtools faidx "$reference_genome" "$region" | seqtk seq -r | sed "1s/.*/>$name/" >> "$tailocin_fasta"
  fi
done < "$tailocin_region"
# Initialize or clear the nocontigs_file and contigmapping_file
> "$nocontigs_file"
> "$contigmapping_file" 

# Tools path
tools=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/shfortailocin/tools





#SPAdes
for bam in "$bam_dir"/*.bam; do
  samplename=$(basename "$bam" .mapped_to_Pseudomonas.dd.q20.markeddup.bam)
  subset_fastq="${trimmed_fastq_dir}/${samplename}_subset.fastq.gz"
  # for 143 samples_cat_merged use --12 
  # --12 <file_name> File with interlaced forward and reverse paired-end reads.
  # --merged <file_name> File with merged paired reads. If the properties of the library permit, overlapping paired-end reads can be merged using special software.
  # Non-empty files with (remaining) unmerged left/right reads (separate or interlaced) must be provided for the same library for SPAdes to correctly detect the original read length.
  # but no additional files, --12 is ok
  # --isolate - isolate (standard) bacterial data;
  # --meta The --meta mode is designed for metagenomic data, which typically involves a heterogeneous mixture of reads. also good for dealing with a mixture of read qualities
  # --only-error-correction Performs read error correction only.
  # --only-assembler Runs assembly module only. if you have high quality reads, but here we have a mix of qualities, better run error correction
  # checked --meta gave more info in low quality reads like in120, it was nothing but here 30 contigs 
  # Determine the correct raw FASTQ file name pattern
  #--careful Reduces mismatches / indels, slightly slower, more fragmented
  #but the PL0042 process is frozen even using -t 2 and require 50Gb, so remove --meta for it, and it works
  #for 109. 120. and so on in total 10 samples who have R1 and R2 raw fastq:
  if [[ "$samplename" == *.* ]]; then
    subset_fastq1="${trimmed_fastq_dir}/${samplename}_subsetR1.fastq.gz"
    subset_fastq2="${trimmed_fastq_dir}/${samplename}_subsetR2.fastq.gz"
    
    spades.py --merge "$subset_fastq" -1 $subset_fastq1 -2 $subset_fastq2  --careful -o "$assembly_dir/$samplename" -k 21,33
  #for HB and PL samples who has only one truncated and merged raw fastq...
  else
  # single end has not meta just use multicell/isolate as default
    spades.py -s "$subset_fastq"  -o "$assembly_dir/$samplename"  --careful -k 21,33
  fi
  # Step 2: Check the number of contigs in each assembly
  contig_file="$assembly_dir/$samplename/contigs.fasta"
#  contig_count=$(grep -c "^>" "$contig_file")
#  echo "$samplename: $contig_count contigs" >> "$contig_stats_dir/contig_counts.txt"
  # Check if the contig file exists
  if [[ ! -f "$contig_file" ]]; then
    echo "$samplename" >> "$nocontigs_file"
    echo "Contig file $contig_file does not exist. Sample name $samplename added to $nocontigs_file."
    continue
  fi



  # Define your output files
  paf_file="${mapping_dir}/${samplename}_mapped.paf"
  #  Run minimap2 to get the mapping
  #sergio's idea: map the reference short chunk of hyplotypes to the long contig, in principle works for both historical and modern samples.
  #$tools/minimap2/minimap2 -cx asm5 "$reference_genome" "$contig_file" > "$paf_file"
  $tools/minimap2/minimap2 -cx asm5 "$contig_file" "$tailocin_fasta" > "$paf_file"


  if [[ $(less "$paf_file" | wc -l) -eq 0 ]]; then
    echo "$samplename nomappedcontig" >> "$nocontigs_file"
    echo "Contig file $contig_file does not exist. Sample name $samplename added to $nocontigs_file."
    continue
  fi
  
done


#then summary like in modern57
