#!/bin/bash -l
#$ -l tmem=80G
#$ -l h_vmem=80G
#$ -l h_rt=2:00:0
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -V
#$ -N h49_localassembly
#$ -t 1-49
#$ -hold_jid PLHB39_reads,lopez10_reads
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/

# ===========================
# Activate environment
# ===========================
source /home/jiajucui/miniconda3/bin/activate phylogeny_snp

# ===========================
# Define directories
# ===========================
bam_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/all_fastq_h49/Ps_bams_maptoOTU5_with_haplotype_h49_dedup/h49_bams_softlink"
output_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.1_HTF_bylocalassembly/historical49paf"

tailocin_region="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/HTFreference/regions.txt"
reference_genome="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/HTFreference/OTU5_ref/Pseudomonas.plate25.C2.pilon.contigs_renamed.with_Tail_Fiber_Haps.fasta"
tools="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/shfortailocin/tools"

tmp_dir="${output_dir}/tmp"
readids_dir="${output_dir}/readids"
trimmed_fastq_dir="${output_dir}/trimmed_tailocin_fastq"
assembly_dir="${output_dir}/assemblies"
mapping_dir="${output_dir}/mappings"

mkdir -p "$tmp_dir" "$readids_dir" "$trimmed_fastq_dir" "$assembly_dir" "$mapping_dir" 

nocontigs_file="${mapping_dir}/nocontigs.txt"
contigmapping_file="${mapping_dir}/contigmapping.txt"
> "$nocontigs_file"
> "$contigmapping_file"

tailocin_fasta="${output_dir}/tailocin_regiondoublecheck.fa"
> "$tailocin_fasta"

# ===========================
# Extract tailocin reference regions
# ===========================
while read -r region name; do
  if [ "$name" == "tailocin" ]; then
    samtools faidx "$reference_genome" "$region" | sed "1s/.*/>$name/" >> "$tailocin_fasta"
  else
    samtools faidx "$reference_genome" "$region" | seqtk seq -r | sed "1s/.*/>$name/" >> "$tailocin_fasta"
  fi
done < "$tailocin_region"

# ===========================
# Select BAM for this array job
# ===========================
bam_list=($(ls "$bam_dir"/*.bam))
bam="${bam_list[$((SGE_TASK_ID-1))]}"
samplename=$(basename "$bam" .mapped_to_Pseudomonas.dd.q20.markeddup.bam)

echo "Processing sample: $samplename"

# ===========================
# SPAdes local assembly
# ===========================
subset_fastq="${trimmed_fastq_dir}/${samplename}_subset.fastq.gz"
subset_fastq1="${trimmed_fastq_dir}/${samplename}_subsetR1.fastq.gz"
subset_fastq2="${trimmed_fastq_dir}/${samplename}_subsetR2.fastq.gz"
assembly_out="${assembly_dir}/${samplename}"
mkdir -p "$assembly_out"

if [[ "$samplename" == *.* ]]; then
  spades.py --merge "$subset_fastq" -1 "$subset_fastq1" -2 "$subset_fastq2" \
            --careful -k 21,33 -o "$assembly_out"
else
  spades.py -s "$subset_fastq" --careful -k 21,33 -o "$assembly_out"
fi

# ===========================
# Map assembled contigs
# ===========================
contig_file="$assembly_out/contigs.fasta"
paf_file="${mapping_dir}/${samplename}_mapped.paf"

if [[ ! -f "$contig_file" ]]; then
  echo "$samplename" >> "$nocontigs_file"
  echo "No contigs found for $samplename"
  exit 0
fi

$tools/minimap2/minimap2 -cx asm5 "$contig_file" "$tailocin_fasta" > "$paf_file"

if [[ $(wc -l < "$paf_file") -eq 0 ]]; then
  echo "$samplename nomappedcontig" >> "$nocontigs_file"
  echo "No mapped contigs found for $samplename"
  exit 0
fi

echo "Completed sample: $samplename"
