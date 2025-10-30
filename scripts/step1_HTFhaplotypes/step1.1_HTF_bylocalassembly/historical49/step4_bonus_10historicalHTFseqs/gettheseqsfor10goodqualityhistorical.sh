#!/usr/bin/env bash
set -euo pipefail
wd=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.1_HTF_bylocalassembly/historical46paf/mappings/
# Input file (the table you pasted with header)
INPUT="$wd/top_matched_haplotypes_by_propcoverage_detailed.txt"
# Output file (sample + best HTF sorted by proportion)
OUTPUT="$wd/getseqs_goodquality_10samples_65percent/results/forseqs_goodqualityselected_samples.tsv"

# 1. Keep header
printf "sample\ttop_HTFhaplotype\tHTF_prop\n" > "$OUTPUT"

# 2. Sort by HTF_prop column (3rd field), numeric descending
#    Keep only rows with prop > 0.65
# i checked prop > 0.5 but those samples cant be recovered
awk 'NR>1 && $3 > 0.65 {print $1, $2, $3}' OFS="\t" "$INPUT" \
  | sort -k3,3nr >> "$OUTPUT"

echo "Done. Sorted list (HTF_prop > 0.65) written to $OUTPUT"


# === Paths ===
wd=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.1_HTF_bylocalassembly/historical46paf/mappings
assemblies_tar=$wd/../assemblies/assemblies.tar.gz
refs=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/HTFreference/forkmer_step0_7refs/HTF_7refs
selected=$wd/getseqs_goodquality_10samples_65percent/results/forseqs_goodqualityselected_samples.tsv
outdir=$wd/getseqs_goodquality_10samples_65percent/results/extracted_seqs
assemblies_tmp=$wd/getseqs_goodquality_10samples_65percent/results/assemblies_tmp

mkdir -p "$outdir" "$assemblies_tmp"

# === Extract contigs.fasta for selected samples from tar.gz ===
tail -n +2 "$selected" | cut -f1 | while read -r sample; do
  echo "Extracting contigs for $sample ..."
  tar -xzf "$assemblies_tar" -C "$assemblies_tmp" "assemblies/${sample}/contigs.fasta" || {
    echo "No contigs found for $sample in tarball"
    continue
  }
done

# === Loop over selected samples ===
tail -n +2 "$selected" | while IFS=$'\t' read -r sample htf prop; do
  echo "Processing $sample ($htf, prop=$prop)..."

  paf="$wd/${sample}_mapped.paf"
  fasta="$assemblies_tmp/assemblies/$sample/contigs.fasta"

  if [[ ! -s "$paf" ]]; then
    echo "Missing PAF for $sample" >&2
    continue
  fi
  if [[ ! -s "$fasta" ]]; then
    echo "Missing contigs for $sample" >&2
    continue
  fi


done

# === Merge into one FASTA ===

# === (Optional) Run MSA ===
# mafft --auto "$outdir/selected10_forMSA.fa" > "$outdir/selected10_forMSA.msa.fa"
# echo "MSA ready: $outdir/selected10_forMSA.msa.fa"
