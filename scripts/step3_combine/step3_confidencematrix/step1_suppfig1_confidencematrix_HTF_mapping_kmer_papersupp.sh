#!/bin/bash
wd="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step3_combine/suppfig1_confidence_matrix_HTF"
mkdir -p $wd
cd $wd
haplo_file="$wd/HP12_combined_lengths_sorted.txt"
output="tmp_htf_supp_with_mapping.tsv"

# Mapping table saved as plain file
mapping_list="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step3_combine/HTF_haplotypes/final_h49_m57_HTFhaplotype_bylocalassembly.txt"
#
# List of all m53 and h43 strains
samples_file="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step3_combine/samples_m53_h43.txt"

# Read sample names into array, skipping header and blank lines
mapfile -t samples < <(awk 'NR>1 && NF > 0' "$samples_file")


echo -e "Strain\tHTF_length_from_bigdataset\tHTF_by_mapping" > "$output"

for s in "${samples[@]}"; do
  # Get HTF length from big dataset (by sample name)
  htf_supp=$(awk -v id="$s" '$1 == id {print $2}' "$haplo_file")
  [[ -z "$htf_supp" ]] && htf_supp="NA"

  # Get HTF haplotype from mapping file
  htf_mapping=$(awk -v id="$s" '$1 == id { $1=""; sub(/^\t/, ""); print }' "$mapping_list")
  [[ -z "$htf_mapping" ]] && htf_mapping="NA"

  echo -e "$s\t$htf_supp\t$htf_mapping" >> "$output"
done

echo "✅ Done: Output saved to $output"


# Input files
supp_file="tmp_htf_supp_with_mapping.tsv"
HTF_FILE=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step3_combine/HTF_haplotypes/form57_h49_HTF_bykmer.txt
#cp $HTF_FILE $wd/HTF_kmer.txt
#kmer_file=$wd/HTF_kmer.txt
kmer_file=$HTF_FILE

#HTF_bykmer_modern_filtered.tsv"
#rm 64.
output_file="$wd/modern_53_h43_htf_supp_mapping_kmer.tsv"

# Extract headers and join by strain
awk 'NR==FNR{
       a[$1] = substr($0, index($0, $2)); next
     }
     FNR==1 {
       print $0 "\tHTF_by_kmer"; next
     }
     {
       kmer = ($1 in a) ? a[$1] : "NA";
       print $0 "\t" kmer
     }' "$kmer_file" "$supp_file" > "$output_file"

echo "✅ Combined file saved to: $output_file"

rm $supp_file


#final extract length:
input="$output_file"
final_output="${output_file%.tsv}_with_lengths.tsv"

awk -F'\t' 'BEGIN{OFS="\t"}
NR == 1 {
  print $0, "mapping_length", "kmer_length"
  next
}
{
  if ($1 == "64.GBR_1933b_S36") {
    # Force all columns to NA
    n = NF + 2  # original columns + 2 new ones
    printf "%s", $1
    for (i = 2; i <= n; i++) printf "\tNA"
    print ""
    next
  }
  # Extract number inside parentheses using match()
  mapping_len = "NA"
  if (match($3, /\(([0-9]+)\)/, arr)) mapping_len = arr[1]

  kmer_len = "NA"
  if (match($4, /\(([0-9]+)\)/, arr)) kmer_len = arr[1]

  print $0, mapping_len, kmer_len
}' "$input" > "$final_output"
rm $output_file

#"$final_output"
#change the line 64.GBR_1933b_S36	NA	 HTF_p21.F9 (1245)	HTF_p21.F9 (1245)	1245	1245 to 64.GBR_1933b_S36	NAforall
echo "✅ Final file with extracted lengths saved to: $final_output"
