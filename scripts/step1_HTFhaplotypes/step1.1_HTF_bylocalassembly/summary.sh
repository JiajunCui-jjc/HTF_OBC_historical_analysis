

cd /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.1_HTF_bylocalassembly/historical49paf/mappings
#cd ../../modern57paf/
#cat top_matched_haplotypes_by_propcoverage_detailed.txt | cut -f1,2 > modern57_HTFhaplotype_bylocalassembly.txt

cat top_matched_haplotypes_by_propcoverage_detailed.txt | cut -f1,2 > h49_HTFhaplotype_bylocalassembly.txt

#cat h40_HTFhaplotype_bylocalassembly.txt ../../modern57paf/modern57_HTFhaplotype_bylocalassembly.txt > h40_m57_HTFhaplotype_bylocalassembly.txt
tm1=h49_HTFhaplotype_bylocalassembly.txt
tm2=../../modern57paf/modern57_HTFhaplotype_bylocalassembly.txt
out=h49_m57_HTFhaplotype_bylocalassembly.txt
{
  # Print header (from the first file)
  head -n 1 "$tm1"

  # Concatenate both files, remove all header lines and blanks
  tail -n +2 "$tm1"
  tail -n +2 "$tm2"
} | grep -v '^$' > "$out"

input=h49_m57_HTFhaplotype_bylocalassembly.txt
input2=../../tmp2h49_m57_HTFhaplotype_bylocalassembly.txt

awk 'BEGIN{OFS="\t"} 
NR==1 {print $0} 
NR>1  {$2="HTF_"$2; print}' "$input" > "$input2"

output=../../final_h49_m57_HTFhaplotype_bylocalassembly.txt

awk 'BEGIN {
  OFS = "\t"
  # Define HTF to length mappings
  map["HTF_p7.G11"]    = "1830"
  map["HTF_p25.A12"]   = "1383"
  map["HTF_p25.C2"]    = "1803"
  map["HTF_p5.D5"]     = "1803"
  map["HTF_p26.D6"]    = "1803"
  map["HTF_p23.B8"]    = "1803"
  map["HTF_p21.F9"]    = "1245"
}
NR == 1 {
  print
  next
}
{
  if ($2 in map) {
    $2 = $2 " (" map[$2] ")"
  }
  print
}' "$input2" > "$output"

rm $input $input2
cp $output /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step3_combine/suppfig1_confidence_matrix_HTF/HTF_mappings.txt 
