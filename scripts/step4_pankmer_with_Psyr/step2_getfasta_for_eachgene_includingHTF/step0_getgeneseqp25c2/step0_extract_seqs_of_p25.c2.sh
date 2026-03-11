ref=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/HTFreference/OTU5_ref
wd=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2026_PNAS_revision/results/step0_refseqs/pergene
#trpE
#utg000001l_p25.C2:2694503-2695984
#trpG
#utg000001l_p25.C2:2717365-2717964
#conda activate phylogeny_snp
cd $wd

awk '
BEGIN {keep=0}

# keep GFF headers
/^#/ {print; next}

# start AFTER trpE (do not print trpE line)
/trpE/ {keep=1; next}

# stop BEFORE trpG (do not print trpG line)
/trpG/ && keep {keep=0; next}

# print only the region between them
keep {print}
' p25c2_headers_counts.gff > p25c2_tailocin_CDS.gff 







gff="p25c2_tailocin_CDS.gff"
fa="$ref/ref.withhyplotype.fasta"
contig="utg000001l_p25.C2"
out="p25c2_tailocin_CDS.fasta"

samtools faidx "$fa"

# Make CDS list but REPLACE seqid with the real contig
awk -v OFS="\t" -v contig="$contig" '
BEGIN{FS="\t"}
/^#/ {next}
$3!="CDS" {next}
{
  start=$4; end=$5; strand=$7; attr=$9;

  id="NA"
  if (match(attr, /ID=[^;]+/)) id=substr(attr, RSTART+3, RLENGTH-3)
  else if (match(attr, /locus_tag=[^;]+/)) id=substr(attr, RSTART+9, RLENGTH-9)

  print contig, start, end, strand, id
}
' "$gff" > "${out}.cds.list.tsv"

: > "$out"
while IFS=$'\t' read -r seqid start end strand id; do
  region="${seqid}:${start}-${end}"
  header="${id}|${region}|strand=${strand}"

  if [[ "$strand" == "-" ]]; then
    samtools faidx "$fa" "$region" \
      | sed "1s/^>.*/>${header}/" \
      | seqtk seq -r - >> "$out"
  else
    samtools faidx "$fa" "$region" \
      | sed "1s/^>.*/>${header}/" >> "$out"
  fi
done < "${out}.cds.list.tsv"

echo "[OK] wrote: $out"
echo "[OK] CDS count: $(grep -c '^>' "$out")"




#also HTF TFA

seqkit grep -r -p 'HTF_|TFA_' $ref/ref.withhyplotype.fasta >> $wd/p25c2_tailocin_CDS.fasta

seqtk seq -l 0 $wd/p25c2_tailocin_CDS.fasta > $wd/p25c2_tailocin_CDS.oneline.fasta

rm $wd/p25c2_tailocin_CDS.fasta


#>BJEIHDPM_02448 is the same as HTF p25.c2 so removed.

#>BJEIHDPM_02449 is the same as TFA p25.c2 so removed.


awk '
BEGIN {skip=0}
/^>/ {
  skip = ($0 ~ /BJEIHDPM_02448|BJEIHDPM_02449/)
}
!skip
' "$wd/p25c2_tailocin_CDS.oneline.fasta" \
> "$wd/p25c2_tailocin_CDS.onelinefinal.fasta"

rm -r $wd/p25c2_tailocin_CDS.oneline.fasta
