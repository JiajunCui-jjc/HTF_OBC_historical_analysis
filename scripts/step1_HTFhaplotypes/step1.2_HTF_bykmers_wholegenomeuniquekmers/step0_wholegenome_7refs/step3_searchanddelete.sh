
wd=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/HTFreference/forkmer_step0_7refs/excludeselfHTF/step3_wg.fa.excludeselfHTF/
#
#cleaned_noHTF/                     step1_oneline.sh
#onelinefasta/                      step2_manuallyserachanddelete.sh
#p13.C7.noselfHTF.fasta             step2_manuallyserachanddelete.txt
for f in $wd/onelinefasta/*.fasta; do
  out="${f%.fasta}.oneline.fasta"
  echo "[*] Converting $f to $out"
  awk '/^>/ {if (seq) print seq; print; seq=""; next} {seq=seq $0} END {print seq}' "$f" > "$out"
done


mkdir -p $wd/cleaned_noHTF

for f in $wd/onelinefasta/*.fasta; do
    sample=$(basename "$f" .oneline.fasta)   # e.g., p13.C7
    htf="$wd/../HTFseqs_samestrandasinwg/sHTF_${sample}.one.fa"
    out="$wd/cleaned_noHTF/${sample}.noselfHTF.fasta"

    echo "[*] Processing $sample ..."
    cat $htf

    awk -v htf="$htf" 'BEGIN {
        while ((getline l < htf) > 0) {
            if (l ~ /^>/) continue    # skip FASTA header
            seq = seq l               # concatenate sequence lines
        }
    }
    /^>/ {print; next}
    {gsub(seq, "", $0); print}
    ' "$f" > "$out"
done
