#!/usr/bin/env bash

# ============================================================
# Hard-coded working directory (as requested)
# ============================================================
WD="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2026_PNAS_revision/results/step2_tailocin_HTF_seqs/tailocin_pergene/final_withrealHTF"
cd "$WD"

# ============================================================
# Inputs
# ============================================================
REAL_HTF="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2026_PNAS_revision/results/step2_tailocin_HTF_seqs/htf/runinlocal_htf_real_length/results/real_trimed_HTF.fasta"

# per-isolate tailocin-per-gene FASTAs live one dir above WD
PERGENE_DIR=".."

# ============================================================
# Safety checks
# ============================================================
[[ -f "$REAL_HTF" ]] || { echo "ERROR: REAL_HTF not found: $REAL_HTF" >&2; exit 1; }
ls "$PERGENE_DIR"/p*__bestHTF_*fa >/dev/null 2>&1 || { echo "ERROR: per-isolate FASTAs not found in $PERGENE_DIR" >&2; exit 1; }

# ============================================================
# Step 1: index real trimmed HTFs -> realHTF.table.tsv
#   columns: sample  nt_len  sequence
#   sample parsed from header: >p8.B3_trimHTF1803  => p8.B3
# ============================================================
awk '
  /^>/{
    if(seq!=""){
      gsub(/^>/,"",hdr)
      split(hdr,a,"_")
      samp=a[1]
      HTF[samp]=seq
      HTFLEN[samp]=length(seq)
    }
    hdr=$0; seq=""
    next
  }
  { gsub(/[ \t\r\n]/,""); seq=seq toupper($0) }
  END{
    if(seq!=""){
      gsub(/^>/,"",hdr)
      split(hdr,a,"_")
      samp=a[1]
      HTF[samp]=seq
      HTFLEN[samp]=length(seq)
    }
    for(s in HTF)
      print s "\t" HTFLEN[s] "\t" HTF[s]
  }
' "$REAL_HTF" > realHTF.table.tsv

# ============================================================
# Step 2: replace old >...|BEST_HTF|... with real trimmed HTF
# Output: final_<sample>_HTFlen_<len>.fa in WD
# ============================================================
for fa in "$PERGENE_DIR"/p*__bestHTF_*fa; do
  samp=$(basename "$fa" | cut -d'_' -f1)
  len=$(awk -v s="$samp" '$1==s{print $2}' realHTF.table.tsv)

  if [[ -z "${len:-}" ]]; then
    echo "WARN: no real trimmed HTF found for $samp; skipping" >&2
    continue
  fi

  awk -v SAMP="$samp" -v TABLE="realHTF.table.tsv" '
    BEGIN{
      FS=OFS="\t"
      while((getline<TABLE)>0){
        hlen[$1]=$2
        hseq[$1]=$3
      }
    }

    # skip the old BEST_HTF record (header + sequence lines) in per-gene fasta
    /^>/{
      if($0 ~ "\\|BEST_HTF\\|"){
        skip=1
        next
      } else {
        skip=0
      }
    }
    skip==0 { print }

    END{
      print ">" SAMP "|BEST_HTF|length=" hlen[SAMP]
      seq=hseq[SAMP]
      for(i=1;i<=length(seq);i+=80)
        print substr(seq,i,80)
    }
  ' "$fa" > "final_${samp}_HTFlen_${len}.fa"
done

echo "Done. Outputs written to:"
echo "  $WD"




# ============================================================
# Step 3: sanity check — filename vs header vs sequence length
# Output: doublecheck_bestHTF.txt
# ============================================================

CHECK_OUT="doublecheck_bestHTF.txt"
: > "$CHECK_OUT"   # truncate if exists

echo -e "file\tfile_HTFlen\theader_HTFlen\tseq_len\tstatus" >> "$CHECK_OUT"

fail=0

for fa in final_*_HTFlen_*.fa; do
  # length from filename
  file_len=$(basename "$fa" | sed -E 's/.*_HTFlen_([0-9]+)\.fa/\1/')

  # length from BEST_HTF header
  header_len=$(awk -F'length=' '
    /^>/ && $0 ~ /\|BEST_HTF\|/ {
      split($2,a,/[^0-9]/)
      print a[1]
      exit
    }
  ' "$fa")

  # actual nucleotide sequence length
  seq_len=$(awk '
    /^>/{next}
    {gsub(/[ \t\r\n]/,""); n+=length($0)}
    END{print n}
  ' "$fa")

  status="OK"
  if [[ "$file_len" != "$header_len" ]]; then
    status="FAIL"
    fail=1
  fi

  echo -e "$(basename "$fa")\t$file_len\t$header_len\t$status" >> "$CHECK_OUT"
done

if [[ "$fail" -eq 0 ]]; then
  echo "All HTF length checks PASSED" >> "$CHECK_OUT"
else
  echo "Some HTF length checks FAILED" >> "$CHECK_OUT"
  exit 2
fi
