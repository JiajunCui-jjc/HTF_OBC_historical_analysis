#!/usr/bin/env bash
set -euo pipefail

wd="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2026_PNAS_revision"

# Multi-FASTA: one tailocin seq per sample (query names in PAF must match these headers)
tailocin_fa="$wd/results/step2_tailocin_HTF_seqs/tailocin/final_tailocin_all53_above98percentcov.fa"

# Per-sample PAF files
paf_dir="$wd/results/step1_tailocin_mapping/pergene/step1_tailocin_betweentrpEG_pergene"

# Summary table for best HTF + length (used for output FILE NAME)
summary="$wd/results/step1_tailocin_mapping/step1_tailocin_betweentrpEG/top_matched_haplotypes_by_propcoverage_detailed.txt"

# Output directory
outdir="$wd/results/step2_tailocin_HTF_seqs/tailocin_pergene"
mkdir -p "$outdir"

command -v samtools >/dev/null 2>&1 || { echo "[ERROR] samtools not found" >&2; exit 1; }
command -v seqtk   >/dev/null 2>&1 || { echo "[ERROR] seqtk not found" >&2; exit 1; }

# Index multi-fasta for faidx
[[ -f "${tailocin_fa}.fai" ]] || samtools faidx "$tailocin_fa"

# Build lookup: sample -> top_HTFhaplotype + HTF_len
tmp_summary="$(mktemp)"
awk -F'\t' 'NR==1{next} $1!=""{print $1"\t"$2"\t"$7}' "$summary" > "$tmp_summary"
# columns: sample  top_HTFhaplotype  HTF_len

check_tsv="$outdir/check_bestHTF_TFA_fromPAF.tsv"
echo -e "sample\tbestHTF_from_summary\tHTF_len_from_summary\tbestHTF_fromPAF\tbestHTF_alnLen\tbestTFA_fromPAF\tbestTFA_alnLen\tentries_written" > "$check_tsv"

for paf in "$paf_dir"/*.paf; do
  [[ -s "$paf" ]] || continue
  sample="$(basename "$paf" .paf)"

  # best HTF + len for FILE NAME
  bestHTF_summary="$(awk -v s="$sample" '$1==s{print $2; exit}' "$tmp_summary" || true)"
  bestHTF_len="$(awk -v s="$sample" '$1==s{print $3; exit}' "$tmp_summary" || true)"
  [[ -n "${bestHTF_summary:-}" ]] || bestHTF_summary="NA"
  [[ -n "${bestHTF_len:-}" ]]     || bestHTF_len="NA"

  # output file name includes best HTF + length (like before)
  outfa="$outdir/${sample}__bestHTF_${bestHTF_summary}__HTFlen_${bestHTF_len}.fa"
  : > "$outfa"

  ###########################################################################
  # A) Best hit per BJEIHDPM_* gene (by largest alnLen, col 11)
  ###########################################################################
  awk -F'\t' '
    $6 ~ /^BJEIHDPM_/ {
      key=$6
      aln=$11+0
      if (!(key in best) || aln > best[key]) { best[key]=aln; line[key]=$0 }
    }
    END { for (k in line) print line[k] }
  ' "$paf" \
  | while IFS=$'\t' read -r qname qlen qstart qend strand tname tlen tstart tend nmatch alnlen mapq rest; do

      # PAF query coords: 0-based [qstart,qend) -> faidx: 1-based inclusive
      s=$((qstart + 1))
      e=$((qend))

      gene="${tname%%|*}"     # BJEIHDPM_024xx
      geneLen="$tlen"         # target length
      alnLen="$alnlen"        # aligned block length

      covProp="$(awk -v a="$alnLen" -v b="$geneLen" 'BEGIN{ if(b>0) printf "%.5f", a/b; else printf "NA" }')"

      # Header: NO bestHTF info, but include gene coverage
      header="${sample}|${gene}|q=${s}-${e}|paf_strand=${strand}|geneLen=${geneLen}|alnLen=${alnLen}|cov=${alnLen}/${geneLen}|covProp=${covProp}"
      region="${qname}:${s}-${e}"

      if [[ "$strand" == "-" ]]; then
        samtools faidx "$tailocin_fa" "$region" \
          | seqtk seq -r -l 0 - \
          | awk -v h="$header" 'NR==1{print ">"h; next}{print}' >> "$outfa"
      else
        samtools faidx "$tailocin_fa" "$region" \
          | seqtk seq -l 0 - \
          | awk -v h="$header" 'NR==1{print ">"h; next}{print}' >> "$outfa"
      fi
    done

  ###########################################################################
  # B) Best HTF_ (by largest alnLen) + extract sequence segment
  ###########################################################################
  bestHTF_line="$(awk -F'\t' '
    $6 ~ /^HTF_/ {
      aln=$11+0
      if (aln > best) {best=aln; line=$0}
    }
    END { if (line!="") print line }
  ' "$paf")"

  bestHTF_fromPAF="NA"; bestHTF_alnLen="0"
  if [[ -n "${bestHTF_line:-}" ]]; then
    IFS=$'\t' read -r qname qlen qstart qend strand tname tlen tstart tend nmatch alnlen mapq rest <<< "$bestHTF_line"

    s=$((qstart + 1)); e=$((qend))
    gene="${tname%%|*}"     # HTF_...
    geneLen="$tlen"
    alnLen="$alnlen"
    bestHTF_fromPAF="$gene"
    bestHTF_alnLen="$alnLen"

    covProp="$(awk -v a="$alnLen" -v b="$geneLen" 'BEGIN{ if(b>0) printf "%.5f", a/b; else printf "NA" }')"
    header="${sample}|BEST_HTF|${gene}|q=${s}-${e}|paf_strand=${strand}|geneLen=${geneLen}|alnLen=${alnLen}|cov=${alnLen}/${geneLen}|covProp=${covProp}"
    region="${qname}:${s}-${e}"

    if [[ "$strand" == "-" ]]; then
      samtools faidx "$tailocin_fa" "$region" \
        | seqtk seq -r -l 0 - \
        | awk -v h="$header" 'NR==1{print ">"h; next}{print}' >> "$outfa"
    else
      samtools faidx "$tailocin_fa" "$region" \
        | seqtk seq -l 0 - \
        | awk -v h="$header" 'NR==1{print ">"h; next}{print}' >> "$outfa"
    fi
  fi

  ###########################################################################
  # C) Best TFA_ (by largest alnLen) + extract sequence segment
  ###########################################################################
  bestTFA_line="$(awk -F'\t' '
    $6 ~ /^TFA_/ {
      aln=$11+0
      if (aln > best) {best=aln; line=$0}
    }
    END { if (line!="") print line }
  ' "$paf")"

  bestTFA_fromPAF="NA"; bestTFA_alnLen="0"
  if [[ -n "${bestTFA_line:-}" ]]; then
    IFS=$'\t' read -r qname qlen qstart qend strand tname tlen tstart tend nmatch alnlen mapq rest <<< "$bestTFA_line"

    s=$((qstart + 1)); e=$((qend))
    gene="${tname%%|*}"     # TFA_...
    geneLen="$tlen"
    alnLen="$alnlen"
    bestTFA_fromPAF="$gene"
    bestTFA_alnLen="$alnLen"

    covProp="$(awk -v a="$alnLen" -v b="$geneLen" 'BEGIN{ if(b>0) printf "%.5f", a/b; else printf "NA" }')"
    header="${sample}|BEST_TFA|${gene}|q=${s}-${e}|paf_strand=${strand}|geneLen=${geneLen}|alnLen=${alnLen}|cov=${alnLen}/${geneLen}|covProp=${covProp}"
    region="${qname}:${s}-${e}"

    if [[ "$strand" == "-" ]]; then
      samtools faidx "$tailocin_fa" "$region" \
        | seqtk seq -r -l 0 - \
        | awk -v h="$header" 'NR==1{print ">"h; next}{print}' >> "$outfa"
    else
      samtools faidx "$tailocin_fa" "$region" \
        | seqtk seq -l 0 - \
        | awk -v h="$header" 'NR==1{print ">"h; next}{print}' >> "$outfa"
    fi
  fi

  # Count written entries reliably
  entries_written="$(grep -c '^>' "$outfa" || true)"

  echo -e "${sample}\t${bestHTF_summary}\t${bestHTF_len}\t${bestHTF_fromPAF}\t${bestHTF_alnLen}\t${bestTFA_fromPAF}\t${bestTFA_alnLen}\t${entries_written}" >> "$check_tsv"

  echo "[OK] ${sample}: wrote ${entries_written} entries -> $outfa"
done

rm -f "$tmp_summary"

echo
echo "[CHECK] QC table (first 80 lines):"
column -t -s $'\t' "$check_tsv" | head -n 80

echo
echo "[CHECK] grep HTF_ hits (first few):"
grep -H $'\tHTF_' "$paf_dir"/*.paf | head
