#!/usr/bin/env bash

# ============================================================
# Usage:
#   bash get_real_trimHTF.sh htf_plusminus200.oneline.fa emboss_transeq.out
#
# Outputs:
#   real_trimed_HTF.summary.tsv
#   real_trimed_HTF.fasta
# ============================================================
cd /Users/cuijiajun/Desktop/others/tmphernan/2026_PNAS_revision/results/htftailocin_betweentrpEandtrpG/htf



NUC_FA="htf_plusminus200.oneline.fa"
AA_FA="emboss_transeq-I20260202-223248-0720-72893955-p1m.out"
mkdir -p ./results
OUT_TSV="./results/real_trimed_HTF.summary.tsv"
OUT_FA="./results/real_trimed_HTF.fasta"

awk -v OUT_TSV="$OUT_TSV" -v OUT_FA="$OUT_FA" '


function revcomp(s,    i,c,rc) {
  rc=""
  for (i=length(s); i>=1; i--) {
    c=toupper(substr(s,i,1))
    if (c=="A") rc=rc "T"
    else if (c=="T") rc=rc "A"
    else if (c=="G") rc=rc "C"
    else if (c=="C") rc=rc "G"
    else rc=rc "N"
  }
  return rc
}


function best_orf(seq,    n,i,j,c,cur,best,bestpos) {
  best=0; bestpos=0
  n=length(seq)

  for (i=1; i<=n-1; i++) {
    # ORF must start with "MD" (M followed by D)
    if (substr(seq,i,2)=="MD") {
      cur=0
      for (j=i; j<=n; j++) {
        c=substr(seq,j,1)
        if (c=="*") break
        cur++
      }
      if (cur>best) { best=cur; bestpos=i }
    }
  }

  _best_len=best
  _best_pos=bestpos
}

function sample_from_header(h,    x) {
  x=h
  sub(/^>/,"",x)
  sub(/_.*/,"",x)
  return x
}

function frame_from_header(h,    x,seg) {
  x=h
  sub(/^>/,"",x)
  if (match(x, /_-?[123]$/)) {
    seg=substr(x, RSTART, RLENGTH)   # like "_1" or "_-2"
    return seg
  }
  return ""
}

# Parse region start/end and strand from nuc header.
# Supports:
#   ..._<start>-<end>_strand=+...
# and fallback:
#   ...:<start>-<end>|strand=+...
function parse_region_and_strand(h,    seg,tmp,a,b) {
  reg_start=""; reg_end=""; strand=""

  if (match(h, /_[0-9]+-[0-9]+_strand=[+-]/)) {
    seg=substr(h, RSTART, RLENGTH)           # "_557550-559180_strand=+"
    sub(/^_/,"",seg)                          # "557550-559180_strand=+"
    split(seg, tmp, "_strand=")               # tmp[1]="557550-559180" tmp[2]="+"
    strand=tmp[2]
    split(tmp[1], a, "-")
    reg_start=a[1]+0
    reg_end=a[2]+0
    return
  }

  if (match(h, /:[0-9]+-[0-9]+\|strand=[+-]/)) {
    seg=substr(h, RSTART, RLENGTH)            # ":557550-559180|strand=+"
    sub(/^:/,"",seg)                           # "557550-559180|strand=+"
    split(seg, tmp, "\\|strand=")              # tmp[1]="557550-559180" tmp[2]="+"
    strand=tmp[2]
    split(tmp[1], b, "-")
    reg_start=b[1]+0
    reg_end=b[2]+0
    return
  }
}

# -------------------------
# PASS 1: Read nuc fasta
# -------------------------
FNR==NR {
  if ($0 ~ /^>/) {
    nuc_sample=sample_from_header($0)
    parse_region_and_strand($0)

    NUC_HDR[nuc_sample]=$0
    NUC_REGSTART[nuc_sample]=reg_start
    NUC_REGEND[nuc_sample]=reg_end
    NUC_STRAND[nuc_sample]=strand
    NUC_SEQ[nuc_sample]=""
    next
  } else {
    gsub(/[ \t\r\n]/,"")
    NUC_SEQ[nuc_sample]=NUC_SEQ[nuc_sample] toupper($0)
    next
  }
}

# -------------------------
# PASS 2: Read AA fasta
# -------------------------
function flush_aa(    samp,fr,bestlen,bestpos) {
  if (aa_id=="") return
  samp=sample_from_header(aa_id)
  fr=frame_from_header(aa_id)

  best_orf(aa_seq)
  bestlen=_best_len
  bestpos=_best_pos

  if (bestlen>0 && (!(samp in BEST_AA_LEN) || bestlen>BEST_AA_LEN[samp])) {
    BEST_AA_LEN[samp]=bestlen
    BEST_AA_POS[samp]=bestpos
    BEST_FRAME[samp]=fr
  }

  aa_id=""
  aa_seq=""
}

FNR!=NR {
  if ($0 ~ /^>/) {
    flush_aa()
    aa_id=$0
    aa_seq=""
    next
  } else {
    gsub(/[ \t\r\n]/,"")
    aa_seq=aa_seq toupper($0)
    next
  }
}

END {
  flush_aa()

  print "sample\tbest_frame\tbest_aa_len\tbest_aa_start\ttrim_nt_len_including_stop\tgenome_nt_start\tgenome_nt_end\tregion_start\tregion_end\tregion_len" > OUT_TSV

  # truncate fasta
  close(OUT_FA)

  for (samp in BEST_AA_LEN) {
    aa_len=BEST_AA_LEN[samp]
    aa_pos=BEST_AA_POS[samp]
    fr=BEST_FRAME[samp]

    nuc=NUC_SEQ[samp]
    if (nuc=="") continue

    reg_start=NUC_REGSTART[samp]
    reg_end=NUC_REGEND[samp]
    L=length(nuc)

    fr2=fr
    sub(/^_/,"",fr2)
    frame=fr2+0
    frame_abs=(frame<0)?-frame:frame

    trim_nt_len=(aa_len+1)*3

    if (frame>0) {
      offset=(frame_abs-1)
      rel_start = offset + (aa_pos-1)*3 + 1
      rel_end   = rel_start + trim_nt_len - 1

      g_start = reg_start + (rel_start-1)
      g_end   = reg_start + (rel_end-1)

      trim = substr(nuc, rel_start, trim_nt_len)
    } else {
      offset=(frame_abs-1)

      rel_startcodon_end   = L - 3*(aa_pos-1) - offset
      rel_startcodon_start = rel_startcodon_end - 2

      stop_aa = aa_pos + aa_len
      rel_stopcodon_end   = L - 3*(stop_aa-1) - offset
      rel_stopcodon_start = rel_stopcodon_end - 2

      rel_low  = rel_stopcodon_start
      rel_high = rel_startcodon_end

      g_start = reg_start + (rel_low-1)
      g_end   = reg_start + (rel_high-1)

      frag = substr(nuc, rel_low, trim_nt_len)
      trim = revcomp(frag)
    }

    # Write table even if trim length mismatched (debug-friendly)
    print samp, fr, aa_len, aa_pos, trim_nt_len, g_start, g_end, reg_start, reg_end, L >> OUT_TSV

    # Only write fasta if exact length matches
    if (length(trim)!=trim_nt_len) continue

    hdr=">" samp "_trimHTF" trim_nt_len
    print hdr >> OUT_FA
    for (i=1; i<=length(trim); i+=80) print substr(trim,i,80) >> OUT_FA
  }

  close(OUT_TSV)
  close(OUT_FA)
}
' "$NUC_FA" "$AA_FA"

awk 'BEGIN{OFS="\t"} NR==1{print "sample","fr","aa_len","trim_nt_len"; next}
     {print $1,$2,$3,$5}' ./results/real_trimed_HTF.summary.tsv > ./results/bestHTF.tsv

echo "Wrote:"
echo "  - $OUT_TSV"
echo "  - $OUT_FA"
