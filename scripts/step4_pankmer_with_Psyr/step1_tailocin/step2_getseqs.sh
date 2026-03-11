#!/bin/bash -l

#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=10:00:0
#$ -wd /SAN/ugi/plant_genom/jiajucui/
#$ -V
#$ -N getseq
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/

source ~/miniconda3/bin/activate phylogeny_snp



###############################################################################
# Extract:
#  (1) Tailocin region sequences for query: utg000001l_p25.C2:2695984-2717365
#  (2) HTF best-haplotype sequences +/-200bp based on TOPTAB, then PAF
#
# Rules:
#  - If multiple hits: choose the one with largest query-aligned span (qend-qstart)
#  - Coverage proportion = (qend-qstart)/qlen
#  - Normalize orientation: if strand is '-', reverse-complement extracted sequence
#
# Outputs:
#  - tailocin/tailocin_regions.fa
#  - htf/htf_plusminus200.fa
#  - summaries/tailocin_besthit_summary.tsv
#  - summaries/htf_besthit_summary.tsv
###############################################################################

# -----------------------------
# USER CONFIG (EDIT THESE)
# -----------------------------
PAF_DIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2026_PNAS_revision/results/step1_tailocin_mapping/step1_tailocin_betweentrpEG"

# Directory containing assemblies for each sample: p12.A11.fasta(.gz) etc.
ASM_DIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2026_PNAS_revision/data/modern53_fasta"

OUTDIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2026_PNAS_revision/results/step2_tailocin_HTF_seqs"

TOPTAB="${PAF_DIR}/top_matched_haplotypes_by_propcoverage_detailed.txt"

TAILOCIN_QNAME="utg000001l_p25.C2:2695984-2717365"
FLANK=200

# -----------------------------
# CHECKS
# -----------------------------
mkdir -p "${OUTDIR}"/{tmp,tailocin,htf,summaries}

for exe in awk sed grep; do
  command -v "$exe" >/dev/null 2>&1 || { echo "[ERROR] missing executable: $exe" >&2; exit 1; }
done
command -v samtools >/dev/null 2>&1 || { echo "[ERROR] samtools not found (needed for faidx)." >&2; exit 1; }

if [[ ! -s "${TOPTAB}" ]]; then
  echo "[ERROR] Missing or empty: ${TOPTAB}" >&2
  exit 1
fi

# -----------------------------
# Helpers
# -----------------------------

# Find assembly file for a sample (tries common extensions)
find_asm() {
  local sample="$1"
  for ext in fasta fa fna; do
    if [[ -s "${ASM_DIR}/${sample}.${ext}" ]]; then
      echo "${ASM_DIR}/${sample}.${ext}"
      return 0
    fi
    if [[ -s "${ASM_DIR}/${sample}.${ext}.gz" ]]; then
      echo "${ASM_DIR}/${sample}.${ext}.gz"
      return 0
    fi
  done
  return 1
}

# Prepare uncompressed fasta + faidx
prepare_faidxable_fasta() {
  local sample="$1"
  local asm_path="$2"
  local tmp_fa="${OUTDIR}/tmp/${sample}.fa"

  if [[ -s "${tmp_fa}" ]]; then
    echo "${tmp_fa}"
    return 0
  fi

  if [[ "${asm_path}" =~ \.gz$ ]]; then
    echo "[INFO] gunzip -> ${tmp_fa}"
    gunzip -c "${asm_path}" > "${tmp_fa}"
  else
    cp -f "${asm_path}" "${tmp_fa}"
  fi

  samtools faidx "${tmp_fa}"
  echo "${tmp_fa}"
}

# Extract region (1-based inclusive)
faidx_region() {
  local fa="$1" contig="$2" start="$3" end="$4"
  samtools faidx "${fa}" "${contig}:${start}-${end}"
}

# Reverse complement FASTA stream (single record)
revcomp_fasta_stream() {
  if command -v seqtk >/dev/null 2>&1; then
    seqtk seq -r -
  else
    awk '
      BEGIN{hdr=""; seq=""}
      /^>/ {hdr=$0; next}
      {gsub(/[ \t\r]/,""); seq=seq $0}
      END{
        n=length(seq); out=""
        for(i=n;i>=1;i--){
          b=toupper(substr(seq,i,1))
          if(b=="A") c="T"
          else if(b=="T") c="A"
          else if(b=="C") c="G"
          else if(b=="G") c="C"
          else c="N"
          out=out c
        }
        print hdr
        for(i=1;i<=length(out);i+=60) print substr(out,i,60)
      }'
  fi
}

# Pick best hit by max(qend-qstart) for a given qname
# Output: sample qname qlen qstart qend strand tname tlen tstart tend mapq dv cm
best_hit_from_paf() {
  local paf="$1"
  local qname="$2"
  local sample="$3"

  awk -v Q="${qname}" -v S="${sample}" '
    BEGIN{best=-1; OFS="\t"}
    $1==Q {
      qlen=$2; qs=$3; qe=$4; strand=$5;
      tname=$6; tlen=$7; ts=$8; te=$9; mapq=$12;

      dv="NA"; cm="NA"
      for(i=13;i<=NF;i++){
        if($i ~ /^dv:f:/){ x=$i; sub(/^dv:f:/,"",x); dv=x }
        if($i ~ /^cm:i:/){ x=$i; sub(/^cm:i:/,"",x); cm=x }
      }

      aln=(qe-qs)
      if(aln>best){
        best=aln
        bestline=S OFS Q OFS qlen OFS qs OFS qe OFS strand OFS tname OFS tlen OFS ts OFS te OFS mapq OFS dv OFS cm
      }
    }
    END{
      if(best<0) exit 2
      print bestline
    }
  ' "${paf}"
}

# Convert PAF target coords (0-based, end-exclusive) to faidx coords (1-based inclusive),
# with flank and clamping.
calc_region_1based() {
  local ts0="$1" te0="$2" tlen="$3" flank="$4"
  local start1=$((ts0 + 1 - flank))
  local end1=$((te0 + flank))
  if (( start1 < 1 )); then start1=1; fi
  if (( end1 > tlen )); then end1="$tlen"; fi
  echo -e "${start1}\t${end1}"
}

###############################################################################
# MAIN
###############################################################################

TAIL_SUM="${OUTDIR}/summaries/tailocin_besthit_summary.tsv"
HTF_SUM="${OUTDIR}/summaries/htf_besthit_summary.tsv"

echo -e "sample\tqname\tqlen\tqstart\tqend\tq_aln\tq_cov\tstrand\ttname\ttlen\ttstart\ttend\tstart1\tend1\tmapq\tdv\tcm" > "${TAIL_SUM}"
echo -e "sample\tbest_HTFhap\tqname\tqlen\tqstart\tqend\tq_aln\tq_cov\tstrand\ttname\ttlen\ttstart\ttend\tstart1\tend1\tmapq\tdv\tcm" > "${HTF_SUM}"

TAIL_FASTA="${OUTDIR}/tailocin/tailocin_regions.fa"
HTF_FASTA="${OUTDIR}/htf/htf_plusminus${FLANK}.fa"
: > "${TAIL_FASTA}"
: > "${HTF_FASTA}"

# sample -> best HTF haplotype (col2 in TOPTAB)
HTF_MAP="${OUTDIR}/tmp/sample_to_bestHTF.tsv"
awk 'BEGIN{FS=OFS="\t"} NR==1{next} {print $1,$2}' "${TOPTAB}" > "${HTF_MAP}"

shopt -s nullglob
pafs=( "${PAF_DIR}"/*.paf )
echo "[INFO] Found ${#pafs[@]} PAF files."

for paf in "${pafs[@]}"; do
  bn=$(basename "${paf}")
  sample="${bn%.paf}"
  echo "[INFO] sample=${sample}"

  # assembly
  if ! asm=$(find_asm "${sample}"); then
    echo "[WARN] No assembly found for ${sample} in ASM_DIR=${ASM_DIR}. Skipping."
    continue
  fi
  fa=$(prepare_faidxable_fasta "${sample}" "${asm}")
  
  # -------------------------
  # (1) Tailocin extraction
  # -------------------------
  if best=$(best_hit_from_paf "${paf}" "${TAILOCIN_QNAME}" "${sample}" 2>/dev/null); then
    IFS=$'\t' read -r S Q qlen qs qe strand tname tlen ts te mapq dv cm <<< "${best}"
    q_aln=$((qe-qs))
    q_cov=$(awk -v a="${q_aln}" -v l="${qlen}" 'BEGIN{printf "%.6f", (l>0 ? a/l : 0)}')

    IFS=$'\t' read -r start1 end1 <<< "$(calc_region_1based "${ts}" "${te}" "${tlen}" 0)"

    echo -e "${S}\t${Q}\t${qlen}\t${qs}\t${qe}\t${q_aln}\t${q_cov}\t${strand}\t${tname}\t${tlen}\t${ts}\t${te}\t${start1}\t${end1}\t${mapq}\t${dv}\t${cm}" >> "${TAIL_SUM}"

    hdr=">${S}|tailocin|${tname}:${start1}-${end1}|strand=${strand}|qcov=${q_cov}|qaln=${q_aln}/${qlen}"
    if [[ "${strand}" == "-" ]]; then
      { echo "${hdr}"; faidx_region "${fa}" "${tname}" "${start1}" "${end1}" | tail -n +2; } \
        | revcomp_fasta_stream >> "${TAIL_FASTA}"
    else
      echo "${hdr}" >> "${TAIL_FASTA}"
      faidx_region "${fa}" "${tname}" "${start1}" "${end1}" | tail -n +2 >> "${TAIL_FASTA}"
    fi
    echo >> "${TAIL_FASTA}"
  else
    echo "[WARN] No tailocin hit in ${bn} for ${TAILOCIN_QNAME}"
  fi

  # -------------------------
  # (2) HTF +/-200 extraction
  # -------------------------
  best_htf=$(awk -v s="${sample}" 'BEGIN{FS="\t"} $1==s {print $2; exit}' "${HTF_MAP}" || true)
  if [[ -z "${best_htf}" ]]; then
    echo "[WARN] No best HTF haplotype for ${sample} in TOPTAB"
    continue
  fi
  HTF_QNAME="HTF_${best_htf}"

  if best2=$(best_hit_from_paf "${paf}" "${HTF_QNAME}" "${sample}" 2>/dev/null); then
    IFS=$'\t' read -r S2 Q2 qlen2 qs2 qe2 strand2 tname2 tlen2 ts2 te2 mapq2 dv2 cm2 <<< "${best2}"
    q_aln2=$((qe2-qs2))
    q_cov2=$(awk -v a="${q_aln2}" -v l="${qlen2}" 'BEGIN{printf "%.6f", (l>0 ? a/l : 0)}')

    IFS=$'\t' read -r start2 end2 <<< "$(calc_region_1based "${ts2}" "${te2}" "${tlen2}" "${FLANK}")"

    echo -e "${S2}\t${best_htf}\t${Q2}\t${qlen2}\t${qs2}\t${qe2}\t${q_aln2}\t${q_cov2}\t${strand2}\t${tname2}\t${tlen2}\t${ts2}\t${te2}\t${start2}\t${end2}\t${mapq2}\t${dv2}\t${cm2}" >> "${HTF_SUM}"

    hdr2=">${S2}|HTF=${best_htf}|${tname2}:${start2}-${end2}|strand=${strand2}|qcov=${q_cov2}|qaln=${q_aln2}/${qlen2}"
    if [[ "${strand2}" == "-" ]]; then
      { echo "${hdr2}"; faidx_region "${fa}" "${tname2}" "${start2}" "${end2}" | tail -n +2; } \
        | revcomp_fasta_stream >> "${HTF_FASTA}"
    else
      echo "${hdr2}" >> "${HTF_FASTA}"
      faidx_region "${fa}" "${tname2}" "${start2}" "${end2}" | tail -n +2 >> "${HTF_FASTA}"
    fi
    echo >> "${HTF_FASTA}"
  else
    echo "[WARN] No HTF hit in ${bn} for ${HTF_QNAME}"
  fi

done

echo "[DONE] Tailocin FASTA : ${TAIL_FASTA}"
echo "[DONE] Tailocin summary: ${TAIL_SUM}"
echo "[DONE] HTF FASTA      : ${HTF_FASTA}"
echo "[DONE] HTF summary    : ${HTF_SUM}"

###############################################################################
# DEBUG: verify faidx extraction really returned sequence
###############################################################################
debug_out="${OUTDIR}/tmp/debug_${sample}_${tname}_${start1}_${end1}.txt"

samtools faidx "${fa}" "${tname}:${start1}-${end1}" > "${debug_out}" 2>&1

echo "[DEBUG] faidx output for ${sample} ${tname}:${start1}-${end1}" >&2
sed -n '1,5p' "${debug_out}" >&2

seq_lines=$(grep -v '^>' "${debug_out}" | tr -d '\n' | wc -c)

if (( seq_lines == 0 )); then
  echo "[ERROR] EMPTY SEQUENCE from faidx!" >&2
  echo "        sample=${sample}" >&2
  echo "        fasta=${fa}" >&2
  echo "        contig=${tname}" >&2
  echo "        region=${start1}-${end1}" >&2
  echo "        FASTA contigs (first 10):" >&2
  grep '^>' "${fa}" | head -n 10 >&2
  echo "        PAF contigs (unique):" >&2
  cut -f6 "${paf}" | sort -u | head -n 10 >&2
  continue
fi
###############################################################################
