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
# USER CONFIG
###############################################################################
ASM_DIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2026_PNAS_revision/data/modern53_fasta"
OUTDIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2026_PNAS_revision/results/step2_tailocin_HTF_seqs/tailocin"

TAILOCIN_QNAME="utg000001l_p25.C2:2695984-2717365"
QLEN=21382   # from your pasted PAF lines

mkdir -p "${OUTDIR}"/{tmp,frags}

OUTFA="${OUTDIR}/tailocin_5samples_concat.fa"
: > "${OUTFA}"

###############################################################################
# Helpers
###############################################################################
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

prepare_faidxable_fasta() {
  local sample="$1"
  local asm_path="$2"
  local tmp_fa="${OUTDIR}/tmp/${sample}.fa"

  if [[ -s "${tmp_fa}" ]]; then
    echo "${tmp_fa}"
    return 0
  fi

  if [[ "${asm_path}" =~ \.gz$ ]]; then
    echo "[INFO] gunzip -> ${tmp_fa}" >&2
    gunzip -c "${asm_path}" > "${tmp_fa}"
  else
    echo "[INFO] copy -> ${tmp_fa}" >&2
    cp -f "${asm_path}" "${tmp_fa}"
  fi

  samtools faidx "${tmp_fa}" >&2
  echo "${tmp_fa}"
}

# Extract contig region (1-based inclusive) as SEQUENCE ONLY (no FASTA header), 1-line
extract_seq_oneline() {
  local fa="$1" contig="$2" start1="$3" end1="$4"
  samtools faidx "${fa}" "${contig}:${start1}-${end1}" \
    | awk 'NR>1{gsub(/[ \t\r]/,""); printf "%s",$0} END{printf "\n"}'
}

# Reverse-complement a one-line DNA sequence
revcomp_oneline() {
  if command -v seqtk >/dev/null 2>&1; then
    # feed as FASTA to seqtk, output seq only
    awk 'BEGIN{print ">x"} {print $0}' \
      | seqtk seq -r - \
      | awk 'NR>1{gsub(/[ \t\r]/,""); printf "%s",$0} END{printf "\n"}'
  else
    awk '
      {
        seq=$0
        n=length(seq)
        out=""
        for(i=n;i>=1;i--){
          b=toupper(substr(seq,i,1))
          if(b=="A") c="T"
          else if(b=="T") c="A"
          else if(b=="C") c="G"
          else if(b=="G") c="C"
          else c="N"
          out=out c
        }
        print out
      }'
  fi
}

###############################################################################
# Hit table (hard-coded from your pasted lines)
#
# Columns:
# sample   qstart  qend   strand  tname  tstart0  tend0
#
# NOTE:
# - PAF tstart/tend are 0-based, end-exclusive.
# - We convert to 1-based inclusive for samtools faidx:
#     start1 = tstart0 + 1
#     end1   = tend0
###############################################################################
HITS_TSV="${OUTDIR}/tmp/5samples_tailocin_hits.tsv"
cat > "${HITS_TSV}" << 'EOF'
p25.B2	7	9392	+	877	14888	24274
p25.B2	9393	11900	+	928	13	2520
p25.B2	11939	21321	+	549	7	9389
p25.D2	7	10872	-	284	6	10872
p25.D2	10915	21368	+	662	60	10513
p26.D6	30	14543	-	524	2	14516
p26.D6	14660	21341	+	403	15	7772
p25.A12	30	10872	-	7	10	9197
p25.A12	11256	21341	+	136	274	10180
p25.C11	30	10934	+	1259	14890	26071
p25.C11	10936	20135	-	178	11	9211
p25.C11	20149	21341	+	1489	5	1016
EOF

###############################################################################
# Main: for each sample, extract fragments in qstart order, normalize strand,
# concatenate with NO Ns, write one FASTA record.
###############################################################################
for sample in p25.B2 p25.D2 p26.D6 p25.A12 p25.C11; do
  echo "[INFO] Processing ${sample}" >&2

  if ! asm=$(find_asm "${sample}"); then
    echo "[ERROR] Assembly not found for ${sample} in ${ASM_DIR}" >&2
    exit 1
  fi
  fa=$(prepare_faidxable_fasta "${sample}" "${asm}")

  # sort fragments by qstart, extract each, RC if needed, append to concat string
  concat_seq=""
  qaln_sum=0

  while IFS=$'\t' read -r s qs qe strand tname ts0 te0; do
    [[ "${s}" == "${sample}" ]] || continue

    hitlen=$((qe-qs))
    qaln_sum=$((qaln_sum + hitlen))

    start1=$((ts0 + 1))
    end1=$((te0))

    # Extract fragment sequence (one-line)
    frag=$(extract_seq_oneline "${fa}" "${tname}" "${start1}" "${end1}")

    # Normalize orientation: if strand is '-', reverse-complement the target fragment
    if [[ "${strand}" == "-" ]]; then
      frag=$(printf "%s\n" "${frag}" | revcomp_oneline)
    fi

    # Concatenate (no Ns)
    concat_seq="${concat_seq}${frag}"

  done < <(awk -v S="${sample}" 'BEGIN{FS=OFS="\t"} $1==S{print}' "${HITS_TSV}" | sort -k2,2n)

  # Compute qcov
  qcov=$(awk -v a="${qaln_sum}" -v l="${QLEN}" 'BEGIN{printf "%.6f", (l>0 ? a/l : 0)}')

  # Header format you requested
  # For multi-fragment, we mark tname/start/end as "concat" and strand as "+"
  hdr=">${sample}|tailocin|concat:1-${QLEN}|strand=+|qcov=${qcov}|qaln=${qaln_sum}/${QLEN}"

  # Write FASTA (sequence on one line)
  {
    echo "${hdr}"
    echo "${concat_seq}"
  } >> "${OUTFA}"

done

echo "[DONE] Wrote multi-FASTA: ${OUTFA}" >&2
echo "[INFO] Records: $(grep -c '^>' "${OUTFA}")" >&2
