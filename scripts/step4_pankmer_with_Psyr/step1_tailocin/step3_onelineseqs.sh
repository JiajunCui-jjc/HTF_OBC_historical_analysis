#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Convert multi-line FASTA to 2-line FASTA (header + sequence on one line)
# Usage:
#   ./fasta_to_oneline.sh <input.fa> <output.fa>
#
# Also includes your two fixed inputs by default (tailocin + HTF) if run with no args.
###############################################################################

oneline() {
  local in_fa="$1"
  local out_fa="$2"

  if [[ ! -s "${in_fa}" ]]; then
    echo "[ERROR] Input FASTA missing/empty: ${in_fa}" >&2
    exit 1
  fi

  echo "[INFO] Converting to one-line FASTA: ${in_fa} -> ${out_fa}" >&2

  if command -v seqtk >/dev/null 2>&1; then
    # -l 0 means no wrapping: sequence becomes one line per record
    seqtk seq -l 0 "${in_fa}" > "${out_fa}"
  else
    # AWK fallback: concatenate sequence lines per record
    awk '
      BEGIN{hdr=""; seq=""}
      /^>/{
        if(hdr!=""){
          print hdr
          print seq
        }
        hdr=$0
        seq=""
        next
      }
      {
        gsub(/[ \t\r]/,"")
        seq=seq $0
      }
      END{
        if(hdr!=""){
          print hdr
          print seq
        }
      }
    ' "${in_fa}" > "${out_fa}"
  fi

  # QC: any header-only records?
  bad=$(awk '
    /^>/{ if(NR>1 && seqlen==0) bad++; seqlen=0; next }
    { gsub(/[ \t\r]/,""); seqlen += length($0) }
    END{ if(NR>0 && seqlen==0) bad++; print bad+0 }
  ' "${out_fa}")

  if [[ "${bad}" != "0" ]]; then
    echo "[WARN] Found ${bad} header-only record(s) in ${out_fa}" >&2
  fi

  echo "[DONE] Wrote: ${out_fa}" >&2
}

###############################################################################
# Main
###############################################################################
if [[ $# -eq 2 ]]; then
  oneline "$1" "$2"
  exit 0
fi

if [[ $# -ne 0 ]]; then
  echo "[USAGE] $0 <input.fa> <output.fa>" >&2
  echo "        Or run with no args to process the default tailocin + HTF files." >&2
  exit 1
fi

BASE="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2026_PNAS_revision/results/step2_tailocin_HTF_seqs"

TAIL_IN="${BASE}/tailocin/tailocin_regions.fa"
HTF_IN="${BASE}/htf/htf_plusminus200.fa"

TAIL_OUT="${BASE}/tailocin/tailocin_regions.oneline.fa"
HTF_OUT="${BASE}/htf/htf_plusminus200.oneline.fa"

oneline "${TAIL_IN}" "${TAIL_OUT}"
oneline "${HTF_IN}"  "${HTF_OUT}"
