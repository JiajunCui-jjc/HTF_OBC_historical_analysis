# temporarily disable nounset, activate, then re-enable
set +u
source /home/jiajucui/miniconda3/bin/activate phylogeny_snp

set -u   # optional: turn nounset back on

# Activate the conda environment


# Usage:
#   ./run_iqtree_simple.sh <aligned_alignment> [MODEL] [THREADS] [SEED]
# Examples:
#   ./run_iqtree_simple.sh clustalo-I20250907-023650-0165-19592007-p1m.aln-fasta
#   ./run_iqtree_simple.sh clustal.aln GTR+F+R 8 12345
#
# Notes:
# - IQ-TREE auto-detects alignment format (FASTA/CLUSTAL/PHYLIP/NEXUS).
# - Outputs: *.contree (tree with supports), *.treefile (best ML tree), *.iqtree (run log).

ALIGN="$1"
MODEL="${2:-MFP}"     # MFP lets IQ-TREE pick best model; or use GTR+F+R for speed
THREADS="${3:-4}"  # or set e.g. 8
SEED="${4:-12345}"


# Run IQ-TREE with ultrafast bootstrap and SH-aLRT supports
iqtree \
  -s "$ALIGN" \
  -m "$MODEL" \
  -B 1000 --alrt 1000 \
  -T "$THREADS" \

echo "Done."
echo "Key outputs next to $ALIGN:"
echo "  $(basename "$ALIGN").contree   # ML tree with UFBoot + SH-aLRT supports"
echo "  $(basename "$ALIGN").treefile  # best ML tree"
echo "  $(basename "$ALIGN").iqtree    # run log + selected model"
