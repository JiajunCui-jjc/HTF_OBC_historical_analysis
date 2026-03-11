#!/bin/bash -l

#$ -l tmem=20G
#$ -l h_vmem=20G
#$ -l h_rt=10:00:0
#$ -wd /SAN/ugi/plant_genom/jiajucui/
#$ -V
#$ -N blastp
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/




# =========================
# User inputs (edit if needed)
# =========================
blastp=/share/apps/ncbi-blast-2.11.0+/bin/blastp
makeblastdb=/share/apps/ncbi-blast-2.11.0+/bin/makeblastdb

queryfasta_53="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2026_PNAS_revision/results/step2_tailocin_HTF_seqs/htf/runinlocal_htf_real_length/results/real_HTF_aa.fasta"
ref_faa="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2026_PNAS_revision/data/ref_fault_T1T2_fromPattricia/ref10.faa"
outdir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2026_PNAS_revision/results/step3_blastp"

threads=4

# =========================
# Type mapping (from your legend)
#   - Blue = Type1 (T1)
#   - Orange = Type2 (T2)
# NOTE: If your ref10.faa contains additional WP_ accessions not listed here,
#       they will be labeled "UNK".
# =========================

mkdir -p "$outdir"
# 1) use node-local scratch if available

# 2) avoid huge LMDB mmap allocation

#"$makeblastdb" -in "$ref" -dbtype prot -out "$dbprefix" -parse_seqids
#fail to alllocate mem so use ref directly



# 2) Run BLASTP (tabular)
# Keep multiple hits; we will sort & take top6 per query later.
# -------------------------
raw="${outdir}/blastp_ref10_vs_queries.raw.tsv"
"$blastp" \
  -query "$queryfasta_53" \
  -subject "$ref_faa" \
  -max_target_seqs 200 \
  -evalue 1e-5 \
  -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore" \
  > "$raw"
#evalue threshold liek Pattricia
# -------------------------
# 3) For each query, keep TOP 6 hits by bitscore desc (then pident desc)



raw="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2026_PNAS_revision/results/step3_blastp/blastp_ref10_vs_queries.raw.tsv"

out="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2026_PNAS_revision/results/step3_blastp/blastp_best_hit.wide.tsv"

T1_LIST="WP_024649699.1 WP_024658146.1 WP_024693886.1 WP_032656593.1"
T2_LIST="WP_005768002.1 WP_024639976.1 WP_024674765.1"
#T3 WP_005896706.1 WP_027901668.1 WP_027901669.1
LC_ALL=C sort -s \
  -k1,1 \
  -k8,8gr \
  -k3,3gr \
  "$raw" \
| awk -v OFS="\t" -v T1_LIST="$T1_LIST" -v T2_LIST="$T2_LIST" '
BEGIN{
  n1=split(T1_LIST,a," "); for(i in a) t1[a[i]]=1
  n2=split(T2_LIST,b," "); for(i in b) t2[b[i]]=1

  print "query","best_ref_type","best_ref_name","aligned_len/ref_len","best_identity","best_bitscore"
}
function type_of(s){
  if (s in t1) return "T1"
  if (s in t2) return "T2"
  return "T3"
}
{
  q=$1
  if (q in seen) next   # first hit per query = best by bitscore
  seen[q]=1

  ref=$2
  pid=$3
  aln=$4
  slen=$6
  bits=$8

  printf "%s\t%s\t%s\t%d/%d\t%.3f\t%.1f\n", \
         q, type_of(ref), ref, aln, slen, pid, bits
}
' > "$out"



echo "Done."
echo "RAW:    $raw"
echo "WIDE:   $out"
