#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   step2_onebyoneretrieve.sh <contigs.fasta> <out.fa> <best_HTF> <contig:start-end:strand> [<contig:start-end:strand> ...]
#
# Example:
#   step2_onebyoneretrieve.sh \
#     results/assemblies_tmp/assemblies/109.NOR_1990/contigs.fasta \
#     ./results/seqs/109.NOR_1990_fulllength_HTFp5.D5.fa \
#     HTF_p5.D5 \
#     NODE_1_length_19856_cov_20.341724:4150-5955:+
###full length: 3 fully mapped, 2 can be recovered:
#all checked by comparing with HTF refs by eye (start and end few bp)

#three full mapped:
#109.NOR_1990    p5.D5   1       longcontig,fully mapped by one contig
#76.LTU_2009_S19 p21.F9  1       longcontig,fully mapped by one contig
#PL0059  p21.F9  1       longcontig,fully mapped by one contig, i guess coinfected two strains red and yellow (p21.F9 1245bp), kmer said its red since deeper kmer read. but assembly show the p21.F9 is complete...
#bash step2_onebyoneretrieve.sh  results/assemblies_tmp/assemblies/109.NOR_1990/contigs.fasta ./results/seqs/109.NOR_1990_fulllength_HTFp5.D5.fa HTF_p5.D5 NODE_1_length_19856_cov_20.341724:4150-5955:+
#bash step2_onebyoneretrieve.sh results/assemblies_tmp/assemblies/76.LTU_2009_S19/contigs.fasta ./results/seqs/76.LTU_2009_S19_full_HTFp21.F9_1245.fa HTF_p21.F9 NODE_1_length_13003_cov_61.548805:202-1449:+
#bash step2_onebyoneretrieve.sh results/assemblies_tmp/assemblies/PL0059/contigs.fasta ./results/seqs/PL0059_fulllength_HTFp21.F9_1245.fa HTF_p21.F9 NODE_1_length_15224_cov_240.485748:13978-15224:-
#15224 is the end so no way to buffer 1bp

#two can be recoverd: 
#PL0042  p25.A12 0.993492        longcontig mappped from 9 to 1383, can be recovered to full length
#bash step2_onebyoneretrieve.sh results/assemblies_tmp/assemblies/PL0042/contigs.fasta ./results/seqs/PL0042_9-1383_p25.A12.fa HTF_p25.A12 NODE_1_length_10427_cov_100.946411:3976-5361:+  
#PL0102  p25.A12 0.993492        1800bpcontig mappped from 9 to 1383, can be recovered to full length
#bash step2_onebyoneretrieve.sh results/assemblies_tmp/assemblies/PL0102/contigs.fasta ./results/seqs/PL0102_9-1383_p25.A12.fa HTF_p25.A12 NODE_1_length_1812_cov_16.274311:319-1704:-

#five that mapped by 2 short contigs, no need to recover since none of them could be fully length recovered so cant check with ORF, keep the mapping region as it was, no buffer so the faidx coordinate is from mapping paf start+1-end
 
#NODE_11_length_843_cov_5.900000 843     8       843
#NODE_8_length_943_cov_4.379121  943     551     943
#PL0240  p21.F9  0.985542        by two contigs, contig1 from 5:840 and contig two from 853 to 1245, the first 5bp could be recovered by extend contig1 but cant get 840-853 but not necessary
#bash step2_onebyoneretrieve.sh results/assemblies_tmp/assemblies/PL0240/contigs.fasta ./results/seqs/PL0240_5-840_853-1245_p21.F9.fa HTF_p21.F9 NODE_11_length_843_cov_5.900000:8-843:- NODE_8_length_943_cov_4.379121:551-943:- 
#PL0131  p5.D5   0.921242        by two contigs, 1 for 0-1257, the contig length is 4000, so we can recover from it the full length,  2 for 1399-1803 only cant be recovered from 1359-1399
#+       NODE_2_length_4225_cov_44.341126        4225    2883    4140 
#HTF_p5.D5       1803    1399    1803    +       NODE_6_length_1002_cov_36.621259        1002    40      444 
#bash step2_onebyoneretrieve.sh results/assemblies_tmp/assemblies/PL0131/contigs.fasta ./results/seqs/PL0131_0-1257_1399-1803_p5.D5.fa HTF_p5.D5 NODE_2_length_4225_cov_44.341126:2884-4140:+ NODE_6_length_1002_cov_36.621259:41-444:+ NODE_2_length_4225_cov_44.341126:2882-4687:+   
#HB0766  p21.F9  0.914859        by two contigs, can recover 40 and 10 in each but not necessary since they are unmapped and there is overlap so we cant convince by ORF
#HTF_p21.F9      1245    586     1245    -       NODE_4_length_1350_cov_5.321944 1350    651     1310    652     659     60      
#HTF_p21.F9      1245    10      490     -       NODE_8_length_490_cov_6.752735  490     10      490     476     480     60      
#bash step2_onebyoneretrieve.sh results/assemblies_tmp/assemblies/HB0766/contigs.fasta ./results/seqs/HB0766_10-490_586-1245_p21.F9.fa HTF_p21.F9 NODE_8_length_490_cov_6.752735:11-490:- NODE_4_length_1350_cov_5.321944:652-1310:-
#PL0046  p21.F9  0.739759        by two contigs, 13:490 and 810:1245, can recover 3:490+22 but not necessary... increasing uncertainty
#HTF_p21.F9      1245    801     1245    +       NODE_6_length_939_cov_13.567329 939     0       444     435     444     60      
#HTF_p21.F9      1245    13      490     -       NODE_12_length_509_cov_16.420168        509     22      499     465     477     60     
#bash step2_onebyoneretrieve.sh results/assemblies_tmp/assemblies/PL0046/contigs.fasta ./results/seqs/PL0046_10-490_801-1245_p21.F9.fa HTF_p21.F9 NODE_12_length_509_cov_16.420168:23-499:- NODE_6_length_939_cov_13.567329:1-444:+
#PL0080  p25.A12 0.715112        by two contigs, 0-543 can extend to 0-543+24, 937-1383 extend to 929-1383 but again not necessary...
#HTF_p25.A12     1383    0       543     -       NODE_4_length_1062_cov_22.946550        1062    24      567     540     543     60      NM:i:3  ms:i:483        AS:i:483
#HTF_p25.A12     1383    937     1383    -       NODE_2_length_4830_cov_25.020221        4830    4371    4817    437     446     60      NM:i:9  ms:i:266        AS:i:266
#bash step2_onebyoneretrieve.sh results/assemblies_tmp/assemblies/PL0080/contigs.fasta ./results/seqs/PL0080_0-543_937-1383_p25.A12.fa HTF_p25.A12 NODE_4_length_1062_cov_22.946550:25-567:- NODE_2_length_4830_cov_25.020221:4372-4817:- 

if [[ $# -lt 4 ]]; then
  echo "Usage: $0 <contigs.fasta> <out.fa> <best_HTF> <contig:start-end:strand> ..." >&2
  exit 1
fi

CONTIGS=$1
OUTFA=$2
BEST_HTF=$3
shift 3  # remaining args are regions

# ---- deps ----
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found in PATH" >&2; exit 1; }; }
need samtools
need seqtk

# ---- index ----
if [[ ! -f "${CONTIGS}.fai" ]]; then
  samtools faidx "$CONTIGS"
fi

# ---- infer sample name from OUTFA (basename without extension; strip trailing _... if you like) ----
bn=$(basename "$OUTFA")
sample=${bn%.*}          # drop .fa/.fasta
# Optional: if you prefer the bit before the first underscore as "sample", uncomment:
# sample=${sample%%_*}

# ---- start fresh output ----
: > "$OUTFA"

for region in "$@"; do
  contig=$(echo "$region" | cut -d: -f1)
  coords=$(echo "$region" | cut -d: -f2)
  strand=$(echo "$region" | cut -d: -f3)

  if [[ -z "${strand:-}" ]]; then
    echo "Missing strand info for region '$region' (must be + or -)" >&2
    exit 1
  fi

  # Build the desired header
  # 1-based inclusive coords are preserved as given in $coords
  header="${sample}|${BEST_HTF}|${contig}:${coords}(${strand})"

  if [[ "$strand" == "+" ]]; then
    samtools faidx "$CONTIGS" "${contig}:${coords}" \
      | awk -v h="$header" 'NR==1{$0=">"h}1' >> "$OUTFA"
  elif [[ "$strand" == "-" ]]; then
    samtools faidx "$CONTIGS" "${contig}:${coords}" \
      | seqtk seq -r - \
      | awk -v h="$header" 'NR==1{$0=">"h}1' >> "$OUTFA"
  else
    echo "Invalid strand '$strand' for region '$region' (must be + or -)" >&2
    exit 1
  fi
done

echo "Extracted regions written to $OUTFA"
