Tailocin region sequences (53 modern assemblies) – between trpE and trpG
======================================================================

Overview
--------
This folder contains extracted tailocin-region DNA sequences from 53 modern assemblies.
The extraction uses minimap2 PAF mappings of a single tailocin reference query
(from p25.C2, region between trpE end and trpG start):

  Query name (minimap2): utg000001l_p25.C2:2695984-2717365
  Query length (qlen):   21382 bp

Each modern assembly was searched for the best match to this query. The resulting
tailocin-region sequences were exported as a multi-FASTA for downstream inspection
and (later) multiple sequence alignment / phylogenetic analyses.


Files
-----
1) final_tailocin_all53_above98percentcov.fa
   - Multi-FASTA containing tailocin-region sequences for all 53 samples.
   - Sequence orientation is normalized:
       * If the minimap2 hit strand is '-', the extracted target segment is reverse-complemented.
       * Otherwise it is kept as-is.
   - Records are unaligned (raw sequences).

2) final_header.txt
   - One-line-per-record headers extracted from the FASTA.
   - Useful for quickly checking coordinates, strand, and coverage without scanning the full FASTA.
   - Header format:

       >SAMPLE|tailocin|TNAME:START-END|strand=STR|qcov=QCOV|qaln=QALN/QLEN

     where:
       SAMPLE  = sample / assembly ID (e.g. p12.A11)
       TNAME   = target contig name in that assembly
       START-END = extracted target coordinates (1-based inclusive)
       STR     = '+' or '-' from minimap2 alignment to the query
       QALN    = query-aligned span (bp), computed as (qend - qstart)
       QLEN    = 21382
       QCOV    = QALN / QLEN

3) readme.txt
   - This document.



How sequences were chosen
-------------------------
For each sample, minimap2 outputs one or more alignments of the query to the assembly.

A) 48 samples: single-contig best hit (quality-filtered)
   - For these samples, a single best minimap2 hit on one contig provides >98% query coverage.
   - “Best hit” is defined as the largest query-aligned span (qend - qstart).
   - These 48 sequences come directly from that single contig interval.

B) 5 samples: multi-hit / multi-contig reconstruction (concatenated fragments)
   - For the following samples, the query coverage is achieved by multiple non-overlapping
     query segments mapping to different contigs/locations.
   - For these samples, fragments were extracted per hit, strand-normalized, ordered by query
     coordinate (qstart), and concatenated with NO padding (no inserted Ns).

   Multi-hit samples summary (hit lengths = qend - qstart; cov_prop = sum/21382):

     sample    hit_sum_expression              cov_prop
     p25.B2    9385 + 9382 + 2507              0.9949
     p25.D2    10865 + 10453                   0.9970
     p26.D6    14513 + 6681                    0.9912
     p25.A12   10085 + 10842                   0.9787
     p25.C11   9199 + 10904 + 1192             0.9959

