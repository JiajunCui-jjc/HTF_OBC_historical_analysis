#!/bin/bash
#ls
#excludeselfHTF  HTF_7refs  wholegenome_7refs
#(base) [jiajucui@pchuckle forkmer_step0_7refs]$ pwd
#/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/HTFreference/forkmer_step0_7refs
wd=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/HTFreference/forkmer_step0_7refs
for genome in $wd/wholegenome_7refs/*.fasta; do
    sample=$(basename "$genome" .fasta)
    # Run minimap2 alignment
    minimap2 -x asm5 "$genome"  $wd/HTF_7refs/HTF_${sample}.fa > $wd/excludeselfHTF/step1_HTFpaf_7wholegenomes/"${sample}.wg.exclude_selfHTF.paf"
done

