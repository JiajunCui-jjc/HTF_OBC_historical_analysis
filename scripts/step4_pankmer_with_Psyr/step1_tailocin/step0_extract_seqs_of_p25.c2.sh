ref=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/HTFreference/OTU5_ref
wd=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2026_PNAS_revision/results/step0_refseqs
#trpE
#utg000001l_p25.C2:2694503-2695984
#trpG
#utg000001l_p25.C2:2717365-2717964
#old
#utg000001l_p25.C2:2696041-2716969
#new
#utg000001l_p25.C2:2695984-2717365
#conda activate phylogeny_snp
cd $wd
samtools faidx $ref/ref.withhyplotype.fasta utg000001l_p25.C2:2696041-2716969 | seqtk seq -r - > $wd/tailocinregion_HTFTFA_p25.c2.fasta
seqkit grep -r -p 'HTF_|TFA_' $ref/ref.withhyplotype.fasta >> $wd/tailocinregion_HTFTFA_p25.c2.fasta


#new between trpEand trpG
samtools faidx $ref/ref.withhyplotype.fasta utg000001l_p25.C2:2695984-2717365 | seqtk seq -r - > $wd/betweentrpEG_tailocinregion_HTFTFA_p25.c2.fasta
seqkit grep -r -p 'HTF_|TFA_' $ref/ref.withhyplotype.fasta >> $wd/betweentrpEG_tailocinregion_HTFTFA_p25.c2.fasta

