#ln -s /SAN/ugi/plant_genom/jiajucui/answer2025_maptoOTU5_nohaplotype_h49/h10_lopez/* /SAN/ugi/plant_genom/jiajucui/answer2025_maptoOTU5_nohaplotype_h49/

wd=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/scripts/step0_datapreprocessing_tailocin_presence/step1_extract_h49genomes_stats/step1_historical49_extractreads_from_metagenome/step2_generate_stat_table
echo 'samplename At_percent At_covered At_depth Ps_percent Ps_covered Ps_depth'  >> $wd/h49dedup_maptoOTU5withouthaplotype_stats.txt 

#samplename At_percent (before dedup) At_covered (after dedup) At_depth (after dedup) Ps_percent (before dedup, exclude haplotypes) Ps_covered (after dedup, exclude haplotypes) Ps_depth (after dedup, exclude haplotypes)
sample=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/historical49.txt
for line in $(cat $sample | sed -n $i'p'| awk '{print $1}');
do
bash  $wd/step1.2_generatetable_examplecode.sh $line >> $wd/h49dedup_maptoOTU5withouthaplotype_stats.txt;
done
cp $wd/h49dedup_maptoOTU5withouthaplotype_stats.txt /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step1_h49metadata/

