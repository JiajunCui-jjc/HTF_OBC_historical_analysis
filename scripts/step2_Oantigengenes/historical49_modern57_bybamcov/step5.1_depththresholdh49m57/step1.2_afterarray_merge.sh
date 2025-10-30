   cd "/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step2_Oantigengenes/historical49_modern57_bybamcov/coverage_h49_m57"
   echo -e "Sample\tMean_total_depth\tSD_total_depth\tMean_minus_025SD\trelative_threshold" > depth_thresholds_mean_minus_025_sd.tsv
   cat tmp/tmp_*.tsv >> depth_thresholds_mean_minus_025_sd.tsv
   rm -r tmp

