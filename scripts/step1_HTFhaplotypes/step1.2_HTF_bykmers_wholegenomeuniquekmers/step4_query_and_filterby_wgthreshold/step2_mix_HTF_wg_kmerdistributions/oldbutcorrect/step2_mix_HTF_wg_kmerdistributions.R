# === Libraries ===
library(data.table)
library(ggplot2)
library(scales)
library(dplyr)

# === Paths ===
htf_data_path <- "/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers/kmer_distribution_mix/HTF_kmer_depth_distributions.tsv"
historical_wg_dir <- "/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step2_isolate_kmers/isolate_jf/dump"
modern_wg_dir <- "/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step2_isolate_kmers/isolate_jf/dump"
out_dir <- "/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers/kmer_distribution_mix/distributions"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# === Thresholds ===
threshold_modern <- fread("/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step2_isolate_kmers/isolate_jf/dump/threshold_m57.tsv", col.names = c("sample", "cutoff"))
threshold_historical <- fread("/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step2_isolate_kmers/isolate_jf/dump/threshold_h49_sd_twosides.tsv", col.names = c("sample", "mean_minus_sd","mean_plus_sd"))
threshold_map_modern <- setNames(threshold_modern$cutoff, threshold_modern$sample)
threshold_map_historical <- setNames(threshold_historical$mean_minus_sd, threshold_historical$sample)

# === Colors ===
htf_colors <- c(
  "HTF_p21.F9" = "#66C2A5",
  "HTF_p23.B8" = "#FFD92F",
  "HTF_p25.A12" = "#A6A6D2",
  "HTF_p25.C2" = "#FC8D62",
  "HTF_p26.D6" = "#8DA0CB",
  "HTF_p5.D5"  = "#E78AC3",
  "HTF_p7.G11" = "#A6D854"
)

# === Load HTF k-mer depth table ===
df_depths <- fread(htf_data_path, col.names = c("sample", "ref", "kmer", "depth"))
df_depths$display_title <- df_depths$sample

# === Process per sample ===
sample_list <- unique(df_depths$sample)

for (sample_name in sample_list) {
  cat("🔧 Processing:", sample_name, "\n")
  
  df_sample <- df_depths[sample == sample_name & depth >= 1]
  if (nrow(df_sample) == 0) {
    cat("⚠️  No HTF data for", sample_name, "\n")
    next
  }
  
  # === Find WG k-mer dump ===
  wg_path <- file.path(historical_wg_dir, paste0(sample_name, "_dump.txt"))
  if (!file.exists(wg_path)) {
    wg_path <- file.path(modern_wg_dir, paste0(sample_name, "_dump_20p_dump.txt"))
  }
  if (!file.exists(wg_path)) {
    cat("❌ No WG file for", sample_name, "\n")
    next
  }
  
  # === Load WG k-mers ===
  # Load WG k-mers
  df_wg <- fread(wg_path, col.names = c("kmer", "count"))
  df_wg <- df_wg[count > 2]  # remove singleton kmers

  if (nrow(df_wg) == 0) {
    cat("❌ WG file has no informative k-mers for", sample_name, "\n")
    next
  }

  # === Right-tail filtering: count ≤ mean + 2×SD ===
  depth_mean <- mean(df_wg$count)
  depth_sd   <- sd(df_wg$count)
  cutoff_max <- depth_mean + 3 * depth_sd
  df_wg <- df_wg[count <= cutoff_max]

  if (nrow(df_wg) == 0) {
    cat("❌ WG file has no k-mers left after filtering for", sample_name, "\n")
    next
  }

  #   Frequency table
  freq_df_wg <- df_wg[, .N, by = count]
  setnames(freq_df_wg, c("depth", "frequency")) 
  # === Determine threshold ===
  if (sample_name %in% names(threshold_map_modern)) {
    threshold_val <- threshold_map_modern[[sample_name]]
  } else if (sample_name %in% names(threshold_map_historical)) {
    threshold_val <- threshold_map_historical[[sample_name]]
  } else {
    cat("❗ No threshold available for", sample_name, "\n")
    threshold_val <- NA
  }
  
  # === Scaling for dual Y-axis ===
  max_wg_y <- max(freq_df_wg$frequency)
  max_htf_y <- max(table(df_sample$depth))
  if (max_htf_y == 0) {
    cat("⚠️  No HTF bins for", sample_name, "\n")
    next
  }
  scale_factor <- max_wg_y / max_htf_y
  
  # === Plot ===
  title <- paste0("HTF & WG k-mer depth for ", unique(df_sample$display_title))
  
  p <- ggplot() +
    # Whole-genome (background grey)
    geom_col(data = freq_df_wg, aes(x = depth, y = frequency),
             fill = "grey50", alpha = 0.3, show.legend = TRUE) +
    
    # HTF overlay (colored)
    geom_histogram(data = df_sample, aes(x = depth, y = ..count.. * scale_factor, fill = ref),
                   binwidth = 1, alpha = 0.8, position = "identity", color = NA) +
    
    # Threshold vertical line
    {if (!is.na(threshold_val)) geom_vline(xintercept = threshold_val, color = "red", linetype = "dashed")} +
    
    scale_y_continuous(
      name = "Whole-genome k-mer frequency",
      sec.axis = sec_axis(~ . / scale_factor, name = "HTF k-mer frequency")
    ) +
    scale_fill_manual(values = htf_colors, drop = FALSE) +
    labs(
      title = title,
      x = "K-mer Depth",
      fill = "HTF Reference"
    ) +
    theme_minimal() +
    theme(
      axis.title.y.right = element_text(color = "black"),
      axis.title.y.left = element_text(color = "grey30"),
      plot.title = element_text(hjust = 0.5)
    )
  
  ggsave(file.path(out_dir, paste0(sample_name, "_merged.pdf")), plot = p, width = 7, height = 5)
}
