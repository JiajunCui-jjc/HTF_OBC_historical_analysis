#!/usr/bin/env Rscript

# =============================================================
# Compare HTF k-mer depth with whole-genome (WG) k-mer depth.
# Modern samples → 1 red dashed line (mean cutoff)
# Historical samples → 1 red dashed lines (mean -  0.5SD)
# =============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(dplyr)
})

# --- Paths ---
htf_data_path <- "/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers/kmer_distribution_mix/HTF_kmer_depth_distributions.tsv"
wg_dir <- "/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step2_isolate_kmers/isolate_jf/dump"
out_dir <- "/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers/kmer_distribution_mix/distributions"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Thresholds ---
threshold_modern <- fread(file.path(wg_dir, "threshold_m57.tsv"), col.names = c("sample", "cutoff"))
threshold_historical <- fread(file.path(wg_dir, "threshold_h49_sd_twosides_0.5sd.tsv"), col.names = c("sample", "mean_minus_sd", "mean_plus_sd"))

threshold_map_modern <- setNames(threshold_modern$cutoff, threshold_modern$sample)
threshold_map_historical_low <- setNames(threshold_historical$mean_minus_sd, threshold_historical$sample)
threshold_map_historical_high <- setNames(threshold_historical$mean_plus_sd, threshold_historical$sample)

# --- Colors ---
htf_colors <- c(
  "HTF_p21.F9" = "#f9d42a",  # yellow (1245 bp)
  "HTF_p23.B8" = "#A30026",  # crimson red
  "HTF_p25.A12" = "#638ccc", # blue (1383 bp)
  "HTF_p25.C2" = "#C51B7D",  # magenta pink
  "HTF_p26.D6" = "#B2182B",  # scarlet red
  "HTF_p5.D5"  = "#D6604D",  # coral red
  "HTF_p7.G11" = "#f6d6ff"   # lilac (1830 bp)
)

# --- Get args ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("❌ Please provide sample name as argument, e.g. Rscript plot_HTF_WG_kmer_depth.R PL0001")
}
sample_name <- args[1]
cat("🔧 Processing:", sample_name, "\n")

# --- Load HTF data ---
df_depths <- fread(htf_data_path, col.names = c("sample", "ref", "kmer", "depth"))
df_sample <- df_depths[sample == sample_name & depth >= 1]
if (nrow(df_sample) == 0) {
  cat("⚠️  No HTF data for", sample_name, "\n")
  quit(status = 0)
}

# --- Load WG dump ---
# --- Load WG dump ---
wg_path1 <- file.path(wg_dir, paste0(sample_name, "_dump.txt"))
wg_path2 <- file.path(wg_dir, paste0(sample_name, "_dump_20p_dump.txt"))

if (file.exists(wg_path1)) {
  wg_path <- wg_path1
} else if (file.exists(wg_path2)) {
  wg_path <- wg_path2
} else {
  cat("❌ No WG file for", sample_name, "\n")
  quit(status = 0)
}

df_wg <- fread(wg_path, col.names = c("kmer", "count"))
df_wg <- df_wg[count > 2]
if (nrow(df_wg) == 0) {
  cat("❌ WG file has no informative k-mers for", sample_name, "\n")
  quit(status = 0)
}


# --- Filter WG to remove high-count tail ---
depth_mean <- mean(df_wg$count)
depth_sd   <- sd(df_wg$count)
cutoff_max <- depth_mean + 3 * depth_sd
df_wg <- df_wg[count <= cutoff_max]
if (nrow(df_wg) == 0) quit(status = 0)

freq_df_wg <- df_wg[, .N, by = count]
setnames(freq_df_wg, c("depth", "frequency"))

# --- Thresholds ---
is_modern <- grepl("^p", sample_name)
is_historical <- !is_modern

if (is_modern && sample_name %in% names(threshold_map_modern)) {
  threshold_val <- threshold_map_modern[[sample_name]]
} else if (is_historical && sample_name %in% names(threshold_map_historical_low)) {
  threshold_low <- threshold_map_historical_low[[sample_name]]
  threshold_high <- threshold_map_historical_high[[sample_name]]
} else {
  cat("❗ No threshold found for", sample_name, "\n")
  quit(status = 0)
}

# --- Scaling for secondary Y axis ---
max_wg_y <- max(freq_df_wg$frequency)
max_htf_y <- max(table(df_sample$depth))
scale_factor <- max_wg_y / max_htf_y

# --- Plot ---
p <- ggplot() +
  geom_col(data = freq_df_wg, aes(x = depth, y = frequency),
           fill = "grey50", alpha = 0.3) +
  geom_histogram(data = df_sample, aes(x = depth, y = ..count.. * scale_factor, fill = ref),
                 binwidth = 1, alpha = 0.8, position = "identity", color = NA) +
  scale_y_continuous(
    name = "Whole-genome k-mer frequency",
    sec.axis = sec_axis(~ . / scale_factor, name = "HTF k-mer frequency")
  ) +
  scale_fill_manual(values = htf_colors, drop = FALSE) +
  labs(title = paste0("HTF & WG k-mer depth: ", sample_name),
       x = "K-mer depth", fill = "HTF reference") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.y.right = element_text(color = "black"),
    axis.title.y.left = element_text(color = "grey30")
  )

# --- Add vertical threshold lines ---
if (is_modern) {
  p <- p + geom_vline(xintercept = threshold_val, color = "red", linetype = "dashed", linewidth = 0.8)
} else {
  p <- p +
    geom_vline(xintercept = threshold_low, color = "red", linetype = "dashed", linewidth = 0.8) }

# --- Save ---
outfile <- file.path(out_dir, paste0(sample_name, "_merged.pdf"))
ggsave(outfile, plot = p, width = 7, height = 5)
cat("✅ Saved:", outfile, "\n")


