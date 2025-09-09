library(data.table)
library(ggplot2)
library(scales)
#  # === Filtered distribution: count >2 and ≤ mean+3SD ===

# === Paths ===
# === Paths ===
#dump_dir <- "/Users/jiajuncui/Desktop/others/tmphernan/2025_June/step9_coinfection/dump/tmp1"
dump_dir <- "/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step2_isolate_kmers/isolate_jf/dump"
#those are after dedup
out_dir <- file.path(dirname(dump_dir), "distribution")
out_dir_filtered <- file.path(out_dir, "filtered")
dir.create(out_dir, showWarnings = FALSE)
dir.create(out_dir_filtered, showWarnings = FALSE)

# === List all *_dump.txt files ===
dump_files <- list.files(dump_dir, pattern = "dump.txt$", full.names = TRUE)

# Assuming dump_files is a list of full paths
for (file in dump_files) {
	# file: process_one_dump.R fragment (drop-in replacement for your per-file code)
sample <- sub("_dump.txt$", "", basename(file))
df <- fread(file, header = FALSE, col.names = c("kmer", "count"))
if (nrow(df) == 0) next

# --- Use raw counts per k-mer (like awk), not just bins ---
counts_raw <- df$count
# Step 1: left-trim like awk: keep only > 2
counts_gt2 <- counts_raw[counts_raw > 2]
if (length(counts_gt2) == 0) {
  message(sample, ": no counts > 2, skipping")
  next
}

# Step 2: compute mean/sd on >2, then threshold
depth_mean <- mean(counts_gt2)
depth_sd   <- sd(counts_gt2)
threshold  <- depth_mean + 3 * depth_sd

# Step 3: filter by <= threshold (same domain as step 1)
counts_trimmed <- counts_gt2[counts_gt2 <= threshold]
if (length(counts_trimmed) == 0) {
  message(sample, ": filtered out (no counts <= mean+3sd). Stats: mean=", 
          signif(depth_mean,4), " sd=", signif(depth_sd,4), " thr=", signif(threshold,4))
  next
}

# Build frequency table for plotting (like your ggplot)
freq_df_filt <- as.data.table(table(counts_trimmed))
setnames(freq_df_filt, c("count", "frequency"))
freq_df_filt[, count := as.numeric(count)]

# Optional: also compute the "awk-like" final mean on trimmed set (rounded)
mean2 <- round(mean(counts_trimmed))
message(sample, ": n_gt2=", length(counts_gt2), 
        " mean=", signif(depth_mean,4), 
        " sd=", signif(depth_sd,4), 
        " thr=", signif(threshold,4),
        " mean_after_trim≈", mean2)

# Plot
p_filt <- ggplot(freq_df_filt, aes(x = count, y = frequency)) +
  geom_col(fill = "#0072B2", alpha = 0.95, width = 0.9) +
  labs(
    title = paste("Filtered (count >2, ≤ mean+3·sd):", sample),
    subtitle = sprintf("mean=%.2f, sd=%.2f, thr=%.2f; final mean≈%d",
                       depth_mean, depth_sd, threshold, mean2),
    x = "K-mer Depth",
    y = "Number of K-mers"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, colour = "grey30")
  )

ggsave(filename = file.path(out_dir_filtered, paste0(sample, "_filtered_kmer_depth_dist.pdf")),
       plot = p_filt, width = 7, height = 5)
}
cat("✅ Done: Plots saved to:\n  -", out_dir, "\n  -", out_dir_filtered, "\n") 
