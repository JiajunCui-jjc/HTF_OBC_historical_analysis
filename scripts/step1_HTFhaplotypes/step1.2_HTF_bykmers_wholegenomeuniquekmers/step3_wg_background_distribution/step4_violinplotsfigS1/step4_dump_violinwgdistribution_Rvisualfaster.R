# ============================================================
# Supplementary Figure — Normalized violin plot of h40 k-mer distributions
# Each isolate scaled by its own mean depth (mean = 1)
# ============================================================

library(data.table)
library(ggplot2)
library(stringr)

# ---- Paths ----
dump_dir <- "/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step2_isolate_kmers/isolate_jf/dump"
out_dir  <- file.path(dirname(dump_dir), "distribution", "summary_violin_norm1")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- List & filter ----
dump_files <- list.files(dump_dir, pattern = "dump.txt$", full.names = TRUE)
dump_files <- dump_files[!grepl("^p", basename(dump_files))]  # only non-p (historical)
cat("✅ Found", length(dump_files), "historical (h409 dump files\n")

# ---- Initialize ----
all_data <- vector("list", length(dump_files))

# ---- Process efficiently ----
for (i in seq_along(dump_files)) {
  file <- dump_files[i]
  sample <- sub("_dump.txt$", "", basename(file))

  df <- fread(file, select = 2, col.names = "count", showProgress = FALSE)
  if (nrow(df) == 0) next

  counts_gt2 <- df$count[df$count > 2]
  rm(df); gc(FALSE)
  if (length(counts_gt2) == 0) next

  depth_mean <- mean(counts_gt2)
  depth_sd   <- sd(counts_gt2)
  threshold  <- depth_mean + 3 * depth_sd
  counts_trimmed <- counts_gt2[counts_gt2 <= threshold]
  rm(counts_gt2); gc(FALSE)
  if (length(counts_trimmed) == 0) next

  # Normalize by mean (so mean depth = 1)
  depth_rel <- counts_trimmed / depth_mean

  all_data[[i]] <- data.table(Sample = sample, DepthRel = depth_rel)

  message(sprintf("[%02d/%02d] %s: mean=%.2f sd=%.2f thr=%.2f n=%d",
                  i, length(dump_files), sample, depth_mean, depth_sd, threshold, length(depth_rel)))

  rm(counts_trimmed, depth_rel); gc(FALSE)
}

# ---- Combine ----
df_all <- rbindlist(all_data, use.names = TRUE)
rm(all_data); gc(FALSE)
cat("✅ Combined normalized data for", length(unique(df_all$Sample)), "samples\n")

# ---- Precompute sample order (median normalized depth) ----
sample_order <- df_all[, .(med = median(DepthRel)), by = Sample][order(med)]$Sample
df_all[, Sample := factor(Sample, levels = sample_order)]

# ---- Violin Plot ----
p_violin <- ggplot(df_all, aes(x = Sample, y = DepthRel)) +
  geom_violin(
    fill = "#e86b7d",
    color = "black",
    scale = "width",
    trim = TRUE,
    linewidth = 0.25,
    bw = 0.3,
    adjust = 1.0,
    alpha = 0.9
  ) +
  geom_boxplot(
    width = 0.08,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.3
  ) +
  geom_hline(yintercept = 1, color = "red", linewidth = 0.4, linetype = "dashed") +
  labs(
    title = "Normalized k-mer depth distribution across historical ATUE5 isolates",
    x = "Sample",
    y = "Relative depth (mean-normalized)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12,face = "bold"),
    axis.text.y = element_text(size = 16,face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, colour = "grey40"),
    axis.title = element_text(size=20,face = "bold")
  )

# ---- Save ----
pdf_file <- file.path(out_dir, "SuppFig_h49_normalizedMean1_violin_kmer_depth.pdf")
png_file <- file.path(out_dir, "SuppFig_h49_normalizedMean1_violin_kmer_depth.png")

ggsave(pdf_file, p_violin, width = 12, height = 6)
ggsave(png_file, p_violin, width = 12, height = 6, dpi = 500)

cat("✅ Saved normalized (mean=1) violin summary figure:\n",
    " -", pdf_file, "\n",
    " -", png_file, "\n")
