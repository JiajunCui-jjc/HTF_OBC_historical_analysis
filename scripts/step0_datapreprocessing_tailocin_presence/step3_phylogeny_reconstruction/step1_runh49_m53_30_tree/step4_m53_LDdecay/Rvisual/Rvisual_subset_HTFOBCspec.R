#!/usr/bin/env Rscript
# ============================================================
# Genome-wide LD decay (short 0–2 kb, long 0–2 Mb) + LD50%
# Combined plot built from a unified dataframe
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

# ------------------------------
# Paths (edit if needed)
# ------------------------------
base_dir <- "/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step3_phylogeny_tree/results_m53_LDdecay/LDcalculation"
out_dir  <- file.path(base_dir, "LDvisual_full_and_HTF_OBC_fixed")
dir.create(out_dir, showWarnings = FALSE)

ld_short_file <- file.path(base_dir, "fullLDdecay_2kbwindow", "LD_2kb.hap.ld")
ld_long_file  <- file.path(base_dir, "subset2_LDdecay2mb_bg", "LD_subset_2pct_2mb.hap.ld")
out_prefix <- file.path(out_dir, "LDresults")

# HTF–OBC distance (bp) for vertical guide in the long panel
HTF <- c(2712429, 2714231)
OBC <- c(1012835, 1032652)
HTF_OBC_dist <- abs(mean(HTF) - mean(OBC))

# ------------------------------
# Helpers
# ------------------------------
read_ld_file <- function(path) {
  # autodetect 5- or 7+ column ld format
  dat <- suppressMessages(read_table(path, comment = "#", col_names = TRUE))
  if (ncol(dat) == 5) {
    colnames(dat) <- c("CHR","POS1","POS2","N_CHR","R2")
  } else if (ncol(dat) >= 7) {
    colnames(dat)[1:7] <- c("CHR","POS1","POS2","N_CHR","R2","D","Dprime")
  } else {
    stop("Unexpected columns in: ", path)
  }
  dat %>%
    mutate(
      POS1 = as.numeric(POS1),
      POS2 = as.numeric(POS2),
      R2   = as.numeric(R2),
      DIST = abs(POS2 - POS1)
    ) %>%
    filter(is.finite(DIST), is.finite(R2))
}

# bin by distance and compute median r^2
ld_decay <- function(df, max_dist, bin_bp) {
  brks <- seq(0, max_dist, by = bin_bp)
  df %>%
    filter(DIST <= max_dist) %>%
    mutate(bin = cut(DIST, breaks = brks, include.lowest = TRUE, right = TRUE)) %>%
    group_by(bin) %>%
    summarise(median_r2 = median(R2, na.rm = TRUE), .groups = "drop") %>%
    mutate(
  	mid = vapply(strsplit(gsub("\\[|\\]|\\(|\\)", "", bin), ","), function(x) {
    vals <- suppressWarnings(as.numeric(x))
    if (length(vals) == 2 && all(!is.na(vals))) mean(vals) else NA_real_
	  }, numeric(1))
	) %>%
	filter(!is.na(mid)) %>%
	arrange(mid)

}

# first distance where median r^2 <= half of max
ld50_first <- function(decay_df) {
  if (nrow(decay_df) == 0) return(NA_real_)
  r2_half <- max(decay_df$median_r2, na.rm = TRUE) * 0.5
  idx <- which(decay_df$median_r2 <= r2_half)[1]
  if (is.na(idx)) return(NA_real_)
  decay_df$mid[idx]
}

# Nice theme
theme_ld <- theme_bw(base_size = 8) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black")
  )

# ------------------------------
# Read data
# ------------------------------
message("→ Reading LD files...")
ld_short_raw <- read_ld_file(ld_short_file)   # 0–2 kb
ld_long_raw  <- read_ld_file(ld_long_file)    # 0–2 Mb (subset background)

# ------------------------------
# Build decays
# ------------------------------
# short: 0–2 kb in 5 bp bins (dense)
decay_short <- ld_decay(ld_short_raw, max_dist = 2000,  bin_bp = 5) %>%
  mutate(Scale = "0–2 kb")

# long: 0–2 Mb in 10 kb bins (smooth)
decay_long  <- ld_decay(ld_long_raw,  max_dist = 2e6, bin_bp = 10000) %>%
  mutate(Scale = "0–2 Mb")

# LD50% (first crossing)
ld50_short_x <- ld50_first(decay_short)
ld50_long_x  <- ld50_first(decay_long)

ld50_short_y <- if (is.finite(ld50_short_x)) max(decay_short$median_r2, na.rm=TRUE) * 0.5 else NA_real_
ld50_long_y  <- if (is.finite(ld50_long_x))  max(decay_long$median_r2,  na.rm=TRUE) * 0.5 else NA_real_

message(sprintf("→ LD50%% short at %.0f bp; long at %.0f bp",
                ld50_short_x, ld50_long_x))

# unify for combined plotting (shared y scale)
decay_all <- bind_rows(decay_short, decay_long)

# same y limits across all plots
y_max <- max(decay_all$median_r2, na.rm = TRUE) * 1.05
y_min <- 0

# ------------------------------
# Plot: short (0–2 kb)
# ------------------------------
p_short <- ggplot(decay_short, aes(x = mid, y = median_r2)) +
  geom_line(linewidth = 0.5) +
  geom_hline(yintercept = ld50_short_y, linetype = "dashed", color = "#de2d26", linewidth = 0.3) +
  geom_vline(xintercept = ld50_short_x, linetype = "dashed", color = "#de2d26", linewidth = 0.3) +
  annotate("text",
           x = min(max(ld50_short_x + 120, 300), 1800),
           y = ld50_short_y,
           label = sprintf("LD50%%  = %.0f bp", ld50_short_x),
           hjust = 0, vjust = -0.6, size = 2.3, color = "#de2d26") +
  scale_x_continuous(limits = c(0, 2000), breaks = seq(0, 2000, 500)) +
  scale_y_continuous(limits = c(y_min, y_max), labels = label_number(accuracy = 0.0005)) +
  labs(x = "Distance between SNPs (bp)", y = expression(median~r^2),
       title = "Genome-wide LD decay (0–2 kb)") +
  theme_ld

ggsave(file.path(out_dir, "LDdecay_short_0_2kb.pdf"),
       p_short, width = 180/25.4, height = 55/25.4, useDingbats = FALSE)

# ------------------------------
# Plot: long (0–2 Mb)
# ------------------------------
p_long <- ggplot(decay_long, aes(x = mid, y = median_r2)) +
  geom_line(linewidth = 0.4) +
  geom_hline(yintercept = ld50_long_y, linetype = "dashed", color = "#de2d26", linewidth = 0.3) +
  geom_vline(xintercept = ld50_long_x, linetype = "dashed", color = "#de2d26", linewidth = 0.3) +
  geom_vline(xintercept = HTF_OBC_dist, linetype = "dotted", color = "#3182bd", linewidth = 0.4) +
  annotate("text",
           x = min(HTF_OBC_dist + 3e4, 2e6 - 5e4),
           y = y_min + 0.75 * (y_max - y_min),
           label = sprintf("HTF–OBC ≈ %.0f kb", HTF_OBC_dist/1e3),
           hjust = 0, vjust = 0, size = 2.3, color = "#3182bd") +
  annotate("text",
           x = min(max(ld50_long_x + 1.5e5, 5e5), 1.8e6),
           y = ld50_long_y,
           label = sprintf("LD50%% = %.0f kb", ld50_long_x/1e3),
           hjust = 0, vjust = -0.6, size = 2.3, color = "#de2d26") +
  scale_x_continuous(limits = c(0, 2e6),
                     breaks = c(0, 5e5, 1e6, 1.5e6, 2e6),
                     labels = function(x) paste0(x/1e3, "k")) +
  scale_y_continuous(limits = c(y_min, y_max), labels = label_number(accuracy = 0.0005)) +
  labs(x = "Distance between SNPs (bp)", y = expression(median~r^2),
       title = "Genome-wide LD decay (0–2 Mb)") +
  theme_ld

ggsave(file.path(out_dir, "LDdecay_long_0_2Mb.pdf"),
       p_long, width = 180/25.4, height = 55/25.4, useDingbats = FALSE)

# ------------------------------
# Combined 2-panel plot from unified df (same y)
# ------------------------------
# We re-draw from the underlying data so both panels share aesthetics & y scale.
# (We still annotate LD50% and HTF–OBC on the correct panel.)

p_short2 <- ggplot(decay_all %>% filter(Scale == "0–2 kb"),
                   aes(x = mid, y = median_r2)) +
  geom_line(linewidth = 0.5) +
  geom_hline(yintercept = ld50_short_y, linetype = "dashed", color = "#de2d26", linewidth = 0.3) +
  geom_vline(xintercept = ld50_short_x, linetype = "dashed", color = "#de2d26", linewidth = 0.3) +
  annotate("text",
           x = min(max(ld50_short_x + 120, 300), 1800),
           y = ld50_short_y,
           label = sprintf("LD50%%  = %.0f bp", ld50_short_x),
           hjust = 0, vjust = -0.6, size = 2.2, color = "#de2d26") +
  scale_x_continuous(limits = c(0, 2000), breaks = seq(0, 2000, 500)) +
  scale_y_continuous(limits = c(y_min, y_max), labels = label_number(accuracy = 0.0005)) +
  labs(x = "Distance between SNPs (bp)", y = expression(median~r^2),
       title = "Genome-wide LD decay (0–2 kb)") +
  theme_ld

p_long2 <- ggplot(decay_all %>% filter(Scale == "0–2 Mb"),
                  aes(x = mid, y = median_r2)) +
  geom_line(linewidth = 0.4) +
  geom_hline(yintercept = ld50_long_y, linetype = "dashed", color = "#de2d26", linewidth = 0.3) +
  geom_vline(xintercept = ld50_long_x, linetype = "dashed", color = "#de2d26", linewidth = 0.3) +
  geom_vline(xintercept = HTF_OBC_dist, linetype = "dotted", color = "#3182bd", linewidth = 0.4) +
  annotate("text",
           x = min(HTF_OBC_dist + 3e4, 2e6 - 5e4),
           y = y_min + 0.75 * (y_max - y_min),
           label = sprintf("HTF–OBC ≈ %.0f kb", HTF_OBC_dist/1e3),
           hjust = 0, vjust = 0, size = 2.2, color = "#3182bd") +
  annotate("text",
           x = min(max(ld50_long_x + 1.5e5, 5e5), 1.8e6),
           y = ld50_long_y,
           label = sprintf("LD50%% = %.0f kb", ld50_long_x/1e3),
           hjust = 0, vjust = -0.6, size = 2.2, color = "#de2d26") +
  scale_x_continuous(limits = c(0, 2e6),
                     breaks = c(0, 5e5, 1e6, 1.5e6, 2e6),
                     labels = function(x) paste0(x/1e3, "k")) +
  scale_y_continuous(limits = c(y_min, y_max), labels = label_number(accuracy = 0.0005)) +
  labs(x = "Distance between SNPs (bp)", y = expression(median~r^2),
       title = "Genome-wide LD decay (0–2 Mb)") +
  theme_ld +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y  = element_blank())

p_combined <- p_short2 + p_long2 + plot_layout(widths = c(1, 1))
ggsave(file.path(out_dir, "LDdecay_combined_short_long.pdf"),
       p_combined, width = 180/25.4, height = 55/25.4, useDingbats = FALSE)

message("\n✔ Saved:")
message("  - ", file.path(out_dir, "LDdecay_short_0_2kb.pdf"))
message("  - ", file.path(out_dir, "LDdecay_long_0_2Mb.pdf"))
message("  - ", file.path(out_dir, "LDdecay_combined_short_long.pdf"))
