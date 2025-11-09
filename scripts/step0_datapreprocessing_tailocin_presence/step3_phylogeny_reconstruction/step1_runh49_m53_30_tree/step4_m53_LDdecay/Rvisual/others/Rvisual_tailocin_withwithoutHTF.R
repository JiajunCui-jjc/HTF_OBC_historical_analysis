#!/usr/bin/env Rscript
# ============================================================
# LD decay within tailocin region — with vs without HTF locus
# Each curve labeled with its own LD₅₀%
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
})

# ------------------------------------------------------------
# Input paths
# ------------------------------------------------------------
ld_path <- "../../../../../../results/step0_datapreprocessing_tailocin_presence/step3_phylogeny_tree/results_m53_LDdecay/LDcalculation/tailocin_only_region/LD_tailocin_20kb.hap.ld"
out_pdf <- gsub(".hap.ld$", "_with_vs_without_HTF_LD50.pdf", ld_path)

# HTF coordinates
HTF <- c(2712429, 2714231)

# ------------------------------------------------------------
# Read LD file
# ------------------------------------------------------------
ld <- suppressMessages(read_table(ld_path, comment = "#", col_names = TRUE))
colnames(ld)[1:7] <- c("CHR","POS1","POS2","N_CHR","R2","D","Dprime")

ld <- ld %>%
  mutate(DIST = abs(POS2 - POS1)) %>%
  filter(!is.na(R2), !is.na(DIST))

# ------------------------------------------------------------
# Split into "with HTF" vs "without HTF"
# ------------------------------------------------------------
ld <- ld %>%
  mutate(Group = ifelse(
    (POS1 >= HTF[1] & POS1 <= HTF[2]) | (POS2 >= HTF[1] & POS2 <= HTF[2]),
    "With HTF", "Without HTF"
  ))

# ------------------------------------------------------------
# Summarise median r² by distance bins
# ------------------------------------------------------------
bin_size <- 100
ld_summary <- ld %>%
  mutate(bin = cut(DIST, breaks = seq(0, max(DIST, na.rm=TRUE)+bin_size, by = bin_size))) %>%
  group_by(Group, bin) %>%
  summarise(median_r2 = median(R2, na.rm = TRUE), .groups = "drop") %>%
  mutate(mid = as.numeric(sub("\\((.*),(.*)]", "\\1", bin)) + bin_size / 2)

# ------------------------------------------------------------
# Compute LD50% (half max r²) for each curve
# ------------------------------------------------------------
ld50_table <- ld_summary %>%
  group_by(Group) %>%
  summarise(
    max_r2 = max(median_r2, na.rm = TRUE),
    r2_half = max_r2 * 0.5,
    x_half = approx(median_r2, mid, xout = max_r2 * 0.5, ties = "ordered")$y,
    .groups = "drop"
  )

print(ld50_table)

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
p <- ggplot(ld_summary, aes(x = mid, y = median_r2, color = Group)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 0.6, alpha = 0.7) +
  # LD50 lines for each group
  geom_vline(data = ld50_table, aes(xintercept = x_half, color = Group),
             linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +
  geom_hline(data = ld50_table, aes(yintercept = r2_half, color = Group),
             linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +
  # LD50 labels
  geom_text(data = ld50_table,
            aes(x = x_half * 1.05, y = r2_half,
                label = sprintf("LD50%% ≈ %.0f bp\nr² = %.2f", x_half, r2_half),
                color = Group),
            hjust = 0, vjust = -0.4, size = 2.5, show.legend = FALSE) +
  scale_color_manual(values = c("With HTF" = "#de2d26", "Without HTF" = "#3182bd")) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(limits = c(0, 0.05), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Distance between SNPs (bp)",
       y = expression(median~r^2),
       title = "LD decay within tailocin region (21 kb length)",
       color = NULL) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2, color = "grey85"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    legend.position = "top",
    legend.text = element_text(face = "bold")
  )

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
ggsave(out_pdf, plot = p, width = 100/25.4, height = 70/25.4, useDingbats = FALSE)
message("✅ Saved: ", out_pdf)




# ============================================================
# (NEW) r² frequency distribution (count / total per group)
# ============================================================
message("→ Plotting r² frequency distribution by group ...")

# ensure ld still exists, otherwise reload as before
if (!exists("ld")) {
  ld <- suppressMessages(read_table(ld_path, comment = "#", col_names = TRUE))
  colnames(ld)[1:7] <- c("CHR","POS1","POS2","N_CHR","R2","D","Dprime")
  ld <- ld %>%
    mutate(DIST = abs(POS2 - POS1)) %>%
    filter(!is.na(R2), !is.na(DIST)) %>%
    mutate(Group = ifelse(
      (POS1 >= HTF[1] & POS1 <= HTF[2]) | (POS2 >= HTF[1] & POS2 <= HTF[2]),
      "With HTF", "Without HTF"
    ))
}

# ------------------------------------------------------------
# Bin R² values and compute proportion per group (robust version)

# ============================================================
# (FINAL) Histogram of r² — proportional counts (0–0.05 range)
# ============================================================
message("→ Plotting refined histogram of r² (proportional counts) ...")

# Filter to relevant range
ld_filtered <- ld %>%
  filter(R2 >= 0, R2 <= 0.05)

# Choose appropriate bin size dynamically based on data density
n_points <- nrow(ld_filtered)
bin_size <- if (n_points > 5e5) 0.0005 else if (n_points > 1e5) 0.001 else 0.002

message(sprintf("Using bin size = %.4f for %d SNP pairs", bin_size, n_points))

# Prepare histogram data: counts → proportions per group
ld_hist <- ld_filtered %>%
  mutate(bin = cut(R2, breaks = seq(0, 0.05, by = bin_size), include.lowest = TRUE)) %>%
  group_by(Group, bin) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  mutate(total = sum(count), prop = count / total) %>%
  ungroup() %>%
  mutate(
    mid = vapply(strsplit(gsub("\\[|\\]|\\(|\\)", "", bin), ","), function(x) {
      if (length(x) == 2 && all(!is.na(as.numeric(x)))) mean(as.numeric(x)) else NA_real_
    }, numeric(1))
  ) %>%
  filter(!is.na(mid))

# ============================================================
# Plot
# ============================================================
p_hist <- ggplot(ld_hist, aes(x = mid, y = prop, fill = Group)) +
  geom_col(position = "identity", alpha = 0.5, width = bin_size) +
  scale_fill_manual(values = c("With HTF" = "#de2d26", "Without HTF" = "#3182bd")) +
  scale_x_continuous(
    limits = c(0, 0.06),
    breaks = seq(0, 0.06, 0.01),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = expression(r^2),
    y = "Proportion of total pairs (%)",
    title = expression("r"^2~"distribution within tailocin region (r2 less than 0.05)"),
    fill = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    legend.position = "top",
    legend.text = element_text(face = "bold")
  )

# ============================================================
# Save
# ============================================================
# Save
# ------------------------------------------------------------
out_pdf_r2dist <- gsub(".hap.ld$", "_r2_frequency_distribution_with_vs_without_HTF.pdf", ld_path)
ggsave(out_pdf_r2dist, plot = p_hist, width = 110/25.4, height = 70/25.4, useDingbats = FALSE)
message("✅ Saved r² frequency distribution → ", out_pdf_r2dist)# ------------------------------------------------------------
