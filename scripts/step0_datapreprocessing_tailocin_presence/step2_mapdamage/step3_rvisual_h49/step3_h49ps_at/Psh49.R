
#!/usr/bin/env Rscript
# -------------------------------------------------------------------------
# Plot 5′ C→T and 3′ G→A misincorporation frequencies for 49 historical samples
# Matching based on partial sample names (e.g., "109.NOR_1990" in full names)
# -------------------------------------------------------------------------

# === Libraries ===
library(ggplot2)
library(dplyr)
library(stringr)

# === Directories ===
indir <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step2_mapdamageplot/tops"
outdir <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step2_mapdamageplot/suppfig3mapdamage"

# === Input files ===
file_CtoT <- file.path(indir, "allinoneCtoT.txt")  # 5′ end
file_GtoA <- file.path(indir, "allinoneGtoA.txt")  # 3′ end
sample_file <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/data/historical49.txt"

# === Load mapDamage summary data ===
CtoT <- read.table(file_CtoT, header = TRUE, sep = "\t")
GtoA <- read.table(file_GtoA, header = TRUE, sep = "\t")

CtoT$End <- "5p"
GtoA$End <- "3p"
df_all <- rbind(CtoT, GtoA)

# === Load sample list (49 historical samples) ===
samples_49 <- read.table(sample_file, header = FALSE)[[1]]

# Combine into regex OR pattern for partial matching
pattern_49 <- paste(samples_49, collapse = "|")

# Filter rows where sample name contains any historical ID substring
df_49 <- df_all[str_detect(df_all$sample, pattern_49), ]

# === Summary of inclusion ===
matched_samples <- unique(df_49$sample)
matched_short <- samples_49[sapply(samples_49, function(x)
  any(str_detect(df_49$sample, fixed(x))))]

n_matched <- length(matched_short)
n_missing <- length(setdiff(samples_49, matched_short))

cat("\n✅ Matched", n_matched, "of 49 historical samples\n")
if (n_missing > 0) {
  cat("⚠️  Missing (not found in mapDamage outputs):\n  ",
      paste(setdiff(samples_49, matched_short), collapse = ", "), "\n")
}

# Summary per end
summary_counts <- df_49 %>%
  group_by(End) %>%
  summarise(
    n_samples = n_distinct(sample),
    n_points = n()
  )
print(summary_counts)

# Split for plotting
df_5p <- df_49 %>% filter(End == "5p")
df_3p <- df_49 %>% filter(End == "3p")

# === Plot: 5′ C→T (red) ===
p1 <- ggplot(df_5p, aes(x = pos_from_end, y = proportion, group = sample)) +
  geom_line(color = "#D04F4F", size = 0.8, alpha = 0.8) +
  labs(
    title = expression("5′ C to T deamination in 46 historical " * italic("Pseudomonas") * " sp. genomes"),
    x = "Distance from 5′ end (bp)",
    y = "C→T frequency"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  ylim(0, 0.04)

# === Plot: 3′ G→A (blue) ===
p2 <- ggplot(df_3p, aes(x = pos_from_end, y = proportion, group = sample)) +
  geom_line(color = "#3A7BDC", size = 0.8, alpha = 0.8) +
  labs(
    title = expression("3′ G to A deamination in 49 historical " * italic("Pseudomonas") * " sp. genomes"),
    x = "Distance from 3′ end (bp)",
    y = "G→A frequency"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  ylim(0, 0.04)

# === Save figures ===
ggsave(file.path(outdir, "CtoT_Ps.png"), p1, width = 9, height = 5, dpi = 600)
ggsave(file.path(outdir, "GtoA_Ps.png"), p2, width = 9, height = 5, dpi = 600)

cat("\n🎯 Saved plots:\n  -", file.path(outdir, "CtoT_At.png"),
    "\n  -", file.path(outdir, "GtoA_At.png"), "\n")
