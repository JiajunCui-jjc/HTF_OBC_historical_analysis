#!/usr/bin/env Rscript
# -------------------------------------------------------------------------
# Summarize and visualize read length distributions (insertsize)
# for Arabidopsis thaliana and Pseudomonas viridiflava
# using 2025_summerpaper_all mapDamage outputs.
# -------------------------------------------------------------------------

# === Libraries ===
library(purrr)
library(dplyr)
library(ggplot2)
library(readr)
library(scales)
library(patchwork)

# === Function to read lgdistribution or length_distribution.txt ===
read_lgdist <- function(file, sample, genome) {
  dat <- suppressWarnings(read_tsv(file, comment = "#", show_col_types = FALSE))
  
  # Handle variable column naming
  if (all(c("Std", "Length", "Occurences") %in% colnames(dat))) {
    dat <- dat %>% select(Length, Occurences)
  } else if (all(c("Length", "Occurences") %in% colnames(dat))) {
    dat <- dat %>% select(Length, Occurences)
  } else {
    dat <- dat[, (ncol(dat)-1):ncol(dat)]
    colnames(dat) <- c("Length", "Occurences")
  }
  
  dat <- dat %>%
    mutate(
      Length = suppressWarnings(as.numeric(Length)),
      Occurences = suppressWarnings(as.numeric(Occurences))
    ) %>%
    group_by(Length) %>%
    summarise(Count = sum(Occurences, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      Sample = sample,
      Genome = genome
    )
  
  return(dat)
}

# === Directories ===
at_dirs <- c(
  "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step2_mapdamageplot/toAt"
)

ps_dirs <- c(
  "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step2_mapdamageplot/tops"
)

# === Collect files ===
at_files <- unlist(lapply(at_dirs, function(d)
  list.files(d, pattern = "(lgdistribution|length_distribution)\\.txt$", recursive = TRUE, full.names = TRUE)
))
ps_files <- unlist(lapply(ps_dirs, function(d)
  list.files(d, pattern = "(lgdistribution|length_distribution)\\.txt$", recursive = TRUE, full.names = TRUE)
))

# === Read and combine ===
at_all <- map_df(at_files, ~read_lgdist(.x, basename(dirname(.x)), "Arabidopsis thaliana"))
ps_all <- map_df(ps_files, ~read_lgdist(.x, basename(dirname(.x)), "Pseudomonas viridiflava"))

# === Keep only historical samples (PL and HB) ===
at_all <- at_all %>% filter(grepl("^(PL|HB)", Sample))
ps_all <- ps_all %>% filter(grepl("^(PL|HB)", Sample))

cat("✅ Loaded:", length(unique(at_all$Sample)), "Arabidopsis samples and",
    length(unique(ps_all$Sample)), "Pseudomonas samples.\n")

# === Convert counts to proportions per sample ===
at_all <- at_all %>%
  group_by(Sample) %>%
  mutate(Prop = Count / sum(Count, na.rm = TRUE)) %>%
  ungroup()

ps_all <- ps_all %>%
  group_by(Sample) %>%
  mutate(Prop = Count / sum(Count, na.rm = TRUE)) %>%
  ungroup()

# === Output directory ===
out_dir <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step2_mapdamageplot/step5_lengthdistribution_insertsize"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# === Plot: Arabidopsis ===
p_at <- ggplot(at_all %>% filter(Length > 0 & Length <= 200),
               aes(x = Length, y = Prop, group = Sample)) +
  geom_line(alpha = 0.8, color = "#d5da6d") +  # olive-green
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.15)) +
  scale_x_continuous(limits = c(0, 200), breaks = seq(0, 200, 50)) +
  labs(x = "Read length (bp)", y = "Reads (%)",
       title = expression(italic("Arabidopsis thaliana"))) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 11, face = "bold")
  )

# === Plot: Pseudomonas ===
p_ps <- ggplot(ps_all %>% filter(Length > 0 & Length <= 200),
               aes(x = Length, y = Prop, group = Sample)) +
  geom_line(alpha = 0.8, color = "#e86b7d") +  # coral
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.15)) +
  scale_x_continuous(limits = c(0, 200), breaks = seq(0, 200, 50)) +
  labs(x = "Read length (bp)", y = "Reads (%)",
       title = expression(italic(Pseudomonas) ~ "sp.")) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 11, face = "bold")
  )

# === Save separate panels ===
ggsave(file.path(out_dir, "ReadLength_Arabidopsis.pdf"), p_at, width = 7, height = 5)
ggsave(file.path(out_dir, "ReadLength_Pseudomonas.pdf"), p_ps, width = 7, height = 5)

# === Combine panels (Supp Fig style) ===
p_combined <- p_at + p_ps +
  plot_layout(ncol = 2) +
  plot_annotation(
    theme = theme(plot.tag = element_text(size = 14, face = "bold"))
  )

ggsave(file.path(out_dir, "SuppFig_ReadLength_AT_PS.pdf"),
       plot = p_combined, width = 12, height = 5, dpi = 600)

# === Weighted mean read length per sample ===
at_meanlen <- at_all %>%
  group_by(Sample) %>%
  summarise(MeanLength = weighted.mean(Length, Count, na.rm = TRUE),
            Genome = unique(Genome)) %>%
  arrange(MeanLength)

ps_meanlen <- ps_all %>%
  group_by(Sample) %>%
  summarise(MeanLength = weighted.mean(Length, Count, na.rm = TRUE),
            Genome = unique(Genome)) %>%
  arrange(MeanLength)

meanlen_all <- bind_rows(at_meanlen, ps_meanlen)

# === Save per-sample mean lengths ===
write.table(meanlen_all, file.path(out_dir, "MeanReadLength_perSample.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# === Summary by sample type (PL vs HB) ===
avg_summary <- meanlen_all %>%
  mutate(
    Type = case_when(
      grepl("^PL", Sample) ~ "PL",
      grepl("^HB", Sample) ~ "HB",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(Type, Genome) %>%
  summarise(
    AvgMeanLength = mean(MeanLength, na.rm = TRUE),
    N_Samples = n(),
    .groups = "drop"
  ) %>%
  arrange(Type, Genome)

write.table(avg_summary, file.path(out_dir, "AvgMeanReadLength_byType.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n✅ Done — Saved all plots and tables to:\n", out_dir, "\n")
print(avg_summary)
