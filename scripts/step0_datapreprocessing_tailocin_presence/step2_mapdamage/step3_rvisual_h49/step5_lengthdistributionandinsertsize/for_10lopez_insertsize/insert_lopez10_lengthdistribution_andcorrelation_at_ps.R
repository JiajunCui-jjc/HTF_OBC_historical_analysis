#!/usr/bin/env Rscript
# ============================================================
# Supp Fig — Fragment Length + Insert Size
# Nature-style Helvetica version (uniform font sizing)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(patchwork)
  library(ggpmisc)
  library(ggrepel)
})

# ---------------------------
# Global Helvetica theme
# ---------------------------
theme_set(
  theme_classic(base_size = 17, base_family = "Helvetica") +
    theme(
      axis.line   = element_line(color = "black", linewidth = 0.6),
      axis.ticks  = element_line(color = "black"),
      axis.text   = element_text(color = "black", size = 14, face = "bold"),
      axis.title  = element_text(size = 16, face = "bold"),
      plot.title  = element_text(size = 18, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 15, face = "bold"),
      legend.text  = element_text(size = 14),
      panel.grid  = element_blank()
    )
)

# ------------------------------------------------------------
# Input paths
# ------------------------------------------------------------
all_lengths_file <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step1_h49metadata/insert_sizes_loepz10/ALL_lengths.tsv"
tmreads_file     <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step1_h49metadata/insert_sizes_loepz10/tmreads_insert_merged.txt"

out_dir <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step2_mapdamageplot/step5_lengthdistribution_insertsize"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Fragment length distribution
# ============================================================
df_len <- read_tsv(all_lengths_file, show_col_types = FALSE) %>%
  group_by(Species, Sample, Length) %>%
  summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  group_by(Species, Sample) %>%
  mutate(Prop = Count / sum(Count, na.rm = TRUE)) %>%
  ungroup()

# ---- Arabidopsis ----
p_arabi_len <- ggplot(df_len %>% filter(Species == "Arabidopsis", Length > 0 & Length <= 500),
                      aes(x = Length, y = Prop, group = Sample)) +
  geom_line(alpha = 0.9, color = "#d5da6d") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.03)) +
  scale_x_continuous(limits = c(0, 500), breaks = seq(0, 500, 100)) +
  labs(x = "Fragment length (bp)", y = "Reads (%)",
       title = expression(italic("Arabidopsis thaliana")))

# ---- Pseudomonas ----
p_ps_len <- ggplot(df_len %>% filter(Species == "Pseudomonas", Length > 0 & Length <= 500),
                   aes(x = Length, y = Prop, group = Sample)) +
  geom_line(alpha = 0.9, color = "#e86b7d") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.03)) +
  scale_x_continuous(limits = c(0, 500), breaks = seq(0, 500, 100)) +
  labs(x = "Fragment length (bp)", y = "Reads (%)",
       title = expression(italic("Pseudomonas")~plain("sp.")))

ggsave(file.path(out_dir, "SuppFig_FragmentLength_Arabidopsis.pdf"), p_arabi_len, width = 7, height = 5, dpi = 600, useDingbats = FALSE)
ggsave(file.path(out_dir, "SuppFig_FragmentLength_Pseudomonas.pdf"), p_ps_len, width = 7, height = 5, dpi = 600, useDingbats = FALSE)

# ============================================================
# 2. Insert size summary (tmreads_insert_merged)
# ============================================================
df_tm <- read.table(tmreads_file) %>% distinct()
colnames(df_tm) <- c("Sample", "CollapsedReads", "R1Reads", "CollapsedProp", "MeanInsert_At", "MeanInsert_Ps")

df_tm <- df_tm %>%
  mutate(across(c(CollapsedReads, R1Reads), ~as.numeric(gsub(",", "", .x)))) %>%
  mutate(across(c(CollapsedProp, MeanInsert_At, MeanInsert_Ps), as.numeric)) %>%
  mutate(tmprop = CollapsedReads / R1Reads) %>%
  filter(is.finite(tmprop), !is.na(MeanInsert_At), !is.na(MeanInsert_Ps))

df_arabi_tm <- df_tm %>% mutate(Species = "Arabidopsis", MeanInsert = MeanInsert_At)
df_ps_tm    <- df_tm %>% mutate(Species = "Pseudomonas", MeanInsert = MeanInsert_Ps)

# ---- Insert-size function (with regression & R²/p) ----
make_insert_plot <- function(df, color, species_label){
  ggplot(df, aes(x = MeanInsert, y = tmprop)) +
    geom_point(color = color, alpha = 0.8, size = 3.2) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.8) +
    ylim(0, 1) +
    stat_poly_eq(
      aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
      formula = y ~ x, parse = TRUE, size = 4.5, vjust = 1.3, hjust = -0.8
    ) +
    labs(x = "Mean insert size (bp)", y = "Collapsed / R1 ratio",
         title = species_label)
}

p_arabi_tm <- make_insert_plot(df_arabi_tm, "#d5da6d", expression(italic("Arabidopsis thaliana")))
p_ps_tm    <- make_insert_plot(df_ps_tm, "#e86b7d", expression(italic("Pseudomonas")~plain("sp.")))

ggsave(file.path(out_dir, "SuppFig_InsertSize_Arabidopsis.pdf"), p_arabi_tm, width = 6, height = 4, dpi = 600, useDingbats = FALSE)
ggsave(file.path(out_dir, "SuppFig_InsertSize_Pseudomonas.pdf"), p_ps_tm, width = 6, height = 4, dpi = 600, useDingbats = FALSE)

# ============================================================
# 3. Combined 4-panel figure + red threshold line (y=0.65)
# ============================================================
p_arabi_tm_ref <- p_arabi_tm +
  geom_hline(yintercept = 0.65, color = "red", linetype = "dashed", linewidth = 0.6) +
  annotate("text", x = Inf, y = 0.65, label = "y = 0.65", hjust = 1.1, vjust = -0.5, color = "red", size = 4.2)

p_ps_tm_ref <- p_ps_tm +
  geom_hline(yintercept = 0.65, color = "red", linetype = "dashed", linewidth = 0.6) +
  annotate("text", x = Inf, y = 0.65, label = "y = 0.65", hjust = 1.1, vjust = -0.5, color = "red", size = 4.2)

p_combined <- (p_arabi_len + p_ps_len) / (p_arabi_tm_ref + p_ps_tm_ref) +
  plot_annotation(
    tag_levels = 'A',
    theme = theme(plot.tag = element_text(size = 18, face = "bold", family = "Helvetica"))
  )

ggsave(file.path(out_dir, "SuppFig_FragmentLength_InsertSize_Combined_withRefLine.pdf"),
       plot = p_combined, width = 10, height = 8, dpi = 600, useDingbats = FALSE)

# ============================================================
# 4. Fragment length colored and labeled by isolate
# ============================================================
df_labels <- df_len %>%
  filter(Length > 0) %>%
  group_by(Species, Sample) %>%
  slice_max(Length, n = 1, with_ties = FALSE) %>%
  ungroup()

make_len_labeled <- function(species_name, color_hex) {
  ggplot(df_len %>% filter(Species == species_name, Length > 0),
         aes(x = Length, y = Prop, group = Sample, color = Sample)) +
    geom_line(linewidth = 0.8, alpha = 0.85) +
    geom_text_repel(
      data = df_labels %>% filter(Species == species_name),
      aes(label = Sample),
      size = 4.5, hjust = 0, nudge_x = 20,
      segment.size = 0.25, show.legend = FALSE
    ) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.03)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = "Fragment length (bp)", y = "Reads (%)",
         title = if (species_name == "Arabidopsis")
           expression(italic("Arabidopsis thaliana")~" fragment length per isolate")
         else
           expression(italic("Pseudomonas")~plain("sp.")~" fragment length per isolate"))
}

p_arabi_len_labeled <- make_len_labeled("Arabidopsis", "#d5da6d")
p_ps_len_labeled    <- make_len_labeled("Pseudomonas", "#e86b7d")

ggsave(file.path(out_dir, "SuppFig_FragmentLength_Arabidopsis_ByIsolate_Labeled.pdf"),
       p_arabi_len_labeled, width = 8, height = 5, dpi = 600, useDingbats = FALSE)
ggsave(file.path(out_dir, "SuppFig_FragmentLength_Pseudomonas_ByIsolate_Labeled.pdf"),
       p_ps_len_labeled, width = 8, height = 5, dpi = 600, useDingbats = FALSE)

# ============================================================
cat("✅ All SuppFig panels exported with Helvetica font (Nature-style) to:\n", out_dir, "\n")
