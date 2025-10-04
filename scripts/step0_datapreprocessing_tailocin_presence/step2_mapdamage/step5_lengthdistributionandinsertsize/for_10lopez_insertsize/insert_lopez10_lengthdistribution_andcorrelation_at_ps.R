# ===============================
# Supp Fig: Fragment Length + Insert Size
# ===============================

library(dplyr)
library(ggplot2)
library(readr)
library(scales)
library(patchwork)

# === Input paths ===
all_lengths_file <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step2_mapdamageplot/step5_lengthdistribution_insertsize/suppfig2_insertsize_10/insert_sizes/ALL_lengths.tsv"
tmreads_file     <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step2_mapdamageplot/step5_lengthdistribution_insertsize/suppfig2_insertsize_10/insert_sizes/tmreads_insert_merged.txt"

# === Output directory ===
out_dir <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step2_mapdamageplot/step5_lengthdistribution_insertsize/suppfig2_insertsize_10/"

# ===============================
# 1. Fragment length distribution
# ===============================
df_len <- read_tsv(all_lengths_file, show_col_types = FALSE)

# combine collapsed + paired, then normalize to proportions
df_len <- df_len %>%
  group_by(Species, Sample, Length) %>%
  summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  group_by(Species, Sample) %>%
  mutate(Prop = Count / sum(Count, na.rm = TRUE)) %>%
  ungroup()

# ---- Arabidopsis ----
p_arabi_len <- ggplot(df_len %>% filter(Species == "Arabidopsis", Length > 0),
                      aes(x = Length, y = Prop, group = Sample)) +
  geom_line(alpha = 0.6, color = "#4575b4") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.03)) +
  labs(x = "Fragment length (bp)", y = "Reads (%)",
       title = expression(italic("Arabidopsis thaliana"))) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11, face = "bold")
  )

# ---- Pseudomonas ----
p_ps_len <- ggplot(df_len %>% filter(Species == "Pseudomonas", Length > 0),
                   aes(x = Length, y = Prop, group = Sample)) +
  geom_line(alpha = 0.6, color = "#d73027") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.03)) +
  labs(x = "Fragment length (bp)", y = "Reads (%)",
       title = expression(italic("Pseudomonas") ~ plain("sp."))) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11, face = "bold")
  )

# Save
ggsave(file.path(out_dir, "SuppFig_FragmentLength_Arabidopsis.pdf"),
       p_arabi_len, width = 7, height = 5)
ggsave(file.path(out_dir, "SuppFig_FragmentLength_Pseudomonas.pdf"),
       p_ps_len, width = 7, height = 5)


# ===============================
# 2. Insert size summary (tmreads_insert_merged)
# ===============================# ===============================
# 2. Insert size summary (tmreads_insert_merged)
# ===============================
df_tm <- read_tsv(tmreads_file, show_col_types = FALSE) %>%
  distinct()  # remove duplicate rows

colnames(df_tm) <- c("Sample", "CollapsedReads", "R1Reads", "CollapsedProp",
                     "MeanInsert_At", "MeanInsert_Ps")

# convert to numeric safely
df_tm <- df_tm %>%
  mutate(
    CollapsedReads = as.numeric(gsub(",", "", CollapsedReads)),
    R1Reads        = as.numeric(gsub(",", "", R1Reads)),
    CollapsedProp  = as.numeric(CollapsedProp),
    MeanInsert_At  = as.numeric(MeanInsert_At),
    MeanInsert_Ps  = as.numeric(MeanInsert_Ps)
  )
# clean NAs or invalid rows before plotting
df_tm <- df_tm %>%
  mutate(tmprop = CollapsedReads / R1Reads) %>%
  filter(is.finite(tmprop), !is.na(MeanInsert_At), !is.na(MeanInsert_Ps))
# calculate tmprop
df_tm <- df_tm %>%
  mutate(tmprop = CollapsedReads / R1Reads)

# split by species
df_arabi_tm <- df_tm %>% mutate(Species = "Arabidopsis", MeanInsert = MeanInsert_At)
df_ps_tm    <- df_tm %>% mutate(Species = "Pseudomonas", MeanInsert = MeanInsert_Ps)

# ---- Arabidopsis insert size ----
p_arabi_tm <- ggplot(df_arabi_tm, aes(x = MeanInsert, y = tmprop)) +
  geom_point(color = "#4575b4", alpha = 0.8, size = 2) +
  
  ylim(0, 1)+
  labs(x = "Mean insert size (bp)", y = "Collapsed / R1 ratio",
       title = expression(italic("Arabidopsis thaliana"))) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11, face = "bold")
  )

# ---- Pseudomonas insert size ----
p_ps_tm <- ggplot(df_ps_tm, aes(x = MeanInsert, y = tmprop)) +
  geom_point(color = "#d73027", alpha = 0.8, size = 2) +
  
  ylim(0, 1)+
  labs(x = "Mean insert size (bp)", y = "Collapsed / R1 ratio",
       title = expression(italic("Pseudomonas") ~ plain("sp."))) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11, face = "bold")
  )

# Save
ggsave(file.path(out_dir, "SuppFig_InsertSize_Arabidopsis.pdf"),
       p_arabi_tm, width = 6, height = 4)
ggsave(file.path(out_dir, "SuppFig_InsertSize_Pseudomonas.pdf"),
       p_ps_tm, width = 6, height = 4)

# ===============================
# Optional combined figure
# ===============================
p_combined <- (p_arabi_len + p_ps_len) / (p_arabi_tm + p_ps_tm) 
#+
#  plot_annotation(tag_levels = 'A')

ggsave(file.path(out_dir, "SuppFig2_InsertSize_FragmentLength_Combined.pdf"),
       p_combined, width = 10, height = 8, dpi = 600)




#with pvalue
# ===============================
# 2. Insert size summary (tmreads_insert_merged)
# ===============================
library(ggpmisc)   # for stat_poly_eq

df_tm <- read_tsv(tmreads_file, show_col_types = FALSE) %>%
  distinct()  # remove duplicate rows

colnames(df_tm) <- c("Sample", "CollapsedReads", "R1Reads", "CollapsedProp",
                     "MeanInsert_At", "MeanInsert_Ps")

# convert to numeric safely
df_tm <- df_tm %>%
  mutate(
    CollapsedReads = as.numeric(gsub(",", "", CollapsedReads)),
    R1Reads        = as.numeric(gsub(",", "", R1Reads)),
    CollapsedProp  = as.numeric(CollapsedProp),
    MeanInsert_At  = as.numeric(MeanInsert_At),
    MeanInsert_Ps  = as.numeric(MeanInsert_Ps)
  ) %>%
  mutate(tmprop = CollapsedReads / R1Reads) %>%
  filter(is.finite(tmprop), !is.na(MeanInsert_At), !is.na(MeanInsert_Ps))

# split by species
df_arabi_tm <- df_tm %>% mutate(Species = "Arabidopsis", MeanInsert = MeanInsert_At)
df_ps_tm    <- df_tm %>% mutate(Species = "Pseudomonas", MeanInsert = MeanInsert_Ps)

# ==== Helper function to build regression annotation ====
make_insert_plot <- function(df, color, species_label){
  ggplot(df, aes(x = MeanInsert, y = tmprop)) +
    geom_point(color = color, alpha = 0.8, size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
    ylim(0, 1) +
    stat_poly_eq(
      aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
      formula = y ~ x,
      parse = TRUE,
      vjust =1.5, hjust = -1,
      size = 3.5
    ) +
    labs(x = "Mean insert size (bp)", y = "Collapsed / R1 ratio",
         title = species_label) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 11, face = "bold")
    )
}

# ==== Build plots ====
p_arabi_tm <- make_insert_plot(df_arabi_tm,
                               "#4575b4",
                               expression(italic("Arabidopsis thaliana")))

p_ps_tm <- make_insert_plot(df_ps_tm,
                            "#d73027",
                              expression(italic("Pseudomonas") ~ plain("sp.")))

# ==== Save ====
ggsave(file.path(out_dir, "SuppFig_withp_InsertSize_Arabidopsis.pdf"),
       p_arabi_tm, width = 6, height = 4)
ggsave(file.path(out_dir, "SuppFig_withP_InsertSize_Pseudomonas.pdf"),
       p_ps_tm, width = 6, height = 4)


# ===============================
# Optional combined figure
# ===============================
p_combined <- (p_arabi_len + p_ps_len) / (p_arabi_tm + p_ps_tm) 
#+
#  plot_annotation(tag_levels = 'A')

ggsave(file.path(out_dir, "SuppFig2_withpvalue_InsertSize_FragmentLength_Combined.pdf"),
       p_combined, width = 10, height = 8, dpi = 600)






#for me to see which isolate:
# ===========================================================
# 3. Fragment length distribution — colored by isolate
# ===========================================================
library(ggrepel)

# ===========================================================
# 3. Fragment length distribution — colored and labeled by isolate
# ===========================================================

# --- helper: take rightmost x-point for label placement ---
df_labels <- df_len %>%
  filter(Length > 0) %>%
  group_by(Species, Sample) %>%
  slice_max(Length, n = 1, with_ties = FALSE) %>%
  ungroup()

# ---- Arabidopsis ----
p_arabi_len_labeled <- ggplot(
  df_len %>% filter(Species == "Arabidopsis", Length > 0),
  aes(x = Length, y = Prop, group = Sample, color = Sample)
) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  geom_text_repel(
    data = df_labels %>% filter(Species == "Arabidopsis"),
    aes(label = Sample),
    size = 3,
    hjust = 0,
    nudge_x = 20,
    segment.size = 0.2,
    show.legend = FALSE
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.03)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +  # room for labels
  labs(
    x = "Fragment length (bp)",
    y = "Reads (%)",
    title = expression(italic("Arabidopsis thaliana") ~ " fragment length per isolate"),
    color = "Isolate"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11, face = "bold"),
    legend.position = "none"
  )

# ---- Pseudomonas ----
p_ps_len_labeled <- ggplot(
  df_len %>% filter(Species == "Pseudomonas", Length > 0),
  aes(x = Length, y = Prop, group = Sample, color = Sample)
) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  geom_text_repel(
    data = df_labels %>% filter(Species == "Pseudomonas"),
    aes(label = Sample),
    size = 3,
    hjust = 0,
    nudge_x = 20,
    segment.size = 0.2,
    show.legend = FALSE
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.03)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = "Fragment length (bp)",
    y = "Reads (%)",
    title = expression(italic("Pseudomonas") ~ plain("sp.") ~ " fragment length per isolate"),
    color = "Isolate"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11, face = "bold"),
    legend.position = "none"
  )

# ---- Save ----
ggsave(file.path(out_dir, "SuppFig_FragmentLength_Arabidopsis_ByIsolate_Labeled.pdf"),
       p_arabi_len_labeled, width = 8, height = 5)

ggsave(file.path(out_dir, "SuppFig_FragmentLength_Pseudomonas_ByIsolate_Labeled.pdf"),
       p_ps_len_labeled, width = 8, height = 5)