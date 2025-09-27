# ---- Libraries ----
library(ggplot2)
library(readr)
library(dplyr)
library(data.table)

# ---- Input ----
infile <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step5_additional_analysis_coinfection_Oantigenbreakingisolates_HTFlengthfreq/HTFlengthfreq/HP12_combined_lengths_sorted.txt"
otu5_file <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step5_additional_analysis_coinfection_Oantigenbreakingisolates_HTFlengthfreq/HTFlengthfreq/1355OTU5andp8A2.txt"

setwd('/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step4_additionalanalysis/step3_HTFfreq/')

# ---- Load data ----
df <- read_tsv(infile, show_col_types = FALSE)
otu5_strains <- fread(otu5_file, header = FALSE)$V1  # strain IDs

# ---- Keep only OTU5 strains ----
df <- df %>%
  filter(Sample %in% otu5_strains)

# ---- Define groups & colors ----
length_colors <- c(
  "1830" = "#f6d6ff",
  "1383" = "#638ccc",
  "1803" = "#800233",
  "1245" = "#f9d42a"
)

length_groups <- c("1803","1245","1383","1830")

# ---- Summarise counts ----
df_summary <- df %>%
  filter(!is.na(Length)) %>%   # remove NA here
  count(Length) %>%
  mutate(
    prop = n / sum(n) * 100,
    prop_label = paste0(round(prop, 1), "%"),
    Length = factor(as.character(Length), levels = rev(length_groups))
  )

df_summary <- df_summary %>% filter(!is.na(Length))
# ---- Annotate O-antigen presence/absence ----



# ---- Plot grouped bars ----# ---- Annotate O-antigen presence/absence ----
df_summary <- df_summary %>%
  mutate(
    Oantigen = case_when(
      Length %in% c("1803", "1245") ~ "O-antigen+",
      Length %in% c("1383", "1830") ~ "O-antigen-",
      TRUE ~ "Unknown"
    )
  )
# ---- Plot grouped bars ----
p <- ggplot(df_summary, aes(x = Oantigen, y = prop, fill = Length)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = NA) +
  geom_text(
    aes(label = prop_label, group = Length),
    position = position_dodge(width = 0.8),
    vjust = -0.5, size = 5, fontface = "bold"
  ) +
  scale_fill_manual(values = length_colors) +
  scale_y_continuous(
    breaks = seq(0, 100, by = 10),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(
    x = "",
    y = "Proportion (%)",
    title = "HTF haplotype distribution"
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_blank(),
    legend.position = "top"
  )

# ---- Save ----
ggsave("HTF_length_barplot_OantigenGroups.pdf", plot = p, width = 5, height = 5, dpi = 600)



#for h34 and m53
# ---- Libraries ----
library(dplyr)
library(ggplot2)
library(binom)

# ---- Working directory ----
setwd('/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step4_additionalanalysis/step3_HTFfreq/')

# ================================
# 1. Load modern + historical prop tables
# ================================
modern_df <- read.table(
  "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers/newsummarymodern57_tailocin_kmer_propnorm.tsv",
  header = TRUE, sep = "\t"
)

# remove problematic samples
modern_df <- modern_df %>%
  filter(!Isolate %in% c("p12.F2", "p13.C7", "p6.A10", "p9.C4"))

historical_df <- read.table(
  "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers/newsummaryhistorical40_tailocin_kmer_propnorm.tsv",
  header = TRUE, sep = "\t"
)

# remove isolate 64 (unknown)
historical_df <- historical_df %>%
  filter(Isolate != "64.GBR_1933b_S36")

# add type
modern_df$IsolateType <- "Modern"
historical_df$IsolateType <- "Historical"

df_all <- bind_rows(modern_df, historical_df)

# ================================
# 2. Exclude low-quality isolates
# ================================
exclude <- c("HB0828", "HB0863", "PL0066", "PL0108", "PL0203", "PL0258",
             "PL0027", "PL0053", "30.ESP_1983b", "PL0026")
df_all <- df_all %>% filter(!Isolate %in% exclude)

# ================================
# 3. Assign length groups
# ================================
df_all <- df_all %>%
  mutate(LengthGroup = case_when(
    Reference == "HTF_p7.G11" ~ "1830",
    Reference == "HTF_p25.A12" ~ "1383",
    Reference %in% c("HTF_p5.D5", "HTF_p23.B8",
                     "HTF_p26.D6", "HTF_p25.C2") ~ "1803",
    Reference == "HTF_p21.F9" ~ "1245",
    TRUE ~ "Unknown"
  )) %>%
  filter(LengthGroup != "Unknown")

# ================================
# 4. Get dominant haplotype per isolate
# ================================
dominant_df <- df_all %>%
  group_by(Isolate) %>%
  filter(PropNorm == max(PropNorm, na.rm = TRUE)) %>%
  ungroup()

# ================================
# 5. Count frequencies and CI
# ================================
total_counts <- dominant_df %>%
  distinct(Isolate, IsolateType) %>%
  count(IsolateType, name = "TotalN")

freq_table_length <- dominant_df %>%
  group_by(LengthGroup, IsolateType) %>%
  summarise(Frequency = n(), .groups = "drop") %>%
  left_join(total_counts, by = "IsolateType") %>%
  mutate(
    Proportion = Frequency / TotalN,
    ci = binom.confint(Frequency, TotalN, method = "wilson"),
    Lower = ci$lower, Upper = ci$upper
  ) %>%
  select(LengthGroup, IsolateType, Frequency, TotalN,
         Proportion, Lower, Upper)

# ================================
# 6. Add O-antigen grouping
# ================================
length_colors <- c("1830"="#f6d6ff", "1383"="#638ccc",
                   "1803"="#800233", "1245"="#f9d42a")
length_groups <- rev(c("1803","1245","1383","1830"))

freq_table_length <- freq_table_length %>%
  mutate(
    Oantigen = case_when(
      LengthGroup %in% c("1803","1245") ~ "O-antigen+",
      LengthGroup %in% c("1383","1830") ~ "O-antigen-"
    ),
    LengthGroup = factor(LengthGroup, levels = length_groups)
  )

# ================================
# 7. Function to plot histograms
# ================================
plot_hist <- function(df, type_label) {
  ggplot(df, aes(x = Oantigen, y = Proportion, fill = LengthGroup)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(
      aes(label = paste0(Frequency, "/", TotalN,
                         " (", round(Proportion*100, 1), "%)"),
          group = LengthGroup),
      position = position_dodge(width = 0.8),
      vjust = -0.5, size = 5, fontface = "bold"
    ) +
    scale_fill_manual(values = length_colors) +
    scale_y_continuous(
      breaks = seq(0, 1, by = 0.1),
      labels = function(x) paste0(x*100, "%"),
      expand = expansion(mult = c(0, 0.12))
    ) +
    labs(
      title = paste("HTF haplotype distribution (", type_label, ")", sep = ""),
      x = "", y = "Proportion of isolates"
    ) +
    theme_classic(base_size = 18) +
    theme(
      axis.text = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 16, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      legend.title = element_blank(),
      legend.position = "top"
    )
}

# ================================
# 8. Generate and save plots
# ================================
for (type in c("Historical","Modern")) {
  df_sub <- freq_table_length %>% filter(IsolateType == type)
  p <- plot_hist(df_sub, type)
  ggsave(paste0("HTF_length_histogram_", type, "_Oantigen.pdf"),
         plot = p, width = 7, height = 5, dpi = 600)
}





#combine h34 and m1350:
# ---- Libraries ----
library(dplyr)
library(ggplot2)
library(patchwork)   # for combining plots

setwd('/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step4_additionalanalysis/step3_HTFfreq/')

# ---- Colors and order ----
length_colors <- c("1830"="#f6d6ff", "1383"="#638ccc",
                   "1803"="#800233", "1245"="#f9d42a")
length_groups <- rev(c("1803","1245","1383","1830"))

# ---- Prepare Modern (df_summary, 1350 OTU5 isolates) ----
df_mod <- df_summary %>%
  mutate(
    Oantigen = case_when(
      Length %in% c("1803","1245") ~ "O-antigen+",
      Length %in% c("1383","1830") ~ "O-antigen-"
    ),
    Length = factor(Length, levels = length_groups),
    IsolateType = "Modern"
  ) %>%
  rename(LengthGroup = Length, Frequency = n, Proportion = prop)

# ---- Prepare Historical (freq_table_length, 34 isolates) ----
freq_table_length_hist <- freq_table_length %>%
  filter(TotalN == 34)
df_hist <- freq_table_length_hist %>%
  mutate(
    Oantigen = case_when(
      LengthGroup %in% c("1803","1245") ~ "O-antigen+",
      LengthGroup %in% c("1383","1830") ~ "O-antigen-"
    ),
    LengthGroup = factor(LengthGroup, levels = length_groups),
    IsolateType = "Historical"
  ) %>%
  mutate(Proportion = Proportion*100)  # convert to percent for consistency

# ---- Function to plot one panel ----
plot_panel <- function(df, title_label, alpha_val = 1) {
  ggplot(df, aes(x = Oantigen, y = Proportion, fill = LengthGroup)) +
    geom_col(position = position_dodge(width = 0.8),
             width = 0.7, alpha = alpha_val) +
    geom_text(
      aes(label = paste0(round(Proportion, 1), "%"), group = LengthGroup),
      position = position_dodge(width = 0.8),
      vjust = -0.5, size = 5, fontface = "bold"
    ) +
    scale_fill_manual(values = length_colors) +
    scale_y_continuous(
      limits = c(0, 45),
      breaks = seq(0, 45, by = 10),
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(title = title_label, x = "", y = "Proportion (%)") +
    theme_classic(base_size = 18) +
    theme(
      axis.text = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 16, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      legend.title = element_blank()
    )
}
library(patchwork)

# ---- Fix proportions ----
# Historical values are fractions (0–1) -> convert once
df_hist <- df_hist %>% mutate(Proportion = Proportion )

# Modern already in percentage from df_summary
df_mod  <- df_mod %>% mutate(Proportion = Proportion)
# ---- Generate panels ----
# ---- Generate panels ----
# ---- Generate panels ----
p_hist <- plot_panel(df_hist, "Historical", alpha_val = 1) +
  guides(fill = guide_legend(title = NULL))

p_mod  <- plot_panel(df_mod, "Modern", alpha_val = 0.7) +
  guides(fill = "none") +
  theme(
    axis.title.y = element_blank(),   # remove y title
    axis.text.y  = element_blank(),   # remove y ticks
    axis.ticks.y = element_blank(),   # remove y tick marks
    axis.line.y  = element_blank()    # remove y axis line
  )

# ---- Combine with shared legend, centered ----
p_combined <- (p_hist + p_mod) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.box.just = "center",
    legend.box = "horizontal"
  )

# ---- Save ----
ggsave("HTF_length_histogram_Historical_vs_Modern_panels_y45.pdf",
       plot = p_combined, width = 12, height = 6, dpi = 600)





#box:
# ---- Libraries ----
library(dplyr)
library(ggplot2)
library(binom)

# ============================
# 1. Historical data (from freq_table_length)
# ============================
freq_table_length_hist <- freq_table_length %>% filter(TotalN == 34)

df_hist <- freq_table_length_hist %>%
  mutate(
    Oantigen = case_when(
      LengthGroup %in% c("1803", "1245") ~ "O-antigen+",
      LengthGroup %in% c("1383", "1830") ~ "O-antigen−"
    )
  ) %>%
  group_by(IsolateType, Oantigen, TotalN) %>%
  summarise(Frequency = sum(Frequency), .groups = "drop") %>%
  mutate(
    Proportion = Frequency / TotalN,
    ci = binom.confint(Frequency, TotalN, method = "wilson"),
    Lower = ci$lower, Upper = ci$upper
  ) %>%
  select(IsolateType, Oantigen, Frequency, TotalN, Proportion, Lower, Upper)

# ============================
# 2. Modern data (from df_summary, 1350 OTU5)
# ============================
df_mod <- df_summary %>%
  mutate(
    Oantigen = case_when(
      Length %in% c("1803", "1245") ~ "O-antigen+",
      Length %in% c("1383", "1830") ~ "O-antigen−"
    )
  ) %>%
  group_by(Oantigen) %>%
  summarise(Frequency = sum(n), TotalN = 1350, .groups = "drop") %>%
  mutate(
    Proportion = Frequency / TotalN,
    ci = binom.confint(Frequency, TotalN, method = "wilson"),
    Lower = ci$lower, Upper = ci$upper,
    IsolateType = "Modern"
  ) %>%
  select(IsolateType, Oantigen, Frequency, TotalN, Proportion, Lower, Upper)

# ============================
# 3. Combine historical + modern
# ============================
df_combined <- bind_rows(df_hist, df_mod) %>%
  mutate(
    IsolateType = factor(IsolateType, levels = c("Historical", "Modern")),
    Oantigen = factor(Oantigen, levels = c("O-antigen−", "O-antigen+")),
    Label = paste0(Frequency, "/", TotalN, " (", round(Proportion*100, 1), "%)")
  )

# ============================
# 4. Plot
# ============================
# ---- Plot with percentages only ----
# Increase dodge width to spread O-antigen− and O-antigen+ apart
dodge <- position_dodge(width = 0.8)

p <- ggplot(df_combined, aes(x = IsolateType, y = Proportion,
                             color = Oantigen)) +
  geom_point(position = dodge,
             size = 3,
             aes(alpha = IsolateType),
             show.legend = TRUE) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                width = 0.1,
                position = dodge,
                show.legend = FALSE) +
  geom_text(aes(label = paste0(round(Proportion*100, 1), "%")),
            position = dodge,
            vjust = -1.2, size = 4, fontface = "bold",
            show.legend = FALSE) +
  scale_color_manual(values = c("O-antigen−" = "#9e9ac8",
                                "O-antigen+" = "#c94f7c")) +
  scale_alpha_manual(values = c("Historical" = 1, "Modern" = 0.7),
                     guide = "none") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(x = "", y = "Proportion (%)") +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 11, face = "bold")
  )

ggsave("Oantigen_freq_hist_vs_modern1350_separated.pdf",
       plot = p, width = 9, height = 3.5, dpi = 600)