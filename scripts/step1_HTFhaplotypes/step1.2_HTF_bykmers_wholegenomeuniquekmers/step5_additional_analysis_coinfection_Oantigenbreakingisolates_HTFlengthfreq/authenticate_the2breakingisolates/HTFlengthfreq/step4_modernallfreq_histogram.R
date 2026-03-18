# ---- Libraries ----
library(ggplot2)
library(readr)
library(dplyr)
library(data.table)
#in total 1350 modern samples in whcih the four length sum up to 1312

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



#for h36 and m53
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
  "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers/newsummaryhistorical49_tailocin_kmer_propnorm_05sd.tsv",
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
exclude <- c("HB0828", "HB0863", "PL0066","PL0108", "PL0203", "PL0258", "64.GBR_1933b_S36", "PL0065","PL0026","PL0027","PL0053")
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





#combine h38 and m1350:
# ---- Libraries ----
library(dplyr)
library(ggplot2)
library(patchwork)   # for combining plots

setwd('/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step4_additionalanalysis/step3_HTFfreq/')

# ---- Colors and order ----
length_colors <- c("1830"="#f6d6ff", "1383"="#638ccc",
                   "1803"="#800233", "1245"="#f9d42a")
length_groups <- rev(c("1803","1245","1383","1830"))
# Clean workspace
#rm(list = ls())

# Reload df_summary creation code first (from your infile)

# Confirm
colnames(df_summary)
df_summary$Length
# Then run this block
df_mod <- df_summary %>%
  mutate(
    Oantigen = case_when(
      .data$Length %in% c("1803", "1245") ~ "O-antigen+",
      .data$Length %in% c("1383", "1830") ~ "O-antigen-"
    ),
    Length = factor(.data$Length, levels = length_groups),
    IsolateType = "Modern"
  ) %>%
  rename(LengthGroup = Length,
         Frequency = n,
         Proportion = prop)

colnames(df_summary)

# ---- Prepare Historical (freq_table_length, 38 isolates) ----
freq_table_length_hist <- freq_table_length %>%
  filter(TotalN == 38)
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
      limits = c(0, 50),
      breaks = seq(0, 50, by = 10),
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

p_mod  <- plot_panel(df_mod, "Modern", alpha_val = 1) +
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




#box four groups:
# ---- Libraries ----
library(dplyr)
library(ggplot2)
library(binom)

length_colors <- c("1830"="#f6d6ff", "1383"="#638ccc",
                   "1803"="#800233", "1245"="#f9d42a")
length_groups <- rev(c("1803","1245","1383","1830"))
# =========================================================
# 1. Historical (36 isolates)
# =========================================================
df_hist4 <- freq_table_length %>%
  filter(TotalN == 38) %>%
  mutate(
    IsolateType = "Historical",
    LengthGroup = factor(LengthGroup, levels = length_groups),
    Proportion = Proportion * 100,
    Lower = Lower * 100,
    Upper = Upper * 100,
    Oantigen = case_when(
      LengthGroup %in% c("1803","1245") ~ "O-antigen+",
      LengthGroup %in% c("1383","1830") ~ "O-antigen-"
    )
  )

# =========================================================
# 2. Modern (1350 but 1312 four legth OTU5 isolates)
# =========================================================
df_mod4 <- df_summary %>%
  mutate(
    LengthGroup = factor(as.character(Length), levels = length_groups),
    IsolateType = "Modern",
    TotalN = 1312,
    Oantigen = case_when(
      LengthGroup %in% c("1803","1245") ~ "O-antigen+",
      LengthGroup %in% c("1383","1830") ~ "O-antigen-"
    )
  ) %>%
  group_by(Oantigen, LengthGroup, IsolateType) %>%
  summarise(Frequency = sum(n), TotalN = unique(TotalN), .groups = "drop") %>%
  mutate(
    ci = binom.confint(Frequency, TotalN, method = "wilson"),
    Proportion = ci$mean * 100,
    Lower = ci$lower * 100,
    Upper = ci$upper * 100
  )

# =========================================================
# 3. Combine
# =========================================================
df_combined4 <- bind_rows(df_hist4, df_mod4) %>%
  mutate(IsolateType = factor(IsolateType, levels = c("Historical","Modern")))

# =========================================================
# 4. Function to make one panel (dot + CI)
# =========================================================
plot_panel <- function(df, panel_title) {
  dodge <- position_dodge(width = 0.8)
  
  ggplot(df, aes(x = Oantigen, y = Proportion, color = LengthGroup)) +
    geom_point(position = position_dodge(width = 0.8),
               width = 0.9, size = 4) +   # slightly larger dots
    geom_errorbar(aes(ymin = Lower, ymax = Upper, group = LengthGroup),
                  width = 0.2, position = dodge, size = 0.7, show.legend = FALSE) + # thicker CI lines
    geom_text(aes(label = paste0(round(Proportion,1), "%"),
                  group = LengthGroup),
              position = dodge, vjust = -2, size = 5, fontface = "bold",color = "black") +
    scale_color_manual(values = length_colors,
                       guide = guide_legend(override.aes = list(size = 5))) + # bigger legend dots
    scale_y_continuous(labels = function(x) paste0(x,"%"), limits = c(0,50)) +
    labs(title = panel_title, x = "", y = "Proportion (%)") +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(size=15,hjust=0.5, face="bold"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.key = element_blank(),
      axis.text = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 14, face = "bold")
      
    )
}

# =========================================================
# 5. Build two panels
# =========================================================
p_hist <- plot_panel(filter(df_combined4, IsolateType=="Historical"), "Historical")
p_mod  <- plot_panel(filter(df_combined4, IsolateType=="Modern"), "Modern") +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank()
  )

# =========================================================
# 6. Combine
# =========================================================
library(patchwork)
p_final <- (p_hist + p_mod) + plot_layout(guides="collect")

# ---- Save (shorter height) ----
ggsave("HTF_length_hist_vs_modern_panels_dotCI_short.pdf",
       plot = p_final, width = 14, height = 4, dpi = 600)









#box two groups:
# ---- Libraries ----
library(dplyr)
library(ggplot2)
library(binom)

# ============================
# 1. Historical data (from freq_table_length)
# ============================
freq_table_length_hist <- freq_table_length %>% filter(TotalN == 38)

df_hist <- freq_table_length_hist %>%
  mutate(
    Oantigen = case_when(
      LengthGroup %in% c("1803", "1245") ~ "O-antigen+",
      LengthGroup %in% c("1383", "1830") ~ "O-antigen-"
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
      Length %in% c("1383", "1830") ~ "O-antigen-"
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
    Oantigen = factor(Oantigen, levels = c("O-antigen-", "O-antigen+")),
    Label = paste0(Frequency, "/", TotalN, " (", round(Proportion*100, 1), "%)")
  )

# ============================
# 4. Plot
# ============================
# ---- Plot with percentages only ----
# Increase dodge width to spread O-antigen- and O-antigen+ apart
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
  geom_text(aes(label = paste0(round(Proportion*100, 1), "%"),
                group = Oantigen),              # <-- crucial for correct position
            position = dodge,
            vjust = -1.0, size = 4, fontface = "bold",
            color = "black", show.legend = FALSE) +
  scale_color_manual(values = c("O-antigen-" = "#9e9ac8",
                                "O-antigen+" = "#c94f7c")) +
  scale_alpha_manual(values = c("Historical" = 1, "Modern" = 1),
                     guide = "none") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(x = "", y = "Proportion (%)") +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 11, face = "bold")
  )

ggsave("Oantigen_freq_hist_vs_modern1350_separated.pdf",
       plot = p, width = 9, height = 3.5, dpi = 600)




#fisher exact test
# build contingency table
mat <- matrix(c(
  df_combined$Frequency[df_combined$IsolateType == "Historical" & df_combined$Oantigen == "O-antigen+"],
  df_combined$Frequency[df_combined$IsolateType == "Historical" & df_combined$Oantigen == "O-antigen-"],
  df_combined$Frequency[df_combined$IsolateType == "Modern" & df_combined$Oantigen == "O-antigen+"],
  df_combined$Frequency[df_combined$IsolateType == "Modern" & df_combined$Oantigen == "O-antigen-"]
), 
nrow = 2, byrow = TRUE)

# add row/col names
rownames(mat) <- c("Historical", "Modern")
colnames(mat) <- c("O-antigen+", "O-antigen-")

# view contingency table
print(mat)

# Fisher's exact test
fisher.test(mat)












# =========================================================
# --- Bootstrap resampling (38 isolates per replicate) ---
# =========================================================
library(dplyr)
library(ggplot2)
library(binom)
library(patchwork)

set.seed(1)

# ---- Define colors ----
length_colors <- c("1830"="#f6d6ff", "1383"="#638ccc",
                   "1803"="#800233", "1245"="#f9d42a")
length_groups <- rev(c("1803","1245","1383","1830"))

# ---- Historical dataset (from previous freq_table_length) ----
df_hist4 <- freq_table_length %>%
  filter(TotalN == 38) %>%
  mutate(
    IsolateType = "Historical",
    LengthGroup = factor(LengthGroup, levels = length_groups),
    Proportion = Proportion * 100,
    Lower = Lower * 100,
    Upper = Upper * 100,
    Oantigen = case_when(
      LengthGroup %in% c("1803","1245") ~ "O-antigen+",
      LengthGroup %in% c("1383","1830") ~ "O-antigen-"
    )
  )

# =========================================================
# --- 1. Bootstrap Modern 1350 → 38 resamples ---
# =========================================================
# df_summary contains per-length n counts
# Expand to simulated isolate-level data
modern_isolates <- df_summary %>%
  mutate(
    LengthGroup = as.character(Length),
    Oantigen = case_when(
      LengthGroup %in% c("1803","1245") ~ "O-antigen+",
      LengthGroup %in% c("1383","1830") ~ "O-antigen-"
    )
  ) %>%
  select(LengthGroup, Oantigen, n)

# Expand to isolate-level (each row = 1 isolate)
modern_expanded <- modern_isolates %>%
  group_by(LengthGroup, Oantigen) %>%
  summarise(
    Isolates = list(rep(LengthGroup, n)),
    .groups = "drop"
  ) %>%
  tidyr::unnest(Isolates)

# ---------------------------------------------------------
# Bootstrap resampling of 38 isolates from 1350
# ---------------------------------------------------------
n_boot <- 1000
sample_size <- 38

boot_results_htf <- replicate(n_boot, {
  sampled <- sample(modern_expanded$LengthGroup, sample_size, replace = TRUE)
  tbl <- table(factor(sampled, levels = length_groups))
  as.numeric(tbl) / sample_size * 100
})

# Convert to tidy dataframe: mean and 95% CI
df_mod4_boot <- as.data.frame(t(boot_results_htf))
colnames(df_mod4_boot) <- length_groups

df_mod4_mean <- df_mod4_boot %>%
  summarise(across(everything(), list(
    mean = mean,
    lower = ~quantile(.x, 0.025),
    upper = ~quantile(.x, 0.975)
  ))) %>%
  tidyr::pivot_longer(cols = everything(),
                      names_to = c("LengthGroup",".value"),
                      names_sep = "_") %>%
  mutate(IsolateType = "Modern")

# Add O-antigen group info
df_mod4_mean <- df_mod4_mean %>%
  mutate(
    Oantigen = case_when(
      LengthGroup %in% c("1803","1245") ~ "O-antigen+",
      LengthGroup %in% c("1383","1830") ~ "O-antigen-"
    ),
    LengthGroup = factor(LengthGroup, levels = length_groups)
  )

# =========================================================
# --- 2. Combine Historical and Modern (bootstrap mean) ---
# =========================================================
df_combined4 <- bind_rows(
  df_hist4 %>%
    select(LengthGroup, Proportion, Lower, Upper, IsolateType, Oantigen),
  df_mod4_mean %>%
    rename(Proportion = mean, Lower = lower, Upper = upper)
) %>%
  mutate(IsolateType = factor(IsolateType, levels = c("Historical","Modern")))

# =========================================================
# --- 3. Plot: four HTF groups ---
# =========================================================
plot_panel_htf <- function(df, panel_title) {
  dodge <- position_dodge(width = 0.8)
  
  ggplot(df, aes(x = Oantigen, y = Proportion, color = LengthGroup)) +
    geom_point(position = dodge, size = 4) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper, group = LengthGroup),
                  width = 0.2, position = dodge, size = 0.7, show.legend = FALSE) +
    geom_text(aes(label = paste0(round(Proportion,1), "%"),
                  group = LengthGroup),
              position = dodge, vjust = -2, size = 5, fontface = "bold", color = "black") +
    scale_color_manual(values = length_colors,
                       guide = guide_legend(override.aes = list(size = 5))) +
    scale_y_continuous(labels = function(x) paste0(x,"%"), limits = c(0,60)) +
    labs(title = panel_title, x = "", y = "Proportion (%)") +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(size=15,hjust=0.5, face="bold"),
      legend.position = "top",
      legend.title = element_blank(),
      axis.text = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 14, face = "bold")
    )
}

p_hist_htf <- plot_panel_htf(filter(df_combined4, IsolateType=="Historical"), "Historical")
p_mod_htf  <- plot_panel_htf(filter(df_combined4, IsolateType=="Modern"), "Modern") +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y  = element_blank())

p_final_htf <- (p_hist_htf + p_mod_htf) + plot_layout(guides="collect")

ggsave("fig1b_a_resamplemodern38_groupbyHTF.pdf",
       plot = p_final_htf, width = 14, height = 4, dpi = 600)


# =========================================================
# --- 4. Collapse into two O-antigen groups ---
# =========================================================
# Historical


df_hist2 <- df_hist4 %>%
  mutate(
    Oantigen = case_when(
      LengthGroup %in% c("1803", "1245") ~ "O-antigen+",
      LengthGroup %in% c("1383", "1830") ~ "O-antigen-"
    )
  ) %>%
  group_by(IsolateType, Oantigen, TotalN) %>%
  summarise(Frequency = sum(Frequency), .groups = "drop") %>%
  mutate(
    Proportion = Frequency / TotalN,
    ci = binom.confint(Frequency, TotalN, method = "wilson"),
    Lower = ci$lower, Upper = ci$upper,
    Proportion = Proportion * 100,
    Lower = Lower * 100,
    Upper = Upper * 100,
  ) %>%
  select(IsolateType, Oantigen, Frequency, TotalN, Proportion, Lower, Upper)

# Modern (bootstrap)
df_mod2_boot <- df_mod4_boot %>%
  apply(2, function(x) NA) # just placeholder to satisfy structure
# recompute from df_mod4_boot
boot_results_oantigen <- apply(boot_results_htf, 2, function(x) NULL)

boot_results_oantigen <- replicate(n_boot, {
  sampled <- sample(modern_expanded$Oantigen, sample_size, replace = TRUE)
  tbl <- table(factor(sampled, levels = c("O-antigen+","O-antigen-")))
  as.numeric(tbl) / sample_size * 100
})

df_mod2_boot <- as.data.frame(t(boot_results_oantigen))
colnames(df_mod2_boot) <- c("O-antigen+","O-antigen-")

df_mod2_mean <- df_mod2_boot %>%
  summarise(across(everything(), list(
    mean = mean,
    lower = ~quantile(.x, 0.025),
    upper = ~quantile(.x, 0.975)
  ))) %>%
  tidyr::pivot_longer(cols = everything(),
                      names_to = c("Oantigen",".value"),
                      names_sep = "_") %>%
  mutate(IsolateType = "Modern")

df_combined2 <- bind_rows(df_hist2,
                          df_mod2_mean %>%
                            rename(Proportion = mean,
                                   Lower = lower,
                                   Upper = upper))

# =========================================================
# --- 5. Plot: two O-antigen groups ---
# =========================================================
dodge <- position_dodge(width = 0.6)

p2 <- ggplot(df_combined2, aes(x = IsolateType, y = Proportion, color = Oantigen)) +
  geom_point(position = dodge, size = 3) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                width = 0.15, position = dodge, size = 0.7) +
  geom_text(aes(label = paste0(round(Proportion,1), "%"), group = Oantigen),
            position = dodge, vjust = -1.2, size = 4, fontface = "bold", color = "black") +
  scale_color_manual(values = c("O-antigen-"="#9e9ac8","O-antigen+"="#c94f7c")) +
  scale_y_continuous(labels = function(x) paste0(x,"%"), limits = c(0,100)) +
  labs(x = "", y = "Proportion (%)", title = "O-antigen groups (resampled Modern n=38)") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size=15,hjust=0.5, face="bold"),
    legend.position = "top",
    legend.title = element_blank(),
    axis.text = element_text(size = 13, face = "bold"),
    axis.title = element_text(size = 13, face = "bold")
  )

ggsave("fig1b_b_resamplemodern38_groupbyOantigen.pdf",
       plot = p2, width = 9, height = 3.5, dpi = 600)









# =========================================================
# --- Country-level bootstrap (top 4 historical countries) ---
# =========================================================
# File with sample names, dates, countries etc.
samples_path <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step5_additional_analysis_coinfection_Oantigenbreakingisolates_HTFlengthfreq/HTFlengthfreq/h49PL0087NAtmpm83_withdatesandlocs_uniq copy.txt"

# File with per-isolate dominant HTF/O-antigen assignments
dominant_path <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step4_additionalanalysis/step3_HTFfreq/step2_3_HTF_oantigen_dominant_table.tsv"


# =========================================================
# --- Bootstrap modern-downsampled (3 countries, no Sweden) ---
# --- Fisher test written to table only ---
# =========================================================
library(dplyr)
library(ggplot2)
library(binom)
library(tidyr)
library(scales)

set.seed(42)
outdir_country <- "top3_bycountry_moderndownsampled_fisher"
if (!dir.exists(outdir_country)) dir.create(outdir_country)

# ---- Load datasets ----
samples <- read.table(samples_path, header = TRUE, sep = "\t")
dominant_df <- read.table(dominant_path, header = TRUE, sep = "\t")

samples$country[samples$country == "United Kingdom"] <- "UK"
samples$country[samples$country == "Congo"] <- "Democratic Republic of the Congo"

dominant_df <- dominant_df %>%
  mutate(LengthGroup = case_when(
    Reference == "HTF_p7.G11" ~ "1830",
    Reference == "HTF_p25.A12" ~ "1383",
    Reference %in% c("HTF_p5.D5","HTF_p23.B8",
                     "HTF_p26.D6","HTF_p25.C2") ~ "1803",
    Reference == "HTF_p21.F9" ~ "1245",
    TRUE ~ NA_character_
  ),
  Oantigen = case_when(
    LengthGroup %in% c("1803","1245") ~ "O-antigen+",
    LengthGroup %in% c("1830","1383") ~ "O-antigen-",
    TRUE ~ NA_character_
  ))

df_hist <- left_join(samples,
                     dominant_df[, c("Isolate","LengthGroup","IsolateType","Oantigen")],
                     by = c("samplename"="Isolate")) %>%
  filter(IsolateType == "Historical", !is.na(LengthGroup))

cat("✅ Historical isolates with HTF:", nrow(df_hist), "\n")

# ---- Identify top countries (exclude Sweden) ----
top_countries <- df_hist %>%
  count(country, name = "n") %>%
  arrange(desc(n)) %>%
  filter(!grepl("Sweden", country, ignore.case = TRUE)) %>%
  slice_head(n = 3)
print(top_countries)

# ---- Expand modern data for bootstrap ----
modern_expanded <- df_summary %>%
  mutate(
    LengthGroup = as.character(Length),
    Oantigen = case_when(
      LengthGroup %in% c("1803", "1245") ~ "O-antigen+",
      LengthGroup %in% c("1383", "1830") ~ "O-antigen-"
    )
  ) %>%
  select(Oantigen, n) %>%
  group_by(Oantigen) %>%
  summarise(Isolates = list(rep(Oantigen, sum(n))), .groups = "drop") %>%
  tidyr::unnest(Isolates)


plot_ci_cache <- list()   # will store per-country CI rows used in plots
# ---- Function to run bootstrap + fisher test + plot ----
plot_country_boot <- function(cname, n_hist) {
  # ---- Historical subset ----
  df_hist_country <- df_hist %>%
    filter(country == cname) %>%
    count(Oantigen, name = "Frequency") %>%
    mutate(
      TotalN = n_hist,
      Proportion = Frequency / TotalN,
      ci = binom.confint(Frequency, TotalN, method = "wilson"),
      Lower = ci$lower * 100,
      Upper = ci$upper * 100,
      Proportion = Proportion * 100,
      IsolateType = "Historical"
    ) %>%
    select(IsolateType, Oantigen, Frequency, TotalN, Proportion, Lower, Upper)
  
  # ---- Bootstrap Modern downsampled to n_hist ----
  n_boot <- 1000
  boot_results <- replicate(n_boot, {
    sampled <- sample(modern_expanded$Isolates, n_hist, replace = TRUE)
    tbl <- table(factor(sampled, levels = c("O-antigen+", "O-antigen-")))
    as.numeric(tbl) / n_hist * 100
  })
  
  df_mod_boot <- as.data.frame(t(boot_results))
  colnames(df_mod_boot) <- c("O-antigen+", "O-antigen-")
  
  df_mod_ci <- df_mod_boot %>%
    summarise(across(everything(), list(
      mean  = mean,
      lower = ~quantile(.x, 0.025),
      upper = ~quantile(.x, 0.975)
    ))) %>%
    tidyr::pivot_longer(cols = everything(),
                        names_to = c("Oantigen", ".value"),
                        names_sep = "_") %>%
    mutate(IsolateType = "Modern") %>%
    transmute(
      IsolateType,
      Oantigen,
      # expected counts from the same % that drive the plot:
      Frequency  = (mean  / 100) * n_hist,
      TotalN     = n_hist,
      Proportion = mean,              # % used in the plot (point)
      Lower      = lower,             # % CI low used in the plot
      Upper      = upper,             # % CI high used in the plot
      # (optional) count-scale CIs, if you want them in the table too:
      CountLower = (lower / 100) * n_hist,
      CountUpper = (upper / 100) * n_hist
    )
  # ---- Combine for plotting ----
  df_combined <- bind_rows(
    df_hist_country,
    df_mod_ci %>% rename(Proportion = Proportion, Lower = Lower, Upper = Upper)
  ) %>%
    mutate(IsolateType = factor(IsolateType, levels = c("Historical", "Modern")))
  
  # ---- CACHE the exact CI rows used in the plot ----
    plot_ci_cache[[cname]] <<- bind_rows(
    mutate(df_hist_country, Country = cname),
    mutate(df_mod_ci,      Country = cname)
  ) %>%
    select(Country, IsolateType, Oantigen, TotalN,
           Frequency, Proportion, Lower, Upper,
           dplyr::any_of(c("CountLower","CountUpper")))
  # ---- Plot (unchanged) ----
  dodge <- position_dodge(width = 0.6)
  p <- ggplot(df_combined, aes(x = IsolateType, y = Proportion, color = Oantigen)) +
    geom_point(position = dodge, size = 3) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper),
                  width = 0.15, position = dodge, size = 0.7) +
    geom_text(aes(label = paste0(round(Proportion, 1), "%"), group = Oantigen),
              position = dodge, vjust = -1.2, size = 4, fontface = "bold", color = "black") +
    scale_color_manual(values = c("O-antigen-" = "#9e9ac8", "O-antigen+" = "#c94f7c")) +
    scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100)) +
    labs(x = "", y = "Proportion (%)",
         title = paste0(cname, " (n=", n_hist, " historical)")) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
      legend.position = "top",
      legend.title = element_blank(),
      axis.text = element_text(size = 13, face = "bold"),
      axis.title = element_text(size = 13, face = "bold")
    )
  
  outfile <- file.path(outdir_country,
                       paste0("Oantigen_freq_", cname,
                              "_hist_vs_resampleModern", n_hist, ".pdf"))
  ggsave(outfile, plot = p, width = 9, height = 3.5, dpi = 600)
  message("✅ Saved: ", outfile)
}
# ---- Run over 3 countries ----
for (i in seq_len(nrow(top_countries))) {
  cname <- top_countries$country[i]
  n_hist <- top_countries$n[i]
  plot_country_boot(cname, n_hist)
}

# =========================================================
# --- Write the exact CI values used in the plots ---
# =========================================================

ci_table <- dplyr::bind_rows(plot_ci_cache) %>%
  dplyr::arrange(Country,
                 factor(IsolateType, levels = c("Historical","Modern")),
                 factor(Oantigen, levels = c("O-antigen+","O-antigen-")))

outfile_ci <- file.path(outdir_country, "Oantigen_plot_CI_top3_exact.txt")
write.table(ci_table, file = outfile_ci, sep = "\t", quote = FALSE, row.names = FALSE)
cat("✅ CI table saved to:\n", outfile_ci, "\n")




# =========================================================
# --- Within-dataset only: Historical (binomial) & Modern (prop.test)
# =========================================================
library(dplyr)
library(tidyr)

set.seed(42)

within_results <- data.frame(
  Country = character(),
  Dataset = character(),
  Plus = numeric(),
  Minus = numeric(),
  Test = character(),
  P_Value = numeric(),
  Significance = character(),
  stringsAsFactors = FALSE
)

sig_code <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  "ns"
}

analyze_country_within <- function(cname, n_hist) {
  # --- Historical counts (ensure both levels exist) ---
  df_hist_country <- df_hist %>%
    filter(country == cname) %>%
    count(Oantigen, name = "Frequency") %>%
    complete(Oantigen = c("O-antigen+","O-antigen-"), fill = list(Frequency = 0))
  
  hist_plus  <- df_hist_country$Frequency[df_hist_country$Oantigen == "O-antigen+" ]
  hist_minus <- df_hist_country$Frequency[df_hist_country$Oantigen == "O-antigen-" ]
  n_hist_tot <- hist_plus + hist_minus
  
  # exact binomial (H)
  bt_hist <- binom.test(hist_plus, n_hist_tot, p = 0.5)
  within_results <<- rbind(within_results, data.frame(
    Country = cname,
    Dataset = "Historical",
    Plus = hist_plus,
    Minus = hist_minus,
    Test = "binomial",
    P_Value = bt_hist$p.value,
    Significance = sig_code(bt_hist$p.value)
  ))
  
  # --- Modern expected counts via bootstrap (keep fractional) ---
  n_boot <- 1000
  boot_results <- replicate(n_boot, {
    sampled <- sample(modern_expanded$Isolates, n_hist, replace = TRUE)
    tbl <- table(factor(sampled, levels = c("O-antigen+","O-antigen-")))
    as.numeric(tbl)
  })
  # means over replicates (fractional)
  df_mod_mean <- as.data.frame(t(boot_results)) %>%
    summarise(plus = mean(V1), minus = mean(V2))
  mod_plus  <- df_mod_mean$plus
  mod_minus <- df_mod_mean$minus
  n_mod_tot <- mod_plus + mod_minus
  
  # prop.test vs 0.5 (Modern, fractional)
  bt_mod <- prop.test(mod_plus, n_mod_tot, p = 0.5, correct = FALSE)
  within_results <<- rbind(within_results, data.frame(
    Country = cname,
    Dataset = "Modern",
    Plus = mod_plus,
    Minus = mod_minus,
    Test = "prop.test",
    P_Value = bt_mod$p.value,
    Significance = sig_code(bt_mod$p.value)
  ))
  
  message("✅ ", cname,
          " | Hist binom p=", signif(bt_hist$p.value,3),
          " | Mod prop.test p=", signif(bt_mod$p.value,3))
}

# ---- Run across the same three countries (top_countries already defined) ----
for (i in seq_len(nrow(top_countries))) {
  cname <- top_countries$country[i]
  n_hist <- top_countries$n[i]
  analyze_country_within(cname, n_hist)
}

# ---- Save within-only summary ----
within_results <- within_results %>%
  mutate(P_Value = signif(P_Value, 6)) %>%
  arrange(Country, Dataset)

outfile_within <- file.path(outdir_country, "Binomial_PropTest_within_summary.txt")
write.table(within_results, file = outfile_within,
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n✅ Within-only results saved to:\n", outfile_within, "\n")
# =========================================================
# --- Between Historical vs Modern (prop.test for floating modern) ---
# =========================================================

between_results <- data.frame(
  Country = character(),
  Hist_Plus = numeric(),
  Hist_Minus = numeric(),
  Mod_Plus = numeric(),
  Mod_Minus = numeric(),
  P_Value = numeric(),
  Significance = character(),
  stringsAsFactors = FALSE
)

sig_code <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  "ns"
}

analyze_between_hist_mod <- function(cname, n_hist) {
  # --- Historical integer counts ---
  df_hist_country <- df_hist %>%
    filter(country == cname) %>%
    count(Oantigen, name = "Frequency") %>%
    tidyr::complete(Oantigen = c("O-antigen+","O-antigen-"), fill = list(Frequency = 0))
  hist_plus  <- df_hist_country$Frequency[df_hist_country$Oantigen == "O-antigen+"]
  hist_minus <- df_hist_country$Frequency[df_hist_country$Oantigen == "O-antigen-"]
  n_hist_tot <- hist_plus + hist_minus
  
  # --- Modern expected counts via bootstrap mean (fractional) ---
  n_boot <- 1000
  boot_results <- replicate(n_boot, {
    sampled <- sample(modern_expanded$Isolates, n_hist, replace = TRUE)
    tbl <- table(factor(sampled, levels = c("O-antigen+","O-antigen-")))
    as.numeric(tbl)
  })
  df_mod_mean <- as.data.frame(t(boot_results)) %>%
    summarise(plus = mean(V1), minus = mean(V2))
  mod_plus  <- df_mod_mean$plus
  mod_minus <- df_mod_mean$minus
  n_mod_tot <- mod_plus + mod_minus
  
  # --- Proportion test between historical and modern ---
  counts <- c(hist_plus, mod_plus)
  totals <- c(n_hist_tot, n_mod_tot)
  bt_between <- prop.test(counts, totals, correct = FALSE)
  
  between_results <<- rbind(between_results, data.frame(
    Country = cname,
    Hist_Plus = hist_plus,
    Hist_Minus = hist_minus,
    Mod_Plus = mod_plus,
    Mod_Minus = mod_minus,
    P_Value = bt_between$p.value,
    Significance = sig_code(bt_between$p.value)
  ))
  
  message("✅ ", cname, 
          " | Between Hist vs Mod prop.test p = ", signif(bt_between$p.value, 3))
}

# ---- Run across same top_countries ----
for (i in seq_len(nrow(top_countries))) {
  cname <- top_countries$country[i]
  n_hist <- top_countries$n[i]
  analyze_between_hist_mod(cname, n_hist)
}

# ---- Save summary ----
between_results <- between_results %>%
  mutate(P_Value = signif(P_Value, 6)) %>%
  arrange(Country)

outfile_between <- file.path(outdir_country, "Between_H_vs_M_propTest_summary.txt")
write.table(between_results, file = outfile_between, sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n✅ Between Historical vs Modern results saved to:\n", outfile_between, "\n")







# do binomial and fisher test per run of 1000 bootstrap and plot the distribution:

# =========================================================
# --- Bootstrap 1000× Binomial + Fisher tests with full outputs ---
# =========================================================
library(dplyr)
library(ggplot2)

set.seed(123)

# ---- Parameters ----
n_boot <- 1000
sample_size <- 8   # e.g. Germany historical n = 8
country_name <- "Germany"

# ---- Historical observed counts ----
hist_plus  <- df_hist %>% filter(country == country_name, Oantigen == "O-antigen+") %>% nrow()
hist_minus <- df_hist %>% filter(country == country_name, Oantigen == "O-antigen-") %>% nrow()
n_hist     <- hist_plus + hist_minus

# ---- Storage ----
results <- data.frame(
  Iteration = integer(),
  Plus_Modern = integer(),
  Minus_Modern = integer(),
  Binomial_p = numeric(),
  Fisher_p = numeric()
)

# ---- Bootstrap loop ----
for (i in seq_len(n_boot)) {
  sampled <- sample(modern_expanded$Isolates, sample_size, replace = TRUE)
  plus  <- sum(sampled == "O-antigen+")
  minus <- sum(sampled == "O-antigen-")
  
  # Binomial test (within modern sample)
  binom_p <- binom.test(plus, sample_size, p = 0.5)$p.value
  
  # Fisher test (historical vs bootstrap modern)
  mat <- matrix(c(hist_plus, hist_minus, plus, minus), nrow = 2, byrow = TRUE)
  fisher_p <- fisher.test(mat)$p.value
  
  results <- rbind(results, data.frame(
    Iteration = i,
    Plus_Modern = plus,
    Minus_Modern = minus,
    Binomial_p = binom_p,
    Fisher_p = fisher_p
  ))
}

# =========================================================
# --- Summary statistics ---
# =========================================================
summary_df <- data.frame(
  Test = c("Binomial", "Fisher"),
  Mean_p = c(mean(results$Binomial_p), mean(results$Fisher_p)),
  Median_p = c(median(results$Binomial_p), median(results$Fisher_p)),
  Prop_p_less_0.05 = c(mean(results$Binomial_p < 0.05), mean(results$Fisher_p < 0.05))
)

# ---- Create output dir ----
outdir <- paste0("bootstrap_pvalue_results_", country_name, "_n", sample_size)
if (!dir.exists(outdir)) dir.create(outdir)

# ---- Write summary + detailed results ----
write.table(summary_df,
            file = file.path(outdir, paste0(country_name, "_pvalue_summary.txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(results,
            file = file.path(outdir, paste0(country_name, "_pvalues_detailed.txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("✅ Saved summary and detailed outputs to", outdir, "\n")

# =========================================================
# --- Plot p-value distributions (separate for each test) ---
# =========================================================

plot_pvalue_dist <- function(pvals, test_name, prop_text, outfile) {
  df_plot <- data.frame(Pvalue = pvals)
  
  p <- ggplot(df_plot, aes(x = Pvalue)) +
    geom_histogram(binwidth = 0.02, fill = "#6699cc", color = "white", alpha = 0.8) +
    geom_vline(xintercept = 0.05, color = "red", linetype = "dashed", linewidth = 1) +
    annotate("text", x = 0.6, y = max(table(cut(pvals, seq(0, 1, 0.02)))) * 0.9,
             label = paste0("Proportion p<0.05 = ", round(prop_text * 100, 1), "%"),
             color = "black", size = 5, fontface = "bold") +
    scale_x_continuous(breaks = seq(0, 1, 0.1)) +
    labs(
      title = paste0(country_name, " (n=", sample_size, ") — ", test_name, " test"),
      x = "p-value",
      y = "Count"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15)
    )
  
  ggsave(outfile, plot = p, width = 7, height = 5, dpi = 600)
  message("✅ Saved: ", outfile)
}

# ---- Plot each test separately ----
plot_pvalue_dist(
  results$Binomial_p,
  "Binomial",
  mean(results$Binomial_p < 0.05),
  file.path(outdir, paste0(country_name, "_Binomial_pvalue_distribution.pdf"))
)

plot_pvalue_dist(
  results$Fisher_p,
  "Fisher",
  mean(results$Fisher_p < 0.05),
  file.path(outdir, paste0(country_name, "_Fisher_pvalue_distribution.pdf"))
)

# =========================================================
# --- Done ---
# =========================================================
cat("✅ Completed bootstrap analysis for", country_name,
    "(", sample_size, " modern resamples, ", n_boot, " iterations )\n")


# =========================================================
# --- Extended summary including mean & SD of sampled counts ---
# =========================================================
summary_df <- data.frame(
  Test = c("Binomial", "Fisher"),
  Mean_p = c(mean(results$Binomial_p), mean(results$Fisher_p)),
  Median_p = c(median(results$Binomial_p), median(results$Fisher_p)),
  Prop_p_less_0.05 = c(mean(results$Binomial_p < 0.05),
                       mean(results$Fisher_p < 0.05)),
  Mean_Plus_Modern = mean(results$Plus_Modern),
  SD_Plus_Modern   = sd(results$Plus_Modern),
  Mean_Minus_Modern = mean(results$Minus_Modern),
  SD_Minus_Modern   = sd(results$Minus_Modern)
)

# Round for cleaner table
summary_df <- summary_df %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

# ---- Save summary ----
outfile_summary <- file.path(outdir, paste0(country_name, "_pvalue_summary.txt"))
write.table(summary_df, file = outfile_summary, sep = "\t", quote = FALSE, row.names = FALSE)

cat("✅ Saved extended summary with mean±SD counts to:\n", outfile_summary, "\n")
