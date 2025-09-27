# Load libraries
library(dplyr)
library(ggplot2)
library(binom)
library(readr)

# === Load data ===
modern_df <- read.table(
  "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers/newsummarymodern57_tailocin_kmer_propnorm.tsv",
  header = TRUE, sep = "\t"
)

# Remove modern samples so 53 modern samples are shown
samples_to_remove <- c("p12.F2", "p13.C7", "p6.A10", "p9.C4")
modern_df <- modern_df[!modern_df$Isolate %in% samples_to_remove, ]

# Historical
historical_df <- read.table(
  "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers/newsummaryhistorical40_tailocin_kmer_propnorm.tsv",
  header = TRUE, sep = "\t"
)
historical_df <- historical_df[!(historical_df$Isolate == "64.GBR_1933b_S36"), ]

# Label isolate type
modern_df$IsolateType <- "Modern"
historical_df$IsolateType <- "Historical"

# Combine and exclude 6 nonOTU5 and 4 samples that have no HTF assigned, as well the 64GBR which has conflicting assignment
df_all <- bind_rows(modern_df, historical_df)
base_exclude <- c("HB0828", "HB0863", "PL0066", "PL0108", "PL0203", "PL0258", "PL0027", "PL0053", "30.ESP_1983b", "PL0026")
excluded_samples <- c(base_exclude, paste0(base_exclude, "|NA"))
df_all <- df_all %>% filter(!(Isolate %in% excluded_samples))

# Annotate HTF length group
df_all <- df_all %>%
  mutate(LengthGroup = case_when(
    Reference == "HTF_p7.G11" ~ "1830",
    Reference == "HTF_p25.A12" ~ "1383",
    Reference %in% c("HTF_p5.D5", "HTF_p23.B8", "HTF_p26.D6", "HTF_p25.C2") ~ "1803",
    Reference == "HTF_p21.F9" ~ "1245",
    TRUE ~ "Unknown"
  )) %>%
  filter(LengthGroup != "Unknown")

# Get dominant HTF per isolate
dominant_df <- df_all %>%
  group_by(Isolate) %>%
  filter(PropNorm == max(PropNorm, na.rm = TRUE)) %>%
  ungroup()

# Define O-antigen status
dominant_df <- dominant_df %>%
  mutate(OAntigenStatus = case_when(
    LengthGroup %in% c("1383", "1830") ~ "HTF group without O-antigen",
    LengthGroup %in% c("1245", "1803") ~ "HTF group with O-antigen",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(OAntigenStatus))

# === Chi-squared test ===
contingency_table <- table(dominant_df$IsolateType, dominant_df$OAntigenStatus)
print(contingency_table)

chisq_result <- chisq.test(contingency_table)
sink("/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step4_additionalanalysis/step3_HTFfreq/chisquare_result.txt")
cat("Chi-squared test p-value:", chisq_result$p.value, "\n")
sink()
# === Proportions for plotting ===
total_counts <- dominant_df %>%
  distinct(Isolate, IsolateType) %>%
  count(IsolateType) %>%
  rename(TotalN = n)

freq_oa_group <- dominant_df %>%
  group_by(OAntigenStatus, IsolateType) %>%
  summarise(Frequency = n(), .groups = "drop") %>%
  left_join(total_counts, by = "IsolateType") %>%
  mutate(
    Proportion = Frequency / TotalN,
    ci = binom.confint(Frequency, TotalN, method = "wilson"),
    Lower = ci$lower,
    Upper = ci$upper
  ) %>%
  select(OAntigenStatus, IsolateType, Frequency, TotalN, Proportion, Lower, Upper)

# === PLOT: Proportion of O-antigen by group ===
freq_oa_group$IsolateType <- factor(freq_oa_group$IsolateType, levels = c("Historical", "Modern"))

p <- ggplot(freq_oa_group, aes(
  x = IsolateType,
  y = Proportion,
  fill = OAntigenStatus
)) +
  geom_point(
    aes(color = OAntigenStatus),
    shape = 21,
    size = 5,
    stroke = 1,
    position = position_dodge(width = 0.6)
  ) +
  geom_errorbar(
    aes(ymin = Lower, ymax = Upper, group = OAntigenStatus),
    width = 0.2,
    position = position_dodge(width = 0.6),
    color = "grey30"
  ) +
  geom_text(
    aes(label = paste0(Frequency, "/", TotalN), group = OAntigenStatus),
    position = position_dodge(width = 0.6),
    vjust = -1.5,
    size = 4.5
  ) +
  scale_fill_manual(
    values = c(
      "HTF group with O-antigen" = "#E69F00",
      "HTF group without O-antigen" = "#56B4E9"
    )
  ) +
  scale_color_manual(
    values = c(
      "HTF group with O-antigen" = "#E69F00",
      "HTF group without O-antigen" = "#56B4E9"
    )
  ) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(
    title = ("HTF Groups in Modern vs Historical dataset"),
    x = NULL,
    y = "Proportion",
    color = "HTF Group",
    fill = "HTF Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12)
  )

ggsave(
  "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step4_additionalanalysis/step3_HTFfreq/step2_1_oantigen_presence_barplot_withchisq.pdf",
  plot = p, width = 8, height = 6
)
cat("✅ Plot saved: oantigen_presence_barplot_withchisq.pdf\n")






#year


# === Add year info ===

# 1. Historical years from text
dates_df <- read_tsv(
  "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step5_additional_analysis_coinfection_Oantigenbreakingisolates_HTFlengthfreq/HTFlengthfreq/allh46_withdatesandlocs_uniq.txt", show_col_types = FALSE)
dates_df <- dates_df %>%
  select(Isolate = `samplename`, Year = year) %>%
  filter(!is.na(Year))

# 2. Modern samples assumed 2018 if isolate starts with 'p'
modern_years <- dominant_df %>%
  filter(IsolateType == "Modern" & grepl("^p", Isolate)) %>%
  mutate(Year = 2018) %>%
  select(Isolate, Year)

# 3. Combine all years
year_df <- bind_rows(dates_df, modern_years)

# Join year into dominant table
dominant_df <- dominant_df %>%
  left_join(year_df, by = "Isolate") %>%
  filter(!is.na(Year))

# Add random Y value for jitter
set.seed(123)  # reproducible jitter
dominant_df$Yjitter <- runif(nrow(dominant_df), min = 0, max = 1)
# === Plot with merged legend ===
p <- ggplot(dominant_df, aes(
  x = Year,
  y = Yjitter,
  shape = OAntigenStatus,
  color = OAntigenStatus
)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_manual(
    values = c(
      "HTF group with O-antigen" = "#E41A1C",
      "HTF group without O-antigen" = "#4DAF4A"
    )
  ) +
  scale_shape_manual(
    values = c(
      "HTF group with O-antigen" = 21,
      "HTF group without O-antigen" = 24
    )
  ) +
  labs(
    title = "Distribution of HTF Groups Over Time",
    x = "Year",
    y = "",
    color = "HTF Group",
    shape = "HTF Group"  # give same name as color
  )  +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "right"
  )

ggsave(
  "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step4_additionalanalysis/step3_HTFfreq/step2_2_1_HTF_oantigen_timeline_jitter.pdf", plot = p, width = 16, height = 3)
cat("✅ Plot saved: HTF_oantigen_timeline_jitter.pdf\n")


write.table(dominant_df,
            file = "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step4_additionalanalysis/step3_HTFfreq/step2_3_HTF_oantigen_dominant_table.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)








library(ggplot2)
library(dplyr)

# === Define colors for the 4 HTF groups ===
length_colors <- c(
  "1830" = "#f6d6ff",
  "1383" = "#638ccc",
  "1803" = "#800233",
  "1245" = "#f9d42a"
)

# === Reorder LengthGroup for y-axis ===
length_groups <- c("1245", "1803", "1383", "1830")

# === Jitter Year for modern samples (narrow spread) ===
set.seed(7)
dominant_df <- dominant_df %>%
  mutate(
    Year = ifelse(IsolateType == "Modern",
                  2018 + runif(n(), -0.2, 0.2),  # jitter around 2018
                  Year),
    LengthGroup = factor(LengthGroup, levels = length_groups),
    Ybase = as.numeric(LengthGroup),
    Yjitter = Ybase + runif(n(), -0.3, 0.3)     # vertical jitter
  )

# === Plot ===
p_htf_timeline <- ggplot(
  dominant_df,
  aes(x = Year, y = Yjitter, color = LengthGroup, alpha = IsolateType)
) +
  geom_point(shape = 16, size = 4) +   # same shape for both, alpha differs
  # Colors: fixed order pink → blue → red → yellow
  scale_color_manual(
    values = length_colors,
    breaks = c("1830","1383","1803","1245")
  ) +
  # Alpha: Modern solid, Historical transparent
  scale_alpha_manual(
    values = c("Historical" = 1, "Modern" = 0.45),guide = "none"   # << suppress legend
  ) +
  facet_wrap(~IsolateType, scales = "free_x") +
  scale_x_continuous(
    breaks = c(1850, 1900, 1950, 2000, 2018),
    labels = c("1850", "1900", "1950", "2000", "2018")
  ) +
  scale_y_continuous(
    breaks = 1:4,
    #labels = length_groups,
    labels = NULL,
    
    expand = expansion(mult = c(0.1, 0.1))
  ) +
  labs(
    #title = "HTF haplotype distribution over time",
    x = "Year",
    #y = "HTF haplotypes (bp)",
    y = NULL,
    color = "HTF haplotypes",
    alpha = "Isolate type"
  ) +
  theme_classic(base_size = 24) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 23),
    #plot.title = element_text(hjust = 0.5, face = "bold")
  )

# === Save figure ===
ggsave(
  "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step4_additionalanalysis/step3_HTFfreq/step2_2_2_HTF_length_timeline_HTFgrouped_alpha.pdf",
  plot = p_htf_timeline,
  width = 20, height = 4
)
