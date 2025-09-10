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

# Remove problematic modern samples
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

# Combine and exclude more problematic samples
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