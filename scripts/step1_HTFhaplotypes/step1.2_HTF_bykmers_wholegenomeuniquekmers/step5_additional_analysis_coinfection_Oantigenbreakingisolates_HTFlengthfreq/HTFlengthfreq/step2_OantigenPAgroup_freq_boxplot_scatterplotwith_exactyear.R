# Load libraries
library(dplyr)
library(ggplot2)
library(binom)
library(readr)

# === Load data ===
# Modern

# ================================
# 1️⃣ Load modern + historical prop tables
# ================================
# Modern
modern_df <- read.table("/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers/newsummarymodern57_tailocin_kmer_propnorm.tsv" ,header = TRUE, sep = "\t"
)
#rm  p12.F2 p13.C7 p6.A10 p9.C4
# List of samples to remove
samples_to_remove <- c("p12.F2", "p13.C7", "p6.A10", "p9.C4")

# Filter out those samples
modern_df <- modern_df[!modern_df$Isolate %in% samples_to_remove, ]
# Historical
historical_df <- read.table(
  "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers/newsummaryhistorical40_tailocin_kmer_propnorm.tsv",
  header = TRUE, sep = "\t"
)

#rm historical_df$Isolate=="64.GBR_1933b_S36" & historical_df$Reference=="HTF_p21.F9"
#for now, we believe 64 is unknown
historical_df <- historical_df[!(historical_df$Isolate == "64.GBR_1933b_S36"), ]
# Label isolate type
modern_df$IsolateType <- "Modern"
historical_df$IsolateType <- "Historical"

# Combine and exclude problematic samples
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

# === Set aesthetics ===
shape_vals <- c(
  "HTF group with O-antigen" = 21,
  "HTF group without O-antigen" = 24
)
color_vals <- c(
  "Modern" = "yellow",        # blue
  "Historical" = "darkgreen"     # red
)

# Count and compute confidence intervals
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


# === PLOT 1: O-antigen proportion barplot ===

# Plot
p1 <- ggplot(freq_oa_group, aes(
  x = OAntigenStatus,
  y = Proportion,
  shape = OAntigenStatus,
  fill = IsolateType
)) +
  geom_point(
    aes(color = IsolateType),
    size = 6,
    position = position_dodge(width = 0.6),
    stroke = 1
  ) +
  geom_errorbar(
    aes(ymin = Lower, ymax = Upper),
    width = 0.2,
    position = position_dodge(width = 0.6)
  ) +
  geom_text(
    aes(label = paste0(Frequency, "/", TotalN)),
    position = position_dodge(width = 0.6),
    vjust = -1.5,
    size = 5
  ) +
  scale_shape_manual(values = shape_vals) +
  scale_color_manual(values = color_vals) +
  scale_fill_manual(values = color_vals) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(
    title = "Proportion of Isolates by O-antigen Status",
    x = "O-antigen Group",
    y = "Proportion",
    shape = "O-antigen Group",
    color = "Isolate Type",
    fill = "Isolate Type"
  ) +
  theme_minimal(base_size = 14)

ggsave(
  "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step4_additionalanalysis/step3_HTFfreq/step2_1_oantigen_presence_barplot_updated.pdf",
  plot = p1, width = 8, height = 6
)
cat("✅ Plot saved: oantigen_presence_barplot.pdf\n")





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
"/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step4_additionalanalysis/step3_HTFfreq/step2_2_HTF_oantigen_timeline_jitter.pdf", plot = p, width = 16, height = 3)
cat("✅ Plot saved: HTF_oantigen_timeline_jitter.pdf\n")


write.table(dominant_df,
            file = "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step4_additionalanalysis/step3_HTFfreq/step2_3_HTF_oantigen_dominant_table.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

