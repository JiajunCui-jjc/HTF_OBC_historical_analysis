#!/usr/bin/env Rscript

# === Libraries ===
library(ggplot2)
library(dplyr)

# --- Working directory ---
setwd('/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers')
#setwd('/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers')

# --- Input file (two-sided threshold results) ---
summary_file <- "summaryhistorical49_tailocin_05sd.tsv"

# --- Load ---
df_all <- read.table(summary_file, header = TRUE, sep = "\t", fill = TRUE)

# === Compute normalized k-mer proportions ===
# 1️⃣ Total matched kmers per isolate
df_total <- df_all %>%
  group_by(Isolate) %>%
  summarise(TotalMatched = sum(Matched_Kmers, na.rm = TRUE), .groups = "drop")

# 2️⃣ Normalize matched kmers (preserve NA)
df_all <- df_all %>%
  left_join(df_total, by = "Isolate") %>%
  mutate(PropNorm = ifelse(!is.na(Matched_Kmers) & TotalMatched > 0,
                           Matched_Kmers / TotalMatched,
                           NA_real_))

# === Identify best-matching HTF per isolate ===
best_refs <- df_all %>%
  group_by(Isolate) %>%
  summarise(
    BestHTF = {
      if (all(is.na(PropNorm))) {
        NA_character_
      } else {
        idx <- which.max(replace(PropNorm, is.na(PropNorm), -Inf))
        Reference[idx]
      }
    },
    .groups = "drop"
  )

# === Assign clean HTF length labels ===
best_refs <- best_refs %>%
  mutate(HTF_by_kmer = case_when(
    is.na(BestHTF)        ~ NA_character_,
    BestHTF == "HTF_p7.G11"  ~ "HTF_p7.G11 (1830)",
    BestHTF == "HTF_p25.A12" ~ "HTF_p25.A12 (1383)",
    BestHTF == "HTF_p25.C2"  ~ "HTF_p25.C2 (1803)",
    BestHTF == "HTF_p5.D5"   ~ "HTF_p5.D5 (1803)",
    BestHTF == "HTF_p26.D6"  ~ "HTF_p26.D6 (1803)",
    BestHTF == "HTF_p23.B8"  ~ "HTF_p23.B8 (1803)",
    BestHTF == "HTF_p21.F9"  ~ "HTF_p21.F9 (1245)",
    TRUE ~ BestHTF
  ))

# ✅ Save clean best-HTF table
write.table(
  best_refs %>% select(Isolate, HTF_by_kmer) %>% arrange(HTF_by_kmer, Isolate),
  "h49sample_bestHTF_05sd.txt",
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
)

# === Join BestHTF info back ===
df_all <- df_all %>%
  left_join(best_refs %>% select(Isolate, BestHTF), by = "Isolate") %>%
  mutate(IsolateLabel = paste(Isolate,
                              ifelse(is.na(BestHTF), "NA", BestHTF),
                              sep = " | "))

# === Ordering for plotting (NA-safe) ===
ordering <- df_all %>%
  group_by(IsolateLabel, BestHTF) %>%
  summarise(
    BestPropNorm = if (all(is.na(PropNorm))) NA_real_ else max(PropNorm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    BestHTF = ifelse(is.na(BestPropNorm), NA, BestHTF)
  ) %>%
  arrange(BestHTF, desc(BestPropNorm)) %>%
  pull(IsolateLabel)

df_all$IsolateLabel <- factor(df_all$IsolateLabel, levels = ordering)

# === Plot ===
pdf("historical49_stacked_PropNorm_byBestHTF_05sd.pdf", width = 12, height = 5)
ggplot(df_all, aes(x = IsolateLabel, y = PropNorm, fill = Reference)) +
  geom_bar(stat = "identity", na.rm = TRUE) +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    legend.position = "right"
  ) +
  labs(
    title = "Normalized K-mer Proportions (Two-Sided Threshold, 49 Historical Isolates)",
    x = "Isolate | Best-Matched HTF",
    y = "Normalized Proportion",
    fill = "HTF Reference"
  )
dev.off()

# === Save updated summary table ===
write.table(
  df_all %>%
    select(Isolate, Reference, Total_Allele_Kmers, Matched_Kmers, PropNorm),
  "newsummaryhistorical49_tailocin_kmer_propnorm_05sd.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
)

cat("✅ Done:\n",
    "- historical49_stacked_PropNorm_byBestHTF_05sd.pdf\n",
    "- h49sample_bestHTF_05sd.txt\n",
    "- newsummaryhistorical49_tailocin_kmer_propnorm_05sd.tsv\n")

