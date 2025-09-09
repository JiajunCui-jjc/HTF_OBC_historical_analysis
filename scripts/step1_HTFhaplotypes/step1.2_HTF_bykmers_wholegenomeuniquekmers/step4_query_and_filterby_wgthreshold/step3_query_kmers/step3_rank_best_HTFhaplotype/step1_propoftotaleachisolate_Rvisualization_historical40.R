library(ggplot2)
library(dplyr)

# --- Set working directory ---
#setwd('/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers')
setwd('/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers')
# --- Input file ---
summary_file <- "summaryhistorical40_tailocin_after0.65wgavgdepth_filter.tsv" 

# --- Samples to exclude ---
base_exclude <- c("HB0828", "HB0863", "PL0066", "PL0108", "PL0203", 
                  "PL0258", "PL0027", "PL0053", "PL0065", "30.ESP_1983b", "PL0026")
excluded_samples <- c(base_exclude, paste0(base_exclude, "|NA"))

# --- Load and filter ---

# --- Load and filter ---
df_all <- read.table(summary_file, header = TRUE, sep = "\t")
df_all <- df_all[!df_all$Isolate %in% excluded_samples, ]

# ✅ Compute total matched kmers per isolate
df_total <- df_all %>%
  group_by(Isolate) %>%
  summarise(TotalMatched = sum(Matched_Kmers, na.rm = TRUE), .groups = "drop")

# ✅ Normalize matched kmers
df_all <- df_all %>%
  left_join(df_total, by = "Isolate") %>%
  mutate(PropNorm = ifelse(TotalMatched > 0, Matched_Kmers / TotalMatched, 0))

# ✅ Identify best HTF match per isolate
best_refs <- df_all %>%
  group_by(Isolate) %>%
  slice_max(order_by = PropNorm, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(BestHTF = Reference)

# ✅ Assign clean HTF label
best_refs <- best_refs %>%
  mutate(HTF_by_kmer = case_when(
    BestHTF == "HTF_p7.G11"  ~ "HTF_p7.G11 (1830)",
    BestHTF == "HTF_p25.A12" ~ "HTF_p25.A12 (1383)",
    BestHTF == "HTF_p25.C2"  ~ "HTF_p25.C2 (1803)",
    BestHTF == "HTF_p5.D5"   ~ "HTF_p5.D5 (1803)",
    BestHTF == "HTF_p26.D6"  ~ "HTF_p26.D6 (1803)",
    BestHTF == "HTF_p23.B8"  ~ "HTF_p23.B8 (1803)",
    BestHTF == "HTF_p21.F9"  ~ "HTF_p21.F9 (1245)",
    TRUE ~ BestHTF
  ))

# ✅ Save clean HTF label per sample
write.table(
  best_refs %>% select(Isolate, HTF_by_kmer) %>% arrange(HTF_by_kmer, Isolate),
  "h40sample_bestHTF.txt",
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
)

# ✅ Join best HTF info back
df_all <- df_all %>%
  left_join(best_refs[, c("Isolate", "BestHTF")], by = "Isolate") %>%
  mutate(IsolateLabel = paste(Isolate, BestHTF, sep = " | "))

# ✅ Set factor level to order bars by BestHTF, then max PropNorm
ordering <- df_all %>%
  group_by(IsolateLabel, BestHTF) %>%
  summarise(BestPropNorm = max(PropNorm), .groups = "drop") %>%
  arrange(BestHTF, desc(BestPropNorm)) %>%
  pull(IsolateLabel)

df_all$IsolateLabel <- factor(df_all$IsolateLabel, levels = ordering)

# ✅ Plot
pdf("historical40_stacked_PropNorm_byBestHTF.pdf", width = 11, height = 5)
ggplot(df_all, aes(x = IsolateLabel, y = PropNorm, fill = Reference)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    legend.position = "right"
  ) +
  labs(title = "Normalized K-mer Proportions Ordered by Best HTF Match",
       x = "Isolate | Best-Matched HTF", y = "Normalized Proportion",
       fill = "HTF Reference")
dev.off()

# ✅ Save updated summary table
write.table(
  df_all %>%
    select(Isolate, Reference, Total_Allele_Kmers, Matched_Kmers, PropNorm),
  "newsummaryhistorical40_tailocin_kmer_propnorm.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
)

cat("✅ Done:\n- historical40_stacked_PropNorm_byBestHTF.pdf\n- h40sample_bestHTF.txt\n- newsummaryhistorical40_tailocin_kmer_propnorm.tsv\n")
