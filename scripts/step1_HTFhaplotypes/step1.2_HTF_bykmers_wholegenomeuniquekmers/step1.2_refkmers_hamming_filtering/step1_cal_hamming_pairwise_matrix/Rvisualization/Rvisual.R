library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)

# === 1️⃣ Set working directory ===
setwd("/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step1_cal_hamming_pairwise_matrix/Rvisualization")
#mkdir -p /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step1_cal_hamming_pairwise_matrix/Rvisualization
# === 2️⃣ Load data ===
df <- read.table("../pairwise_hamming_distance_longformatbypy.tsv", header=TRUE, sep="\t")

# ✅ Check data
print(head(df))
# === 3️⃣ Per-pair distribution plot (log10 y-axis) ===
p_perpair <- df %>%
  group_by(QueryRef, TargetRef, HammingDist) %>%
  summarise(Count=n(), .groups="drop") %>%
  group_by(QueryRef, TargetRef) %>%
  mutate(Proportion = Count/sum(Count)) %>%
  ggplot(aes(x=HammingDist, y=Proportion, fill=QueryRef)) +
  geom_bar(stat="identity", position="dodge") +
  scale_y_sqrt() +  # ✅ square-root scale
  facet_grid(QueryRef ~ TargetRef) +
  labs(title="Hamming distance per reference pair (sqrt y-scale)",
       x="Hamming distance", y="Proportion (sqrt)") +
  theme_minimal()

ggsave("hamming_distance_perpair_sqrt.pdf", p_perpair, width=12, height=8)

# === 4️⃣ Overall distribution plot ===
p_total <- df %>%
  group_by(HammingDist) %>%
  summarise(Count=n(), .groups="drop") %>%
  mutate(Proportion = Count/sum(Count)) %>%
  ggplot(aes(x=HammingDist, y=Proportion)) +
  geom_bar(stat="identity", fill="steelblue", color="black") +
  labs(title="Overall hamming distance", x="Hamming distance", y="Proportion") +
  theme_minimal()

ggsave("hamming_distance_total.pdf", p_total, width=8, height=6)

# === 5️⃣ Minimum distance matrix ===
min_dist_table <- df %>%
  group_by(QueryRef, TargetRef) %>%
  summarise(MinDistance = min(HammingDist), .groups="drop") %>%
  pivot_wider(names_from=TargetRef, values_from=MinDistance)

write.table(min_dist_table, "hamming_distance_min_matrix.tsv", sep="\t", row.names=FALSE, quote=FALSE)

# === 6️⃣ Median distance matrix ===
median_dist_table <- df %>%
  group_by(QueryRef, TargetRef) %>%
  summarise(MedianDistance = median(HammingDist), .groups="drop") %>%
  pivot_wider(names_from=TargetRef, values_from=MedianDistance)

write.table(median_dist_table, "hamming_distance_median_matrix.tsv", sep="\t", row.names=FALSE, quote=FALSE)

# ✅ Print summary
cat("✅ Minimum and median distance tables saved.\n")

# === 7️⃣ Convert to matrix for heatmap ===
rownames_min <- min_dist_table$QueryRef
min_matrix <- as.matrix(min_dist_table[,-1])
rownames(min_matrix) <- rownames_min

rownames_median <- median_dist_table$QueryRef
median_matrix <- as.matrix(median_dist_table[,-1])
rownames(median_matrix) <- rownames_median

# === 8️⃣ Plot heatmaps ===
pheatmap(min_matrix,
         display_numbers = TRUE,
         number_format = "%.0f",
         color = colorRampPalette(c("white", "orange", "red"))(50),
         main = "Minimum Hamming Distance",
         filename = "hamming_distance_min_heatmap.pdf",
         width=8, height=8)

pheatmap(median_matrix,
         display_numbers = TRUE,
         number_format = "%.0f",
         color = colorRampPalette(c("white", "lightblue", "blue"))(50),
         main = "Median Hamming Distance",
         filename = "hamming_distance_median_heatmap.pdf",
         width=8, height=8)

cat("✅ Heatmaps saved.\n")
