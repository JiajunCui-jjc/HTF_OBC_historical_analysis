# ============================
# Supp Fig 1: Endogenous composition + Coverage/Depth
# Nature-style: larger readable text (base 16–18)
# ============================

library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(ggtext)
#dir.create("visualization_suppfig1", showWarnings = FALSE)

setwd("/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step1_h49metadata/visualization_suppfig1")
infile <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step1_h49metadata/SuppTable_h49.txt"
df <- read.table(infile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Convert % to numeric fractions
cols_to_clean <- c("At_proportion", "Ps_proportion","At_covered_proportion","Ps_covered_proportion")
df[cols_to_clean] <- lapply(df[cols_to_clean], function(x) as.numeric(sub("%", "", x)) / 100)
# ---- Colors ----
col_at    <- "#d5da6d"
col_ps    <- "#e86b7d"
col_other <- "grey85"

# ============================
# (A) Endogenous composition stacked barplot
# ============================

df_long <- df %>%
  mutate(Other_proportion = 1 - (At_proportion + Ps_proportion)) %>%
  arrange(Ps_proportion) %>%
  mutate(SAMPLE = factor(SAMPLE, levels = rev(SAMPLE))) %>%
  select(SAMPLE, At_proportion, Ps_proportion, Other_proportion) %>%
  pivot_longer(cols = c("Other_proportion","At_proportion","Ps_proportion"),
               names_to = "Species", values_to = "Prop") %>%
  mutate(Species = factor(Species,
                          levels = c("Other_proportion","At_proportion","Ps_proportion"),
                          labels = c("Other",
                                     "*Arabidopsis thaliana*",
                                     "*Pseudomonas* sp.")))

p_a <- ggplot(df_long, aes(x = SAMPLE, y = Prop, fill = Species)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = c("Other" = col_other,
                               "*Arabidopsis thaliana*" = col_at,
                               "*Pseudomonas* sp." = col_ps)) +
  labs(x = "Sample", y = "Proportion of endogenous reads") +
  theme_classic(base_size = 17) +    # ↑ main text size increased
  theme(
    axis.text.x = element_text(size = 13, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_markdown(size = 15, face = "bold")
  )

# ============================
# (B) Covered proportion vs mean depth
# ============================

df_scatter <- df %>%
  select(SAMPLE, At_cov = At_covered_proportion, At_depth=At_average_depth,
         Ps_cov = Ps_covered_proportion, Ps_depth=Ps_average_depth) %>%
  pivot_longer(cols = everything(),
               names_to = c("Species",".value"),
               names_pattern = "(.*)_(.*)") %>%
  mutate(Species = recode(Species,
                          "At"="*Arabidopsis thaliana*",
                          "Ps"="*Pseudomonas* sp.")) %>%
  filter(!is.na(depth) & !is.na(cov))

max_depth_ps <- df_scatter %>% filter(Species == "*Pseudomonas* sp.") %>% pull(depth) %>% max(na.rm = TRUE)
max_depth_at <- df_scatter %>% filter(Species == "*Arabidopsis thaliana*") %>% pull(depth) %>% max(na.rm = TRUE)
xmax_ps <- ceiling(max_depth_ps / 10) * 10 + 10
xmax_at <- ceiling(max_depth_at / 5) * 5

# --- Arabidopsis ---
p_b1 <- df_scatter %>%
  filter(Species == "*Arabidopsis thaliana*") %>%
  ggplot(aes(x = depth, y = cov)) +
  geom_point(color = col_at, size = 4, alpha = 0.9) +
  scale_x_continuous(breaks = seq(0, xmax_at, by = 5), limits = c(0, xmax_at)) +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Mean depth", y = "Covered proportion",
       title = expression(italic("Arabidopsis thaliana"))) +
  theme_classic(base_size = 17) +
  theme(
    plot.title = element_text(size = 17, face = "bold.italic", hjust = 0.5),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14, face = "bold")
  )

# --- Pseudomonas ---
p_b2 <- df_scatter %>%
  filter(Species == "*Pseudomonas* sp.") %>%
  ggplot(aes(x = depth, y = cov)) +
  geom_point(color = col_ps, size = 4, alpha = 0.9) +
  scale_x_continuous(breaks = seq(0, xmax_ps, by = 25), limits = c(0, xmax_ps)) +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Mean depth", y = "Covered proportion",
       title = expression(italic("Pseudomonas")~"sp.")) +
  theme_classic(base_size = 17) +
  theme(
    plot.title = element_text(size = 17, face = "bold.italic", hjust = 0.5),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14, face = "bold")
  )

# ============================
# Combine panels
# ============================

p_final <- (p_a / (p_b1 + p_b2)) +
  plot_annotation(tag_levels = 'A', theme = theme(plot.title = element_text(size=18, face="bold"))) +
  plot_layout(heights = c(1.3, 1))

# ---- Save ----
ggsave("SuppFig1_Endogenous_CoverageDepth_final_Nature.pdf",
       plot = p_final, width = 12, height = 8, dpi = 600, useDingbats = FALSE)

ggsave("SuppFig1a_Endogenous_sortedByPs_Nature.pdf",
       plot = p_a, width = 9, height = 5.5, dpi = 600)

ggsave("SuppFig1b_At_Ps_separate_Nature.pdf",
       plot = p_b1 + p_b2 + plot_layout(ncol = 2),
       width = 11, height = 5, dpi = 600)

cat(sprintf("✅ X-axis auto-set: At ≤ %.1f, Ps ≤ %.1f\n", xmax_at, xmax_ps))
