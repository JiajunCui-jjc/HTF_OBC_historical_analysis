#!/usr/bin/env Rscript
# ===========================================================
# Nature-style HTF & O-antigen country summary (v4)
# Symmetric 95 % CI, baseline n=0, proper denominators,
# top-4 + all-country variants, auto color index
# ===========================================================

library(dplyr)
library(ggplot2)
library(readr)
library(binom)
library(patchwork)
library(tidyr)

# -----------------------------------------------------------
# Input paths
# -----------------------------------------------------------
samples_path <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step5_additional_analysis_coinfection_Oantigenbreakingisolates_HTFlengthfreq/HTFlengthfreq/h49PL0087NAtmpm83_withdatesandlocs_uniq copy.txt"
dominant_path <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step4_additionalanalysis/step3_HTFfreq/step2_3_HTF_oantigen_dominant_table.tsv"
dominant_df <- dominant_df %>%
  mutate(Isolate = ifelse(grepl("^p", Isolate),
                          gsub("_", ".", Isolate),
                          Isolate))
outdir <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step4_additionalanalysis/step3_HTFfreq/country_boxplot"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
setwd(outdir)

# -----------------------------------------------------------
# Load + clean
# -----------------------------------------------------------
samples <- read.table(samples_path, header = TRUE, sep = "\t")
dominant_df <- read.table(dominant_path, header = TRUE, sep = "\t")

samples$country[samples$country == "Congo"] <- "Democratic Republic of the Congo"
samples$country[samples$country == "United Kingdom"] <- "UK"

dominant_df <- dominant_df %>%
  mutate(LengthGroup = case_when(
    Reference == "HTF_p7.G11" ~ "1830",
    Reference == "HTF_p25.A12" ~ "1383",
    Reference %in% c("HTF_p5.D5","HTF_p23.B8","HTF_p26.D6","HTF_p25.C2") ~ "1803",
    Reference == "HTF_p21.F9" ~ "1245"
  ),
  Oantigen = case_when(
    LengthGroup %in% c("1803","1245") ~ "O-antigen+",
    LengthGroup %in% c("1830","1383") ~ "O-antigen-"
  )) %>%
  filter(IsolateType == "Historical", !is.na(LengthGroup))

df_hist <- left_join(samples,
                     dominant_df[, c("Isolate","LengthGroup","IsolateType","Oantigen")],
                     by = c("samplename"="Isolate")) %>%
  filter(IsolateType == "Historical", !is.na(LengthGroup))
cat("✅ Historical isolates with HTF:", nrow(df_hist), "\n")

# -----------------------------------------------------------
# Colors
# -----------------------------------------------------------
length_colors <- c("1830"="#f6d6ff","1383"="#638ccc",
                   "1803"="#800233","1245"="#f9d42a")
oantigen_colors <- c("O-antigen-"="#9e9ac8","O-antigen+"="#c94f7c")

# -----------------------------------------------------------
# Theme
# -----------------------------------------------------------
nature_theme <- theme_minimal(base_size = 15) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.text.x = element_text(angle=30,hjust=1,face="bold"),
    axis.text.y = element_text(face="bold"),
    axis.title = element_text(face="bold",size=15),
    plot.title = element_text(face="bold",hjust=0.5,size=16),
    legend.position="none",
    plot.margin=margin(20,20,20,20)
  )

# -----------------------------------------------------------
# Helper: symmetric 95 % CI on counts
# -----------------------------------------------------------
get_country_ci <- function(df, group_cols){
  df %>%
    count(across(all_of(group_cols)), name="Count") %>%
    group_by(across(all_of(setdiff(group_cols,"country")))) %>%
    mutate(Total = sum(Count)) %>%
    ungroup() %>%
    mutate(ci = binom.confint(Count, Total, method="wilson"),
           Lower = ci$lower * Total,
           Upper = ci$upper * Total)
}

# ===========================================================
# 1️⃣ All countries (simple count)
# ===========================================================
country_counts <- df_hist %>% count(country, name="Count") %>% arrange(desc(Count))
p1 <- ggplot(country_counts, aes(x=reorder(country,-Count), y=Count)) +
  geom_point(size=4, color="black") +
  geom_text(aes(label=Count), vjust=-1, size=4.5, fontface="bold") +
  labs(x="", y="Isolate count", title="Historical isolates with HTF by country") +
  nature_theme +
  expand_limits(y = max(country_counts$Count) * 1.15)
ggsave("1_allcountries_points.pdf", p1, width=8, height=5, dpi=600)

# ===========================================================
# 2️⃣ Top 4 countries (simple)
# ===========================================================
top4 <- head(country_counts$country, 4)
top3 <- head(country_counts$country, 3)

p2 <- ggplot(filter(country_counts, country %in% top4),
             aes(x=reorder(country,-Count), y=Count)) +
  geom_point(size=4, color="black") +
  geom_text(aes(label=Count), vjust=-1, size=4.5, fontface="bold") +
  labs(x="", y="Isolate count", title="Top 4 countries (historical HTF isolates)") +
  nature_theme +
  expand_limits(y = max(country_counts$Count) * 1.15)
ggsave("2_top4countries_points.pdf", p2, width=6, height=4, dpi=600)

# ===========================================================
# 3️⃣ HTF haplotypes × top 4 countries (fixed CI)
# ===========================================================
df_htf_top4 <- df_hist %>%
  filter(country %in% top4) %>%
  count(country, LengthGroup, name="Count") %>%
  group_by(LengthGroup) %>%
  mutate(Total = sum(Count)) %>%
  ungroup() %>%
  mutate(ci = binom.confint(Count, Total, method="wilson"),
         Lower = ci$lower * Total,
         Upper = ci$upper * Total)

p3 <- ggplot(df_htf_top4, aes(x=reorder(country,-Count), y=Count, color=LengthGroup)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0.15, size=0.7) +
  geom_point(size=4) +
  geom_text(aes(label=Count), vjust=-1, size=4, fontface="bold", color="black") +
  facet_wrap(~LengthGroup, ncol=2, scales="free_y") +
  scale_color_manual(values=length_colors) +
  labs(x="", y="Isolate count", title="HTF haplotype distribution (top 4 countries)") +
  nature_theme +
  expand_limits(y = max(df_htf_top4$Upper, na.rm=TRUE) * 1.15)
ggsave("3_top4countries_HTF_panels_points.pdf", p3, width=8, height=6, dpi=600)

# ===========================================================
# 4️⃣ HTF haplotypes × all countries
# ===========================================================
df_all_htf <- get_country_ci(df_hist, c("country","LengthGroup"))
p4 <- ggplot(df_all_htf, aes(x=reorder(country,-Count), y=Count, color=LengthGroup)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0.15, size=0.7) +
  geom_point(size=3.5) +
  geom_text(aes(label=Count), vjust=-1, size=3.5, fontface="bold", color="black") +
  facet_wrap(~LengthGroup, ncol=2, scales="free_y") +
  scale_color_manual(values=length_colors) +
  labs(x="", y="Isolate count", title="HTF haplotype distribution (all countries)") +
  nature_theme +
  expand_limits(y = max(df_all_htf$Upper, na.rm=TRUE) * 1.15)
ggsave("4_allcountries_HTF_panels_points.pdf", p4, width=10, height=6, dpi=600)

# ===========================================================
# 5️⃣ O-antigen ± × all countries
# ===========================================================
df_oant <- get_country_ci(df_hist, c("country","Oantigen"))
p5 <- ggplot(df_oant, aes(x=reorder(country,-Count), y=Count, color=Oantigen)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0.15, size=0.7) +
  geom_point(size=4) +
  geom_text(aes(label=Count), vjust=-1, fontface="bold", size=4.5, color="black") +
  facet_wrap(~Oantigen, ncol=2, scales="free_y") +
  scale_color_manual(values=oantigen_colors) +
  labs(x="", y="Isolate count", title="O-antigen group distribution (all countries)") +
  nature_theme +
  expand_limits(y = max(df_oant$Upper, na.rm=TRUE) * 1.15)
ggsave("5_Oantigen_groups_allcountries_NatureStyle.pdf", p5, width=15, height=5, dpi=600)



# ===========================================================
# 6️⃣ O-antigen ± × top 4 countries — proportion (n = 38 total)
# ===========================================================
library(ggplot2)
library(dplyr)
library(binom)
library(scales)

total_hist <- 38  # total historical isolates with HTF

df_oant_top4_prop <- df_hist %>%
  filter(country %in% top4) %>%
  count(country, Oantigen, name = "Count") %>%
  mutate(
    Total = total_hist,
    Proportion = Count / Total,
    ci = binom.confint(Count, Total, method = "wilson"),
    Lower = ci$lower,
    Upper = ci$upper
  )

# ---- Plot with percentages only ----
dodge <- position_dodge(width = 0.8)

p6_prop <- ggplot(df_oant_top4_prop,
                  aes(x = reorder(country, -Proportion),
                      y = Proportion,
                      color = Oantigen)) +
  geom_point(position = dodge, size = 3, show.legend = TRUE) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                width = 0.1, position = dodge, show.legend = FALSE) +
  geom_text(aes(label = paste0(round(Proportion * 100, 1), "%"),
                group = Oantigen),
            position = dodge,
            vjust = -1.0, size = 4, fontface = "bold",
            color = "black", show.legend = FALSE) +
  facet_wrap(~Oantigen, ncol = 2, scales = "free_y") +
  scale_color_manual(values = c("O-antigen-" = "#9e9ac8",
                                "O-antigen+" = "#c94f7c")) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, max(df_oant_top4_prop$Upper, na.rm = TRUE) * 1.15),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = "", y = "Proportion (%)") +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 11, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.y.right = element_blank(),
    axis.title.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right = element_blank()
  )

ggsave("7_Oantigen_groups_top4countries_Proportion_n38_classic.pdf",
       plot = p6_prop, width = 9, height = 3.5, dpi = 600)


# ===========================================================
# 6️⃣ O-antigen ± × top 4 countries — merged single panel
# ===========================================================
library(ggplot2)
library(dplyr)
library(binom)
library(scales)

total_hist <- 38  # total historical isolates with HTF

df_oant_top4_prop <- df_hist %>%
  filter(country %in% top4) %>%
  count(country, Oantigen, name = "Count") %>%
  mutate(
    Total = total_hist,
    Proportion = Count / Total,
    ci = binom.confint(Count, Total, method = "wilson"),
    Lower = ci$lower,
    Upper = ci$upper
  )

# ---- Plot with both O-antigen ± in one panel ----
dodge <- position_dodge(width = 0.6)

p6_prop_merged <- ggplot(df_oant_top4_prop,
                         aes(x = reorder(country, -Proportion),
                             y = Proportion,
                             color = Oantigen)) +
  geom_point(position = dodge, size = 3, show.legend = TRUE) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                width = 0.15, position = dodge, show.legend = FALSE) +
  geom_text(aes(label = paste0(round(Proportion * 100, 1), "%"),
                group = Oantigen),
            position = position_dodge(width = 0.6),
            vjust = -0.8, size = 4, fontface = "bold",
            color = "black", show.legend = FALSE) +
  scale_color_manual(values = c("O-antigen-" = "#9e9ac8",
                                "O-antigen+" = "#c94f7c")) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, max(df_oant_top4_prop$Upper, na.rm = TRUE) * 1.15),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = "", y = "Proportion (%)") +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 11, face = "bold")
  )

ggsave("8_Oantigen_groups_top4countries_Proportion_merged_n38.pdf",
       plot = p6_prop_merged, width = 9, height = 4, dpi = 600)

# ===========================================================
# 8️⃣ O-antigen ± × top 4 countries — merged single panel (m + h)
# ===========================================================
library(ggplot2)
library(dplyr)
library(binom)
library(scales)
dominant_df <- read.table(dominant_path, header = TRUE, sep = "\t")

dominant_df <- dominant_df %>%
  mutate(LengthGroup = case_when(
    Reference == "HTF_p7.G11" ~ "1830",
    Reference == "HTF_p25.A12" ~ "1383",
    Reference %in% c("HTF_p5.D5","HTF_p23.B8","HTF_p26.D6","HTF_p25.C2") ~ "1803",
    Reference == "HTF_p21.F9" ~ "1245"
  ),
  Oantigen = case_when(
    LengthGroup %in% c("1803","1245") ~ "O-antigen+",
    LengthGroup %in% c("1830","1383") ~ "O-antigen-"
  )) 
samples <- read.table(samples_path, header = TRUE, sep = "\t")
samples$country[samples$country == "United Kingdom"] <- "UK"
samples$samplename <- ifelse(grepl("^p", samples$samplename),
                          gsub("_", ".", as.character(samples$samplename)),
                          as.character(samples$samplename))


# include both historical + modern
df_oant_all <- right_join(
  samples[, c("samplename", "country")],
  dominant_df[, c("Isolate", "LengthGroup", "IsolateType", "Oantigen")],
  by = c("samplename" = "Isolate")
)

cat("✅ Total isolates with HTF (m + h):", nrow(df_oant_all), "\n")

total_all <- nrow(df_oant_all)

# Use same top 4 countries (based on historical)
df_oant_top4_prop_all <- df_oant_all %>%
  filter(country %in% top4) %>%
  count(country, Oantigen, name = "Count") %>%
  mutate(
    Total = total_all,
    Proportion = Count / Total,
    ci = binom.confint(Count, Total, method = "wilson"),
    Lower = ci$lower,
    Upper = ci$upper
  )

# ---- Plot both O-antigen ± in one panel ----
dodge <- position_dodge(width = 0.6)

p8_prop_merged_all <- ggplot(df_oant_top4_prop_all,
                             aes(x = reorder(country, -Proportion),
                                 y = Proportion,
                                 color = Oantigen)) +
  geom_point(position = dodge, size = 3, show.legend = TRUE) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                width = 0.15, position = dodge, show.legend = FALSE) +
  geom_text(aes(label = paste0(round(Proportion * 100, 1), "%"),
                group = Oantigen),
            position = position_dodge(width = 0.6),
            vjust = -0.8, size = 4, fontface = "bold",
            color = "black", show.legend = FALSE) +
  scale_color_manual(values = c("O-antigen-" = "#9e9ac8",
                                "O-antigen+" = "#c94f7c")) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, max(df_oant_top4_prop_all$Upper, na.rm = TRUE) * 1.15),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = "", y = "Proportion (%)",
       title = paste0("O-antigen ± × top 4 countries (m + h, n=", total_all, ")")) +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 11, face = "bold")
  )

ggsave("9mandh_Oantigen_groups_top4countries_Proportion_merged_mh.pdf",
       plot = p8_prop_merged_all, width = 9, height = 4, dpi = 600)
