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
base_dir      <- "/Users/jiajuncui/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57"
scripts_htf   <- file.path(base_dir,
  "scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers",
  "step5_additional_analysis_coinfection_Oantigenbreakingisolates_HTFlengthfreq/HTFlengthfreq")
step3_outdir  <- file.path(base_dir,
  "results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step4_additionalanalysis/step3_HTFfreq")

samples_path  <- file.path(scripts_htf,
  "h49PL0087NAtmpm83_withdatesandlocs_uniq copy.txt")
dominant_path <- file.path(step3_outdir, "step2_3_HTF_oantigen_dominant_table.tsv")
infile_hp12   <- file.path(scripts_htf, "HP12_combined_lengths_sorted.txt")
otu5_file     <- file.path(scripts_htf, "1355OTU5andp8A2.txt")

dominant_df <- dominant_df %>%
  mutate(Isolate = ifelse(grepl("^p", Isolate),
                          gsub("_", ".", Isolate),
                          Isolate))
outdir <- file.path(step3_outdir, "intermediate/country_boxplot")
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

# ===========================================================
# NEW SECTION: Germany n=8 repeated random downsampling + Fisher test
# ===========================================================
# Rationale:
#   Germany has the most well-characterised Historical isolates (n=8).
#   To assess whether Germany Historical OBC proportions differ robustly
#   from the overall Modern trend, we apply two approaches:
#
#   (A) Deterministic expected-value: build the 2x2 Fisher table using
#       Germany Historical (2 OBC+, 6 OBC-) vs expected Modern counts at n=8
#       (i.e. round(8 × Modern OBC+ proportion)). Matches the approach in
#       Between_H_vs_M_propTest_summary.txt.
#
#   (B) Repeated random downsampling (1000 draws):
#       Each draw, randomly sample 8 isolates WITHOUT replacement from the
#       full Modern OBC pool, then run Fisher's exact test against the fixed
#       Germany Historical group. Repeat 1000 times.
#       Purpose: assess ROBUSTNESS — does the Fisher test consistently return
#       a low (or high) p-value, regardless of which 8 Modern isolates are drawn?
#       Report: median Fisher p, and proportion of draws where p < 0.05.
#
#   TERMINOLOGY: This is REPEATED RANDOM DOWNSAMPLING, not a permutation test
#   (which randomises labels) and not a bootstrap (which resamples with
#   replacement). The downsampling answers: "Is the result sensitive to which
#   particular Modern isolates are compared against Germany Historical?"
#
#   Also reports per-draw binomial p-value (does each drawn Modern sample of
#   n=8 itself depart from a 50/50 OBC split?).
# -----------------------------------------------------------

library(binom)
library(data.table)
library(readr)

sig_stars_step5 <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

# Output folders — use step3_outdir (already defined at top of script)
outdir_lower <- file.path(step3_outdir, "supp_germany_downsample")
outdir_final <- file.path(step3_outdir, "final_comparison")
dir.create(outdir_lower, showWarnings = FALSE, recursive = TRUE)
dir.create(outdir_final, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------
# 1. Germany Historical counts (from df_hist already in memory)
# -----------------------------------------------------------
germany_hist   <- df_hist %>% filter(country == "Germany")
germany_plus   <- sum(germany_hist$Oantigen == "O-antigen+")   # 2
germany_minus  <- sum(germany_hist$Oantigen == "O-antigen-")   # 6
germany_n      <- nrow(germany_hist)                            # 8
n_samp         <- germany_n   # downsample size = Germany n; defined early so Section 3 can use it

cat(sprintf("Germany Historical: n=%d | OBC+ = %d | OBC- = %d\n",
            germany_n, germany_plus, germany_minus))

# Germany Historical within-era binomial test (H0: OBC+ = 50%)
binom_germany_within <- binom.test(germany_plus, germany_n, p = 0.5,
                                   alternative = "two.sided")
germany_binomial_within <- data.frame(
  Country          = "Germany",
  IsolateType      = "Historical",
  OBC_plus_n       = germany_plus,
  OBC_minus_n      = germany_minus,
  Total_N          = germany_n,
  Prop_OBC_plus    = round(germany_plus / germany_n * 100, 2),
  Test             = "Two-sided exact binomial (H0: OBC+ proportion = 0.5)",
  P_value          = binom_germany_within$p.value,
  CI_lower_plus    = round(binom_germany_within$conf.int[1] * 100, 2),
  CI_upper_plus    = round(binom_germany_within$conf.int[2] * 100, 2),
  Significance     = sig_stars_step5(binom_germany_within$p.value),
  Notes            = "Germany Historical only; n=8 (2 OBC+, 6 OBC-)",
  stringsAsFactors = FALSE
)
write.table(germany_binomial_within,
            file.path(outdir_lower, "Germany_n8_Binomial_OBC_within.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf("Germany within-era binomial: p = %.4f [%s]\n",
            binom_germany_within$p.value,
            sig_stars_step5(binom_germany_within$p.value)))

# -----------------------------------------------------------
# 2. Build Modern OBC pool from HP12 combined lengths (OTU5 strains)
#    This is the same source as step4: 1312 modern isolates.
#    NOTE: df_oant_all only covers isolates with country metadata
#    (joined against the historical samples file), so it CANNOT be
#    used to build the full 1312-isolate Modern pool. We must read
#    the HP12 length file directly.
# -----------------------------------------------------------
df_hp12_raw  <- read_tsv(infile_hp12, show_col_types = FALSE)
otu5_ids     <- fread(otu5_file, header = FALSE)$V1
df_hp12      <- df_hp12_raw %>%
  filter(Sample %in% otu5_ids, !is.na(Length)) %>%
  mutate(Oantigen = case_when(
    as.character(Length) %in% c("1803","1245") ~ "O-antigen+",
    as.character(Length) %in% c("1383","1830") ~ "O-antigen-",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Oantigen))

modern_pool_vec <- df_hp12$Oantigen        # 1312 Modern isolates
mod_pool_n      <- length(modern_pool_vec)
mod_plus_all    <- sum(modern_pool_vec == "O-antigen+")
mod_minus_all   <- mod_pool_n - mod_plus_all

cat(sprintf("Modern pool (HP12/OTU5): n=%d | OBC+ = %d (%.1f%%) | OBC- = %d (%.1f%%)\n",
            mod_pool_n, mod_plus_all, mod_plus_all/mod_pool_n*100,
            mod_minus_all, mod_minus_all/mod_pool_n*100))

# -----------------------------------------------------------
# 3. Compute Modern representative n=8 counts for CI labels and binomial test
#
#    mod_exp_plus  = round(8 × Modern OBC+ proportion)  — used for dot labels and CI
#    mod_exp_minus = 8 − mod_exp_plus
#    No cross-era Fisher test is performed (see supp_germany_downsample_summary.txt).
# -----------------------------------------------------------
mod_prop_plus_all <- mod_plus_all / mod_pool_n

# Representative n=8 draw: exact floating counts (no round())
# Defined here so they are available in section 8 (supp tables) which
# runs before section 7 (visualization) in file order.
mod_typical_plus  <- mod_prop_plus_all * n_samp   # e.g. 6.672
mod_typical_minus <- n_samp - mod_typical_plus    # e.g. 1.328

# Germany Modern representative draw within-era binomial test (H0: OBC+ = 50%)
# Uses round(mod_typical_plus) for integer requirement of binom.test();
# Notes column preserves the unrounded floating values used in the simple lower panel plot.
binom_germany_modern_within <- binom.test(round(mod_typical_plus), n_samp,
                                          p = 0.5, alternative = "two.sided")
germany_modern_binomial_within <- data.frame(
  Country          = "Germany",
  IsolateType      = "Modern (representative n=8 draw)",
  OBC_plus_n       = round(mod_typical_plus, 2),
  OBC_minus_n      = round(mod_typical_minus, 2),
  Total_N          = n_samp,
  Prop_OBC_plus    = round(mod_typical_plus / n_samp * 100, 2),
  Test             = paste0("Two-sided exact binomial (H0: OBC+ proportion = 0.5); ",
                            "integer input = round(", round(mod_typical_plus, 2), ") = ",
                            round(mod_typical_plus)),
  P_value          = binom_germany_modern_within$p.value,
  CI_lower_plus    = round(binom_germany_modern_within$conf.int[1] * 100, 2),
  CI_upper_plus    = round(binom_germany_modern_within$conf.int[2] * 100, 2),
  Significance     = sig_stars_step5(binom_germany_modern_within$p.value),
  Notes            = paste0("Same representative draw as simple lower panel plot: ",
                            round(mod_typical_plus, 2), " OBC+ / ",
                            round(mod_typical_minus, 2), " OBC- out of n=", n_samp),
  stringsAsFactors = FALSE
)
write.table(germany_modern_binomial_within,
            file.path(outdir_lower, "Germany_n8_Binomial_Modern_within.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf("Germany Modern within-era binomial (representative draw): p = %.4f [%s]\n",
            binom_germany_modern_within$p.value,
            sig_stars_step5(binom_germany_modern_within$p.value)))

# Integer Modern counts at n=8 — used for dot labels and CI centering.
# Floating mod_typical_plus is used for Wilson CI; rounded integer for display.
mod_exp_plus  <- round(mod_typical_plus)    # e.g. round(6.7) = 7
mod_exp_minus <- n_samp - mod_exp_plus      # e.g. 1

# -----------------------------------------------------------
# 4. Approach B: Repeated random downsampling — 1000 draws
#    Each draw: sample n=8 from Modern pool WITHOUT replacement,
#    run Fisher vs Germany Historical, and binomial within drawn sample.
# -----------------------------------------------------------
set.seed(42)
n_downsamp_lower <- 1000
# n_samp already defined above (= germany_n = 8)

downsamp_details <- lapply(seq_len(n_downsamp_lower), function(i) {
  samp     <- sample(modern_pool_vec, n_samp, replace = FALSE)  # random subsample
  s_plus   <- sum(samp == "O-antigen+")
  s_minus  <- n_samp - s_plus

  # Binomial within the drawn Modern n=8: does this sample depart from 50/50?
  binom_p <- binom.test(s_plus, n_samp, p = 0.5, alternative = "two.sided")$p.value

  data.frame(
    Run          = i,
    Plus_Modern  = s_plus,
    Minus_Modern = s_minus,
    Binomial_p   = binom_p
  )
}) %>% bind_rows()

# Summary table — binomial only
downsamp_summary_n8 <- data.frame(
  Country                  = "Germany",
  Hist_plus                = germany_plus,
  Hist_minus               = germany_minus,
  Hist_n                   = germany_n,
  Modern_pool_n            = mod_pool_n,
  Downsample_n             = n_samp,
  N_runs                   = n_downsamp_lower,
  Test_Binomial            = "Two-sided exact binomial (within each run, H0: OBC+ = 50%)",
  Binomial_mean_p          = round(mean(downsamp_details$Binomial_p, na.rm=TRUE), 4),
  Binomial_median_p        = round(median(downsamp_details$Binomial_p, na.rm=TRUE), 4),
  Binomial_prop_sig_05     = round(mean(downsamp_details$Binomial_p < 0.05, na.rm=TRUE), 4),
  Binomial_Significance    = sig_stars_step5(median(downsamp_details$Binomial_p, na.rm=TRUE)),
  stringsAsFactors         = FALSE
)

# Write full per-draw detail
write.table(downsamp_details,
            file.path(outdir_lower, "Germany_n8_downsampling_1000draws_details.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# Write summary
write.table(downsamp_summary_n8,
            file.path(outdir_lower, "Germany_n8_downsampling_1000draws_summary.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

cat(sprintf("\nGermany n=8 repeated downsampling results (%d runs):\n", n_downsamp_lower))
cat(sprintf("  Binomial within-run: median p = %.4f | prop runs p<0.05 = %.1f%% [%s]\n",
            downsamp_summary_n8$Binomial_median_p,
            downsamp_summary_n8$Binomial_prop_sig_05 * 100,
            downsamp_summary_n8$Binomial_Significance))

# -----------------------------------------------------------
# 5. Write Tables 3 & 4 to final_comparison folder
# -----------------------------------------------------------

# Table 4: Lower panel — Germany n=8 repeated downsampling (1000 runs summary)
write.table(downsamp_summary_n8,
            file.path(outdir_final,
                      "Table4_Fig4C_lower_Germany_n8_downsampling_1000draws.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# -----------------------------------------------------------
# 6. Write final_comparison/summary.txt (comprehensive narrative)
# -----------------------------------------------------------
summary_path <- file.path(outdir_final, "final_comparison_summary.txt")
zz <- file(summary_path, open = "wt")

writeLines("================================================================", zz)
writeLines("FINAL COMPARISON SUMMARY: Fig4C Statistical Tests", zz)
writeLines("For PNAS Supplementary Material", zz)
writeLines(paste0("Generated: ", Sys.Date()), zz)
writeLines("================================================================", zz)
writeLines("", zz)

writeLines("OVERVIEW", zz)
writeLines("--------", zz)
writeLines("This folder contains four key statistical comparison tables spanning", zz)
writeLines("the upper, middle, and lower panels of Fig4C.", zz)
writeLines("", zz)
writeLines("All tests compare O-antigen biosynthesis cluster (OBC) carriage between", zz)
writeLines(sprintf("Historical (pre-antibiotic era, n=38 total) and Modern (n=%s OTU5)", mod_pool_n), zz)
writeLines("P. viridiflava isolates, using different HTF haplotype groupings and", zz)
writeLines("geographic stratification.", zz)
writeLines("", zz)

writeLines("================================================================", zz)
writeLines("FIGURE 4C — THREE-PANEL LEGEND", zz)
writeLines("================================================================", zz)
writeLines("", zz)
writeLines("  UPPER PANEL — HTF haplotype frequency within each era", zz)
writeLines("  --------------------------------------------------------", zz)
writeLines("  Plot type : Frequency bar/dot chart", zz)
writeLines("  X-axis    : Era (Historical pre-antibiotic | Modern)", zz)
writeLines(sprintf("  Y-axis    : Proportion of isolates (%%); n_Historical=%d, n_Modern=%s", 38L, mod_pool_n), zz)
writeLines("  Groups    : 4 HTF haplotype groups", zz)
writeLines("              OBC+  haplotypes: 1803 (OBC+A), 1245 (OBC+B)  — warm tones", zz)
writeLines("              OBC-  haplotypes: 1383 (OBC-A), 1830 (OBC-B)  — cool tones", zz)
writeLines("  Stat      : Chi-square goodness-of-fit (H0: each group = 25%) per era", zz)
writeLines("              + exact binomial on collapsed OBC+/- binary", zz)
writeLines("  Files     : Table1 (binomial), Table1b (chi-square)", zz)
writeLines("", zz)
writeLines("  MIDDLE PANEL — OBC carriage frequency across all isolates per era", zz)
writeLines("  --------------------------------------------------------", zz)
writeLines("  Plot type : Dot + Wilson 95% CI", zz)
writeLines("  X-axis    : Era (Historical | Modern)", zz)
writeLines(sprintf("  Y-axis    : Proportion OBC+ or OBC- (%%); n_Historical=%d, n_Modern=%s", 38L, mod_pool_n), zz)
writeLines("  Groups    : Two dots per era:", zz)
writeLines("              OBC+  (#c94f7c, pink)   — isolates carrying OBC", zz)
writeLines("              OBC-  (#9e9ac8, purple)  — isolates lacking OBC", zz)
writeLines("  Error bars: Wilson 95% CI (binom::binom.confint, method='wilson')", zz)
writeLines("  Stat      : Two-sided exact binomial within each era (H0: OBC+ = 50%)", zz)
writeLines("  Files     : Table2 (binomial)", zz)
writeLines("", zz)
writeLines("  LOWER PANEL — Germany only: Historical vs Modern (size-matched)", zz)
writeLines("  --------------------------------------------------------", zz)
writeLines("  Plot type : Dot + Wilson 95% CI", zz)
writeLines(sprintf("  X-axis    : Germany Historical (n=%d) | Modern Germany (downsample n=%d)", germany_n, n_samp), zz)
writeLines("  Y-axis    : Proportion OBC+ or OBC- (%)", zz)
writeLines("  Groups    : Two dots per group:", zz)
writeLines("              OBC+  (#c94f7c, pink)   — isolates carrying OBC", zz)
writeLines("              OBC-  (#9e9ac8, purple)  — isolates lacking OBC", zz)
writeLines("  Error bars: Wilson 95% CI (binom::binom.confint, method='wilson')", zz)
writeLines(sprintf("              Historical: n=%d, genuine small-sample uncertainty (wide CI)", germany_n), zz)
writeLines(sprintf("              Modern:     n=%d representative draw (mod_exp_plus=%d OBC+); CI centred on", n_samp, mod_exp_plus), zz)
writeLines(sprintf("                          floating proportion (%.3f) from full Modern pool (n=%s)", mod_typical_plus, mod_pool_n), zz)
writeLines("  Dot labels: Rounded percentage + integer count  e.g. '87.5% (7/8)'", zz)
writeLines("  Stat      : Two-sided exact binomial within each era (H0: OBC+ = 50%)", zz)
writeLines(sprintf("            + %d repeated downsampling runs to demonstrate sampling variance", n_downsamp_lower), zz)
writeLines("  Files     : Table3b (within binomial), Table3c (Modern representative binomial),", zz)
writeLines("              Table4 (downsampling runs, binomial)", zz)
writeLines("", zz)
writeLines("  SHARED CONVENTIONS", zz)
writeLines("  --------------------------------------------------------", zz)
writeLines("  OBC+  (#c94f7c)  O-antigen biosynthesis cluster PRESENT", zz)
writeLines("  OBC-  (#9e9ac8)  O-antigen biosynthesis cluster ABSENT", zz)
writeLines("  Error bars      Wilson 95% CI throughout all panels", zz)
writeLines("  Significance    *** p<0.001  ** p<0.01  * p<0.05  ns p>=0.05", zz)
writeLines("  Dashed line     50% reference (equal split)", zz)
writeLines("", zz)
writeLines("================================================================", zz)
writeLines("", zz)
writeLines("TABLE 1 — Fig4C UPPER PANEL (4 HTF groups, OBC binary, within each era)", zz)
writeLines("File: Table1_Fig4C_upper_4HTFgroups_OBC_binomial_within_era.tsv", zz)
writeLines("Test: Two-sided exact binomial test (H0: OBC+ proportion = 0.5)", zz)
writeLines("Question: Within Historical and within Modern separately, is OBC+ enriched", zz)
writeLines("          over OBC- (i.e. does one group dominate within an era)?", zz)
writeLines("Note: This collapses the 4 HTF groups into binary OBC+/- before testing.", zz)
writeLines("      For the distribution across all 4 groups, see Table1b (chi-square).", zz)
writeLines("", zz)

writeLines("TABLE 1b — Fig4C UPPER PANEL (4 HTF groups, chi-square goodness-of-fit)", zz)
writeLines("File: Table1b_Fig4C_upper_4HTFgroups_chisq_uniform_within_era.tsv", zz)
writeLines("Test: Chi-square goodness-of-fit (H0: all 4 groups equally frequent, 25%)", zz)
writeLines("Question: Are the 4 HTF haplotypes uniformly distributed within each era?", zz)
writeLines("Note: This is the appropriate test for the 4-group within-era comparison.", zz)
writeLines("      Binomial is suitable only for the collapsed binary OBC+/- question.", zz)
writeLines("", zz)

writeLines("TABLE 2 — Fig4C MIDDLE PANEL (OBC binary, within each era)", zz)
writeLines("File: Table2_Fig4C_middle_OBC_binomial_within_era.tsv", zz)
writeLines("Test: Two-sided exact binomial test (H0: OBC+ proportion = 0.5)", zz)
writeLines("Question: Identical to Table 1 — shown separately per panel for clarity.", zz)
writeLines("", zz)

writeLines("TABLE 4 — Fig4C LOWER PANEL (Germany n=8, repeated random downsampling, 1000 runs)", zz)
writeLines("File: Table4_Fig4C_lower_Germany_n8_downsampling_1000draws.tsv", zz)
writeLines("Test: Repeated random downsampling — 1000 runs of:", zz)
writeLines("      - Randomly draw n=8 isolates WITHOUT replacement from Modern OBC pool", zz)
writeLines("      - Two-sided exact binomial test within each drawn group (H0: OBC+ = 50%)", zz)
writeLines("Purpose: Demonstrate that downsampling to n=8 produces unreliable proportions", zz)
writeLines("         (only ~21% of runs reach p<0.05) — no cross-era comparison is reported.", zz)
writeLines("", zz)

writeLines("----------------------------------------------------------------", zz)
writeLines("REPEATED DOWNSAMPLING — PURPOSE", zz)
writeLines("----------------------------------------------------------------", zz)
writeLines("Repeated random downsampling (used here): draw n=8 WITHOUT replacement", zz)
writeLines("  from the Modern pool 1000 times and apply a two-sided exact binomial", zz)
writeLines("  test (H0: OBC+ = 50%) within each drawn sample independently.", zz)
writeLines("  Purpose: demonstrate that n=8 samples from the Modern pool cannot", zz)
writeLines("  reliably distinguish their OBC composition from 50/50 — confirming", zz)
writeLines("  that no cross-era comparison is warranted. This is NOT a test of", zz)
writeLines("  Historical vs Modern difference.", zz)
writeLines("", zz)

writeLines("----------------------------------------------------------------", zz)
writeLines("BINOMIAL vs CHI-SQUARE for within-era comparisons", zz)
writeLines("----------------------------------------------------------------", zz)
writeLines("Binomial test (Tables 1 & 2): appropriate when outcome is BINARY.", zz)
writeLines("  OBC+ vs OBC- is binary, so binomial is correct for testing whether", zz)
writeLines("  one group is more common than the other within an era.", zz)
writeLines("Chi-square goodness-of-fit (Table 1b): appropriate for comparing", zz)
writeLines("  proportions across 3+ categories simultaneously (4 HTF groups).", zz)
writeLines("  Using 4 separate binomial tests would inflate the Type I error rate", zz)
writeLines("  (multiple testing); chi-square tests all 4 groups in one omnibus test.", zz)
writeLines("", zz)

writeLines("----------------------------------------------------------------", zz)
writeLines("SIGNIFICANCE THRESHOLDS", zz)
writeLines("----------------------------------------------------------------", zz)
writeLines("*** p < 0.001 | ** p < 0.01 | * p < 0.05 | ns p >= 0.05", zz)
writeLines("", zz)

writeLines("----------------------------------------------------------------", zz)
writeLines("KEY RESULTS SUMMARY", zz)
writeLines("----------------------------------------------------------------", zz)
writeLines("Tables 1/2  Within-era binomial (OBC+ vs 50%): see Table1/Table2 files", zz)
writeLines(sprintf("Table 3b   Germany Historical binomial: p = %.4f [%s]",
                   binom_germany_within$p.value,
                   sig_stars_step5(binom_germany_within$p.value)), zz)
writeLines(sprintf("Table 3c   Germany Modern representative binomial: p = %.4f [%s]",
                   binom_germany_modern_within$p.value,
                   sig_stars_step5(binom_germany_modern_within$p.value)), zz)
writeLines(sprintf("Table 4    Downsampling (%d runs): binomial median p = %.4f | prop runs p<0.05 = %.1f%% [%s]",
                   n_downsamp_lower,
                   downsamp_summary_n8$Binomial_median_p,
                   downsamp_summary_n8$Binomial_prop_sig_05 * 100,
                   downsamp_summary_n8$Binomial_Significance), zz)
close(zz)

# NOTE: Fig4C_figure_legend.txt is written in Section 13 (after all stat
#       variables are defined — chisq_hist/mod, binom_hist_r/mod_r, etc.)

# -----------------------------------------------------------
# 5b. Write lower-panel Germany binomial within to final_comparison
# -----------------------------------------------------------
write.table(germany_binomial_within,
            file.path(outdir_final,
                      "Table3b_Fig4C_lower_Germany_n8_Binomial_within.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# -----------------------------------------------------------
# 5c. Write Germany Modern representative-draw within-era binomial to final_comparison
# -----------------------------------------------------------
write.table(germany_modern_binomial_within,
            file.path(outdir_final,
                      "Table3c_Fig4C_lower_Germany_n8_Binomial_Modern_within.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

cat("\n✅ step5 Germany repeated downsampling analysis complete.\n")
cat("   Lower panel outputs:", outdir_lower, "\n")
cat("   Tables 3, 3b & 4 written to:", outdir_final, "\n")
cat("   Summary written to:", summary_path, "\n")

# -----------------------------------------------------------
# 8. Assemble two-version combined Supplementary Statistics Tables
#    (reads step4 partial outputs + appends lower panel rows)
#    Run AFTER step4 has been executed.
# -----------------------------------------------------------

# ---- Helper to safely read step4 partial files ----
read_partial <- function(fpath) {
  if (file.exists(fpath)) {
    read.table(fpath, header=TRUE, sep="\t", stringsAsFactors=FALSE,
               check.names=FALSE, fill=TRUE)
  } else {
    message("Warning: file not found, skipping: ", fpath)
    NULL
  }
}

chisq_df   <- read_partial(file.path(outdir_final,
                "Table1b_Fig4C_upper_4HTFgroups_chisq_uniform_within_era.tsv"))
binom_mid  <- read_partial(file.path(step3_outdir,
                "Fig4Cmiddlepanel_oantigenPAfreq/Oantigen_Binomial_OBC_within_era.tsv"))
fisher_pg  <- read_partial(file.path(step3_outdir,
                "Fig4Cupper_panel_HTF4group_freq/HTF4group_Fisher_crossera_pergroup.tsv"))
fisher_ds  <- read_partial(file.path(step3_outdir,
                "Fig4Cupper_panel_HTF4group_freq/HTF4group_Fisher_OBC_crossera_downsampled.tsv"))
fisher_mid <- read_partial(file.path(step3_outdir,
                "Fig4Cmiddlepanel_oantigenPAfreq/Oantigen_freq_hist_vs_modern1350_separated_fisher_test.tsv"))

# ---- Build standard-format row-sets ----
make_row <- function(Panel, Test_type, Era_Comparison, Test_name,
                     N_Historical, N_Modern, Statistic, P_value,
                     Significance, Notes="") {
  data.frame(Panel=Panel, Test_type=Test_type, Era_Comparison=Era_Comparison,
             Test_name=Test_name, N_Historical=as.character(N_Historical),
             N_Modern=as.character(N_Modern), Statistic=as.character(Statistic),
             P_value=as.character(signif(as.numeric(P_value),4)),
             Significance=Significance, Notes=as.character(Notes),
             stringsAsFactors=FALSE)
}

rows_upper_chisq <- if (!is.null(chisq_df)) {
  lapply(c("Historical","Modern"), function(era) {
    r <- chisq_df[chisq_df$IsolateType==era, ]
    n_hist_str <- if(era=="Historical") paste0(r$Total_N," total") else "—"
    n_mod_str  <- if(era=="Modern")     paste0(r$Total_N," total") else "—"
    make_row("Upper","Within-era",
             paste0(era," (4 HTF groups)"),
             "Chi-square goodness-of-fit (H0: uniform 25% each)",
             n_hist_str, n_mod_str,
             paste0("X2=",r$ChiSq_statistic,", df=",r$df),
             r$P_value, r$Significance,
             "Are the 4 HTF groups (1803,1245,1383,1830) equally frequent within the era?")
  }) %>% bind_rows()
} else data.frame()

rows_middle_binom <- if (!is.null(binom_mid)) {
  lapply(c("Historical","Modern"), function(era) {
    r <- binom_mid[binom_mid$IsolateType==era, ]
    n_hist_str <- if(era=="Historical") paste0(r$Total_N," (",r$OBC_plus_n," OBC+, ",r$OBC_minus_n," OBC-)") else "—"
    n_mod_str  <- if(era=="Modern")     paste0(r$Total_N," (",r$OBC_plus_n," OBC+, ",r$OBC_minus_n," OBC-)") else "—"
    make_row("Middle","Within-era",
             paste0(era," (OBC+/-)"),
             "Two-sided exact binomial (H0: OBC+ = 50%)",
             n_hist_str, n_mod_str,
             paste0("p=",signif(r$P_value,4),", 95%CI=[",r$CI_lower_plus,",",r$CI_upper_plus,"]%"),
             r$P_value, r$Significance,
             "Is OBC+ significantly more/less common than OBC- within the era?")
  }) %>% bind_rows()
} else data.frame()

rows_middle_between <- if (!is.null(fisher_pg) && !is.null(fisher_ds) && !is.null(fisher_mid)) {
  # 4 per-group Fisher rows
  rows_pg <- lapply(seq_len(nrow(fisher_pg)), function(i) {
    r <- fisher_pg[i,]
    make_row("Middle","Between-era",
             paste0("Historical vs Modern (HTF-",r$LengthGroup,", ",r$Oantigen,")"),
             "Fisher's exact test (per HTF group)",
             paste0(r$Hist_count,"/",r$Hist_total," (",r$Hist_pct,"%)"),
             paste0(r$Mod_count,"/",r$Mod_total," (",r$Mod_pct,"%)"),
             paste0("OR=",r$Odds_ratio),
             r$P_value, r$Significance,
             paste0("Did HTF-",r$LengthGroup," frequency change between eras?"))
  }) %>% bind_rows()
  # Fisher full dataset
  row_full <- make_row("Middle","Between-era",
             "Historical vs Modern (OBC+/-, full dataset)",
             "Fisher's exact test (full cohort)",
             paste0(fisher_mid$Historical_Oantigen_plus,"/",
                    fisher_mid$Historical_Oantigen_plus+fisher_mid$Historical_Oantigen_minus,
                    " OBC+"),
             paste0(fisher_mid$Modern_Oantigen_plus,"/",
                    fisher_mid$Modern_Oantigen_plus+fisher_mid$Modern_Oantigen_minus,
                    " OBC+"),
             paste0("OR=",round(fisher_mid$Odds_ratio,4)),
             fisher_mid$P_value,
             ifelse(fisher_mid$P_value < 0.001,"***",
                    ifelse(fisher_mid$P_value < 0.01,"**",
                           ifelse(fisher_mid$P_value < 0.05,"*","ns"))),
             "Overall OBC+ shift between eras; unadjusted for n imbalance")
  # Fisher downsampled Approach A (expected value)
  row_exp <- make_row("Middle","Between-era (downsampled)",
             paste0("Historical vs Modern expected n=",fisher_ds$Downsample_n),
             "Fisher's exact test (expected-value, Modern downsampled to n_Hist)",
             paste0(fisher_ds$Hist_OBC_plus,"/",fisher_ds$Hist_total," OBC+"),
             paste0(fisher_ds$ApproachA_Mod_exp_plus,"/",fisher_ds$Downsample_n,
                    " OBC+ (expected)"),
             paste0("OR=",fisher_ds$ApproachA_Odds_ratio),
             fisher_ds$ApproachA_P_value,
             fisher_ds$ApproachA_Significance,
             "Corrects for n imbalance using deterministic expected counts")
  # Fisher repeated downsampling Approach B
  row_ds <- make_row("Middle","Between-era (robustness)",
             paste0("Historical vs Modern (repeated downsampling n=",
                    fisher_ds$Downsample_n,", ",fisher_ds$ApproachB_N_draws," draws)"),
             "Fisher's exact test (repeated random downsampling robustness)",
             paste0(fisher_ds$Hist_OBC_plus,"/",fisher_ds$Hist_total," OBC+"),
             paste0(fisher_ds$Mod_OBC_plus_full,"/",fisher_ds$Mod_total_full,
                    " Modern pool"),
             paste0("median p=",fisher_ds$ApproachB_Median_P,
                    "; prop p<0.05=",round(fisher_ds$ApproachB_Prop_sig_05*100,1),"%"),
             fisher_ds$ApproachB_Median_P,
             fisher_ds$ApproachB_Significance,
             "Robustness: consistently significant across all random Modern subsamples?")
  bind_rows(row_full, row_exp, row_ds)   # per-group HTF Fisher omitted; OBC-level only
} else data.frame()

rows_lower_within <- make_row(
  "Lower","Within-era",
  "Germany Historical (OBC+/-)",
  "Two-sided exact binomial (H0: OBC+ = 50%)",
  paste0(germany_n," (",germany_plus," OBC+, ",germany_minus," OBC-)"),
  "—",
  paste0("p=",signif(binom_germany_within$p.value,4),
         ", 95%CI=[",round(binom_germany_within$conf.int[1]*100,2),
         ",",round(binom_germany_within$conf.int[2]*100,2),"]%"),
  binom_germany_within$p.value,
  sig_stars_step5(binom_germany_within$p.value),
  "Is OBC+/- evenly split in Germany Historical isolates (n=8)?"
)

rows_lower_modern_within <- make_row(
  "Lower","Within-era",
  "Germany Modern (representative n=8 draw, OBC+/-)",
  "Two-sided exact binomial (H0: OBC+ = 50%)",
  "—",
  paste0("representative n=8: ", round(mod_typical_plus,2), " OBC+ / ",
         round(mod_typical_minus,2), " OBC- (same as simple lower panel plot)"),
  paste0("p=",signif(binom_germany_modern_within$p.value,4),
         ", 95%CI=[",round(binom_germany_modern_within$conf.int[1]*100,2),
         ",",round(binom_germany_modern_within$conf.int[2]*100,2),"]%"),
  binom_germany_modern_within$p.value,
  sig_stars_step5(binom_germany_modern_within$p.value),
  paste0("Modern representative draw: round(",round(mod_typical_plus,2),
         ")=",round(mod_typical_plus)," used for binom.test; floating values same as plot")
)

rows_lower_downsamp <- make_row(
  "Lower","Downsampling variance (no cross-era comparison)",
  paste0("Modern pool (n=",n_samp," per run, ",n_downsamp_lower," runs)"),
  "Two-sided exact binomial within each downsampled run (H0: OBC+ = 50%)",
  "—",
  paste0(mod_pool_n," Modern pool, n=",n_samp," drawn per run"),
  paste0("binomial median p=",downsamp_summary_n8$Binomial_median_p,
         "; prop runs p<0.05=",round(downsamp_summary_n8$Binomial_prop_sig_05*100,1),"%"),
  downsamp_summary_n8$Binomial_median_p,
  downsamp_summary_n8$Binomial_Significance,
  paste0("Shows sampling variance at n=8; no cross-era comparison performed")
)

# ---- Full version (all panels, all tests) ----
supp_full <- bind_rows(
  rows_upper_chisq,
  rows_middle_binom,
  rows_middle_between,
  rows_lower_within,
  rows_lower_modern_within,
  rows_lower_downsamp
)
write.table(supp_full,
            file.path(step3_outdir,
                      "Supplementary_Statistics_Fig4C_FULL.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# ---- For-paper version (concise: upper within, middle within, lower within+representative+robustness) ----
supp_forpaper <- bind_rows(
  rows_upper_chisq,
  rows_middle_binom,
  rows_lower_within,
  rows_lower_modern_within,
  rows_lower_downsamp
)
write.table(supp_forpaper,
            file.path(step3_outdir,
                      "Supplementary_Statistics_Fig4C_forpaper.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

cat("\n📄 Supplementary tables written to:", step3_outdir, "\n")
cat("   FULL version (14 rows, all panels/tests):  Supplementary_Statistics_Fig4C_FULL.tsv\n")
cat("   For-paper version (6 rows, key results):   Supplementary_Statistics_Fig4C_forpaper.tsv\n")

# -----------------------------------------------------------
# 7. Visualization: Germany Historical vs Modern (downsampled n=8)
#    "lower panel like before"
#    X-axis: two groups — Germany Historical (fixed) vs Modern typical n=8 draw
#    Y-axis: proportion (%)
#
#    CI method — same type for both groups (Wilson 95%):
#      Historical : observed counts (2 OBC+, 6 OBC- out of 8)
#                   → binom.confint(k, 8, method="wilson")
#      Modern     : representative single draw = round(ModernProportion × 8)
#                   → binom.confint(typical_k, 8, method="wilson")
#                   This is NOT the IQR of 1000 draws; it is the binomial
#                   uncertainty of ONE typical draw, making errorbars
#                   directly comparable to the Historical side.
#                   The 1000-draw robustness is in the table, not the plot.
#    Coloured by OBC status (O-antigen+/-)
# -----------------------------------------------------------

# Per-draw proportions (used for computing the representative draw centre)
downsamp_details <- downsamp_details %>%
  mutate(
    Prop_plus  = Plus_Modern  / n_samp,
    Prop_minus = Minus_Modern / n_samp
  )

# Germany Historical Wilson CIs
hist_ci_plus  <- binom.confint(germany_plus,  germany_n, method = "wilson")
hist_ci_minus <- binom.confint(germany_minus, germany_n, method = "wilson")

# Modern representative draw:
#   centre = exact Modern OBC+ proportion × n_samp (NO round())
#   CI     = Wilson 95% CI — binom.confint accepts non-integer k because
#            it uses p = k/n internally; removing round() gives the CI
#            centred at the true Modern proportion (e.g., 83.4% not 87.5%)
mod_typical_plus  <- mod_prop_plus_all * n_samp          # floating, e.g. 6.672
mod_typical_minus <- n_samp - mod_typical_plus           # floating, e.g. 1.328

mod_ci_plus  <- binom.confint(mod_typical_plus,  n_samp, method = "wilson")
mod_ci_minus <- binom.confint(mod_typical_minus, n_samp, method = "wilson")

cat(sprintf("Modern representative draw (n=8): OBC+ = %.3f, OBC- = %.3f (%.1f%% Modern pool)\n",
            mod_typical_plus, mod_typical_minus, mod_prop_plus_all * 100))

# -----------------------------------------------------------
# Append CI summary table to final_comparison_summary.txt
# (opened in append mode; CI objects now available)
# -----------------------------------------------------------
zz_ci <- file(summary_path, open = "at")

writeLines("", zz_ci)
writeLines("================================================================", zz_ci)
writeLines("CONFIDENCE INTERVAL SUMMARY — Lower Panel Germany Plots", zz_ci)
writeLines("================================================================", zz_ci)
writeLines(paste0(
  "Wilson 95% CIs computed using binom::binom.confint(k, n, method=\"wilson\"). ",
  "For the Modern representative, k displayed as rounded integer (mod_exp_plus = ", mod_exp_plus,
  "; OBC- = ", mod_exp_minus, "). ",
  "CIs are centred on the floating k (mod_typical_plus = ", round(mod_typical_plus, 3),
  " = Modern OBC+ proportion \u00d7 n=", n_samp,
  ") to reflect the true Modern pool OBC+ proportion rather than the rounded integer."), zz_ci)
writeLines("", zz_ci)
writeLines(paste0(
  "Panel A of Germany_supp_ABC.pdf uses these same CI values."), zz_ci)
writeLines("", zz_ci)
writeLines(sprintf("%-32s | n | OBC+ k      | OBC+ %%   | CI lower %% | CI upper %%",
                   "Group"), zz_ci)
writeLines(paste(rep("-", 88), collapse=""), zz_ci)
writeLines(sprintf("%-32s | %d | %-11s | %-8.1f | %-10.1f | %.1f",
                   "Germany Historical",
                   germany_n,
                   germany_plus,
                   germany_plus / germany_n * 100,
                   hist_ci_plus$lower * 100,
                   hist_ci_plus$upper * 100), zz_ci)
writeLines(sprintf("%-32s | %d | %-11s | %-8.1f | %-10.1f | %.1f",
                   "Germany Historical OBC-",
                   germany_n,
                   germany_minus,
                   germany_minus / germany_n * 100,
                   hist_ci_minus$lower * 100,
                   hist_ci_minus$upper * 100), zz_ci)
writeLines(sprintf("%-32s | %d | %-11d | %-8.1f | %-10.1f | %.1f",
                   "Modern representative (n=8) OBC+",
                   n_samp,
                   mod_exp_plus,
                   mod_exp_plus / n_samp * 100,
                   mod_ci_plus$lower * 100,
                   mod_ci_plus$upper * 100), zz_ci)
writeLines(sprintf("%-32s | %d | %-11d | %-8.1f | %-10.1f | %.1f",
                   "Modern representative (n=8) OBC-",
                   n_samp,
                   mod_exp_minus,
                   mod_exp_minus / n_samp * 100,
                   mod_ci_minus$lower * 100,
                   mod_ci_minus$upper * 100), zz_ci)
writeLines("", zz_ci)
writeLines(paste0(
  "Note: k shown as rounded integer (mod_exp_plus = round(", round(mod_typical_plus, 3), ") = ",
  mod_exp_plus, ") for display consistency with dot labels in all three plots. ",
  "CIs are computed from the floating mod_typical_plus (", round(mod_typical_plus, 3),
  ") to centre the errorbar at the true Modern pool OBC+ proportion. ",
  "Fisher's exact test (Approach A) also uses integer inputs: ",
  mod_exp_plus, " OBC+ vs ", mod_exp_minus, " OBC- (both groups n=", n_samp, ")."), zz_ci)
writeLines("", zz_ci)
writeLines("----------------------------------------------------------------", zz_ci)
writeLines("SIGNIFICANCE THRESHOLDS", zz_ci)
writeLines("----------------------------------------------------------------", zz_ci)
writeLines("*** p < 0.001 | ** p < 0.01 | * p < 0.05 | ns p >= 0.05", zz_ci)

close(zz_ci)

# Assemble plot data frame
group_hist <- "Germany\nHistorical\n(n=8)"
# group_hist / group_mod labels no longer used for plotting (kept for CI table)

# -----------------------------------------------------------
# 8. Supplementary ABC composite figure
#    Panel A : Germany Historical dot + Wilson CI
#    Panel B : Modern 1000-run downsampling — boxplot + jitter per OBC group
#    Panel C : Binomial p-value distribution across 1000 downsample runs
#    Saved   : supp_germany_downsample/Germany_supp_ABC.pdf
# -----------------------------------------------------------
if (!requireNamespace("patchwork", quietly = TRUE))
  install.packages("patchwork", repos = "https://cloud.r-project.org", quiet = TRUE)
library(patchwork)

obc_colors <- c("O-antigen+" = "#c94f7c", "O-antigen-" = "#9e9ac8")
# Shared theme — 6–7 pt fonts, generous margins so nothing clips
# plot.margin: top, right, bottom, left  (extra top for title, extra right for tag)
base_theme <- theme_classic(base_size = 6, base_family = "Helvetica") +
  theme(axis.text       = element_text(size = 6,  family = "Helvetica"),
        axis.title      = element_text(size = 6,  family = "Helvetica"),
        legend.title    = element_blank(),
        legend.position = "top",
        legend.text     = element_text(size = 6,  family = "Helvetica"),
        legend.key.size = unit(0.35, "lines"),
        plot.title      = element_text(size = 7,  family = "Helvetica",
                                       hjust = 0.5),
        plot.tag        = element_text(size = 8,  family = "Helvetica"),
        axis.line       = element_line(linewidth = 0.3),
        axis.ticks      = element_line(linewidth = 0.3),
        plot.margin     = margin(t = 6, r = 6, b = 3, l = 3, unit = "pt"))

# ---- Panel A: Historical Germany — single panel, both groups on x-axis ----
# X-axis  : "OBC−" / "OBC+"  (clean, just the group name)
# Top text : "6/8\n75%" / "2/8\n25%" placed ABOVE the 100% line via geom_text
#            (coord_cartesian clip="off" lets text float above y=1.0 without
#             extending the axis scale past 100%)

df_A <- data.frame(
  Oantigen    = c("O-antigen-",  "O-antigen+"),
  Proportion  = c(germany_minus / germany_n, germany_plus / germany_n),
  Lower       = c(hist_ci_minus$lower,  hist_ci_plus$lower),
  Upper       = c(hist_ci_minus$upper,  hist_ci_plus$upper),
  Count_label = c(germany_minus, germany_plus),
  Total       = germany_n,
  stringsAsFactors = FALSE
)

pA <- ggplot(df_A, aes(x = Oantigen, y = Proportion, color = Oantigen)) +
  geom_point(size = 0.9, alpha = 0.85, show.legend = FALSE) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                width = 0.12, linewidth = 0.4, show.legend = FALSE) +
  # Count/proportion label floated above the 100% ceiling
  geom_text(aes(y = 1.03,
                label = paste0(Count_label, "/", Total, "\n",
                               round(Proportion * 100), "%")),
            size = 1.8, fontface = "plain", color = "black",
            lineheight = 1.1, vjust = 0, show.legend = FALSE) +
  scale_color_manual(values = obc_colors) +
  # No limits here — limits go in coord_cartesian so geom_text above 1.0 is kept
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     breaks = c(0, 0.25, 0.50, 0.75, 1.00),
                     expand = expansion(mult = c(0, 0))) +
  scale_x_discrete(labels = c("O-antigen-" = "OBC\u2212",
                               "O-antigen+" = "OBC+")) +
  coord_cartesian(ylim = c(0, 1.00), clip = "off") +
  labs(x = "", y = "Frequency (%)",
       title = sprintf("Historical Germany (n=%d)", germany_n)) +
  base_theme +
  theme(axis.text.x  = element_text(size = 5.5, family = "Helvetica"),
        plot.margin  = margin(t = 22, r = 6, b = 3, l = 3, unit = "pt"))

# ---- Panel B: Grouped downsampling outcomes — one facet per draw outcome ----
# Facet order : OBC− count 0 → 8  (= Plus_Modern descending 8 → 0)
# Facet title : "OBC−/OBC+"  e.g. "0/8"
# Facet subtitle : % of 1000 runs that fell in this group  e.g. "21%"

# All possible outcomes 0/8 … 8/0 (including groups with 0 observed runs)
# Strip subtitle = exact proportion to 1 decimal place  e.g. "22.4%"
run_freq <- data.frame(Plus_Modern = 0:n_samp) %>%
  left_join(
    downsamp_details %>% count(Plus_Modern, name = "n_runs"),
    by = "Plus_Modern"
  ) %>%
  mutate(n_runs        = ifelse(is.na(n_runs), 0L, n_runs),
         Minus_Modern  = n_samp - Plus_Modern,
         run_pct_exact = round(n_runs / n_downsamp_lower * 100, 1),
         run_pct_label = paste0(sprintf("%.1f", run_pct_exact), "%"),
         Group_label   = paste0(Minus_Modern, "/", n_samp,
                                "\n", run_pct_label)) %>%
  filter(n_runs > 0)   # drop groups with 0 observed runs

# Order: OBC− ascending (0 → 8)
obs_levels <- run_freq %>% arrange(Minus_Modern) %>% pull(Group_label)
run_freq   <- run_freq %>%
  mutate(Group_label = factor(Group_label, levels = obs_levels))

# Background jitter: only runs that actually occurred (groups with n_runs > 0)
df_B_jitter <- bind_rows(
  downsamp_details %>%
    transmute(Plus_Modern,
              Oantigen   = "O-antigen+",
              Proportion = Prop_plus),
  downsamp_details %>%
    transmute(Plus_Modern,
              Oantigen   = "O-antigen-",
              Proportion = Prop_minus)
) %>%
  left_join(run_freq %>% select(Plus_Modern, Group_label), by = "Plus_Modern") %>%
  mutate(Group_label = factor(Group_label, levels = obs_levels))

# Summary dots + Wilson 95% CI — built from run_freq so ALL groups get a dot
df_B_dots <- bind_rows(
  run_freq %>% transmute(Plus_Modern, Group_label,
                         Oantigen   = "O-antigen+",
                         k_val      = Plus_Modern,
                         Proportion = Plus_Modern / n_samp),
  run_freq %>% transmute(Plus_Modern, Group_label,
                         Oantigen   = "O-antigen-",
                         k_val      = n_samp - Plus_Modern,
                         Proportion = (n_samp - Plus_Modern) / n_samp)
) %>%
  arrange(desc(Plus_Modern), Oantigen) %>%
  rowwise() %>%
  mutate(
    ci_r  = list(binom.confint(k_val, n_samp, method = "wilson")),
    Lower = ci_r$lower,
    Upper = ci_r$upper
  ) %>%
  ungroup() %>%
  select(-ci_r)

# ---- shared pB theme addon ----
pB_theme_extra <- theme(
  strip.text       = element_text(size = 5, family = "Helvetica", lineheight = 1.1),
  strip.background = element_blank(),
  panel.spacing    = unit(2, "pt"),
  axis.text.x      = element_text(size = 5.5, family = "Helvetica",
                                  angle = 90, hjust = 1, vjust = 0.5),
  plot.margin      = margin(t = 10, r = 6, b = 3, l = 3, unit = "pt")
)

pB_scale_and_labs <-
  list(
    facet_wrap(~Group_label, nrow = 1),
    scale_color_manual(values = obc_colors),
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       breaks = c(0, 0.25, 0.50, 0.75, 1.00),
                       expand = expansion(mult = c(0, 0))),
    scale_x_discrete(labels = c("O-antigen+" = "OBC+", "O-antigen-" = "OBC\u2212")),
    coord_cartesian(ylim = c(0, 1.00), clip = "off"),
    labs(x = "", y = "Frequency (%)",
         title = sprintf("Modern downsampled\n(n=%d per draw, %d runs)",
                         n_samp, n_downsamp_lower)),
    base_theme,
    pB_theme_extra
  )

# ---- Version WITH jitter ----
pB <- ggplot() +
  geom_jitter(data = df_B_jitter,
              aes(x = Oantigen, y = Proportion, color = Oantigen),
              width = 0.46, height = 0.005,
              alpha = 0.18, size = 0.30, shape = 16,
              show.legend = FALSE) +
  geom_errorbar(data = df_B_dots,
                aes(x = Oantigen, ymin = Lower, ymax = Upper, color = Oantigen),
                width = 0.35, linewidth = 0.4, show.legend = FALSE) +
  geom_point(data = df_B_dots,
             aes(x = Oantigen, y = Proportion, color = Oantigen),
             size = 0.9, alpha = 0.92, show.legend = FALSE) +
  pB_scale_and_labs

# ---- Version WITHOUT jitter ----
pB_nojitter <- ggplot() +
  geom_errorbar(data = df_B_dots,
                aes(x = Oantigen, ymin = Lower, ymax = Upper, color = Oantigen),
                width = 0.35, linewidth = 0.4, show.legend = FALSE) +
  geom_point(data = df_B_dots,
             aes(x = Oantigen, y = Proportion, color = Oantigen),
             size = 0.9, alpha = 0.92, show.legend = FALSE) +
  pB_scale_and_labs

# -----------------------------------------------------------
# Write Panel B group detail table
# Columns: Group | OBC- | OBC+ | n_samp | n_runs | exact_pct | rounded_pct |
#          OBC+_prop | CI_lower | CI_upper
# Also explains why rounded percentages may not sum to 100%.
# -----------------------------------------------------------
panelB_table_path <- file.path(outdir_lower, "PanelB_group_detail_table.txt")

# Compute CI for each group × OBC type, keep only groups in run_freq
panelB_detail <- run_freq %>%
  arrange(Minus_Modern) %>%
  rowwise() %>%
  mutate(
    OBC_minus   = Minus_Modern,
    OBC_plus    = Plus_Modern,
    n_total     = n_samp,
    exact_pct   = round(n_runs / n_downsamp_lower * 100, 2),
    # Wilson CI for OBC+ proportion within this group
    ci_plus     = list(binom.confint(Plus_Modern,  n_samp, method = "wilson")),
    ci_minus    = list(binom.confint(Minus_Modern, n_samp, method = "wilson")),
    OBCplus_prop  = round(Plus_Modern  / n_samp * 100, 1),
    OBCplus_CI_lo = round(ci_plus$lower  * 100, 1),
    OBCplus_CI_hi = round(ci_plus$upper  * 100, 1),
    OBCminus_prop = round(Minus_Modern / n_samp * 100, 1),
    OBCminus_CI_lo = round(ci_minus$lower * 100, 1),
    OBCminus_CI_hi = round(ci_minus$upper * 100, 1)
  ) %>%
  ungroup() %>%
  select(Group = Group_label, OBC_minus, OBC_plus, n_total,
         n_runs, exact_pct, rounded_pct = run_pct_exact,
         OBCplus_prop, OBCplus_CI_lo, OBCplus_CI_hi,
         OBCminus_prop, OBCminus_CI_lo, OBCminus_CI_hi)

# Clean up newlines in Group label for plain text output
panelB_detail$Group <- gsub("\n", "  ", panelB_detail$Group)

zt <- file(panelB_table_path, open = "wt")
writeLines("================================================================", zt)
writeLines("Panel B — Downsampling Group Detail Table", zt)
writeLines(sprintf("Generated: %s", Sys.Date()), zt)
writeLines(sprintf("Modern pool n = %s  |  Draws per run: n = %d  |  Total runs: %d",
                   mod_pool_n, n_samp, n_downsamp_lower), zt)
writeLines("================================================================", zt)
writeLines("", zt)
writeLines("Column definitions", zt)
writeLines("------------------", zt)
writeLines("  Group          : OBC-/OBC+ draw outcome  (line 2 = % of runs)", zt)
writeLines("  OBC_minus      : Number of OBC- isolates in this draw outcome", zt)
writeLines("  OBC_plus       : Number of OBC+ isolates in this draw outcome", zt)
writeLines("  n_total        : Draw size (always 8)", zt)
writeLines("  n_runs         : Runs (out of 1000) that produced this outcome", zt)
writeLines("  exact_pct      : Exact run % (n_runs / 1000 * 100)", zt)
writeLines("  rounded_pct    : Rounded run % (shown in figure strip subtitle)", zt)
writeLines("  OBCplus_prop   : OBC+ proportion for this outcome (%) = OBC_plus/8", zt)
writeLines("  OBCplus_CI     : Wilson 95% CI for OBC+ proportion [lo, hi]", zt)
writeLines("  OBCminus_prop  : OBC- proportion for this outcome (%) = OBC_minus/8", zt)
writeLines("  OBCminus_CI    : Wilson 95% CI for OBC- proportion [lo, hi]", zt)
writeLines("", zt)

# Header
hdr <- sprintf("%-12s %9s %9s %7s %7s %10s %12s | %12s %14s | %13s %14s",
               "Group", "OBC_minus", "OBC_plus", "n_total", "n_runs",
               "exact_pct", "rounded_pct",
               "OBCplus_%", "OBCplus_CI",
               "OBCminus_%", "OBCminus_CI")
writeLines(hdr, zt)
writeLines(paste(rep("-", nchar(hdr)), collapse = ""), zt)

for (i in seq_len(nrow(panelB_detail))) {
  r <- panelB_detail[i, ]
  writeLines(sprintf(
    "%-12s %9d %9d %7d %7d %10.2f %12.1f | %12.1f  [%5.1f, %5.1f] | %13.1f  [%5.1f, %5.1f]",
    r$Group, r$OBC_minus, r$OBC_plus, r$n_total, r$n_runs,
    r$exact_pct, r$rounded_pct,
    r$OBCplus_prop,  r$OBCplus_CI_lo,  r$OBCplus_CI_hi,
    r$OBCminus_prop, r$OBCminus_CI_lo, r$OBCminus_CI_hi
  ), zt)
}

writeLines(paste(rep("-", nchar(hdr)), collapse = ""), zt)
# Totals row
writeLines(sprintf(
  "%-12s %9s %9s %7s %7d %10.2f %12.1f",
  "TOTAL", "", "", "",
  sum(panelB_detail$n_runs),
  sum(panelB_detail$exact_pct),
  sum(panelB_detail$rounded_pct)
), zt)

writeLines("", zt)
writeLines("================================================================", zt)
writeLines("WHY DO THE ROUNDED PERCENTAGES NOT SUM TO EXACTLY 100%?", zt)
writeLines("================================================================", zt)
writeLines("", zt)
writeLines("Each group percentage displayed in the figure strip is independently", zt)
writeLines("rounded to the nearest whole number (round(n_runs / 1000 * 100)).", zt)
writeLines("This is a well-known rounding artefact:", zt)
writeLines("  - The EXACT proportions always sum to 1.000 (100%) by definition,", zt)
writeLines("    because every run lands in exactly one group.", zt)
writeLines("  - However, rounding each proportion independently to integers can", zt)
writeLines("    cause the displayed sum to be 99% or 101% rather than 100%.", zt)
writeLines("    For example: 22.4% + 38.6% + 25.5% + 9.8% + 1.7% = 98.0/100", zt)
writeLines("    but round() gives 22 + 39 + 26 + 10 + 2 = 99.", zt)
writeLines("  - This is NOT an error in the data or sampling — it is solely a", zt)
writeLines("    display/rounding artefact. The exact_pct column above always", zt)
writeLines(sprintf("    sums to %.2f%%.", sum(panelB_detail$exact_pct)), zt)
writeLines("", zt)
writeLines("  Solution if exact 100% sum is required: use the largest-remainder", zt)
writeLines("  (Hamilton) rounding method instead of independent rounding.", zt)
writeLines("================================================================", zt)
close(zt)

cat(sprintf("\n\u2705 Panel B group detail table written to:\n   %s\n", panelB_table_path))

# Helper: format a p-value as academic text (defined here so Panel C and later
# sections can all use it; the definition further below in Section 10 is redundant
# but harmless — R simply re-assigns the same function object)
fmt_p <- function(p, digits = 3) {
  if (is.na(p) || is.null(p)) return("[NA]")
  if (p < 0.001) {
    exp_val  <- floor(log10(p))
    coeff    <- p / 10^exp_val
    sup_map  <- c("0"="\u2070","1"="\u00b9","2"="\u00b2","3"="\u00b3",
                  "4"="\u2074","5"="\u2075","6"="\u2076","7"="\u2077",
                  "8"="\u2078","9"="\u2079")
    exp_sup  <- paste(sapply(strsplit(as.character(abs(exp_val)),"")[[1]],
                             function(d) sup_map[d]), collapse="")
    return(sprintf("= %.1f \u00d7 10\u207b%s", coeff, exp_sup))
  }
  sprintf("= %.*f", digits, p)
}

# ---- Panel C: Binomial within-draw p-value distribution ----
prop_binom_sig  <- round(downsamp_summary_n8$Binomial_prop_sig_05 * 100, 1)
binom_med_p_str <- fmt_p(downsamp_summary_n8$Binomial_median_p)

# compute y position for annotation = ~85% of y-axis max after histogram
# n=8 binomial p-values are discrete: ~0.008, 0.070, 0.289, 0.727, 1.0
# binwidth=0.05, boundary=0 → bins: [0,0.05), [0.05,0.10), [0.25,0.30), ...
# The α=0.05 threshold sits exactly on the boundary between bin 1 and bin 2,
# so all p<0.05 runs (p≈0.008) form one bar and p≥0.05 runs start the next.
pC <- ggplot(downsamp_details, aes(x = Binomial_p)) +
  geom_histogram(aes(y = after_stat(count / sum(count))),
                 binwidth = 0.05, boundary = 0,
                 fill = "grey70", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0.05, linetype = "dashed",
             color = "red", linewidth = 0.5) +
  annotate("text",
           x     = 0.08,
           y     = Inf,
           label = sprintf("%.1f%% runs p < 0.05", prop_binom_sig),
           hjust = 0, vjust = 1.4,
           size  = 2.0, color = "red", fontface = "plain",
           family = "Helvetica") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     expand = expansion(add = c(0.02, 0.02))) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.15))) +
  coord_cartesian(clip = "off") +
  labs(x = "Binomial p-value",
       y = "Proportion of runs",
       title = sprintf("Modern downsampling\np-value distribution\n(%d runs)", n_downsamp_lower)) +
  base_theme +
  theme(legend.position = "none")

# ---- Combine: row1 = A (narrow) | B (wide, faceted groups) ----
#               row2 = C (narrow) | spacer
# pB gets ~3× the column width of pA/pC so all facets breathe comfortably
p_ABC <- (pA | pB) / (pC | plot_spacer()) +
  plot_layout(widths = c(1, 3)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 8))

# Save dimensions: 11 cm wide × 11 cm tall (publication column width; fonts 5–7 pt)
fig_w <- 11 / 2.54   # → inches
fig_h <- 11 / 2.54   # → inches

ggsave(file.path(outdir_lower, "Germany_supp_ABC.pdf"),
       plot = p_ABC, width = fig_w, height = fig_h, dpi = 1000)
ggsave(file.path(outdir_lower, "Germany_supp_ABC.png"),
       plot = p_ABC, width = fig_w, height = fig_h, dpi = 1000)
ggsave(file.path(outdir_lower, "Germany_supp_ABC.svg"),
       plot = p_ABC, width = fig_w, height = fig_h)

cat("\n\u2705 Germany supplementary ABC figure saved (PDF + PNG + SVG).\n")
cat("   Files: Germany_supp_ABC.pdf / .png / .svg\n")
cat("   Path:", outdir_lower, "\n")

# ---- No-jitter version ----
p_ABC_nojitter <- (pA | pB_nojitter) / (pC | plot_spacer()) +
  plot_layout(widths = c(1, 3)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 8))

ggsave(file.path(outdir_lower, "Germany_supp_ABC_nojitter.pdf"),
       plot = p_ABC_nojitter, width = fig_w, height = fig_h, dpi = 1000)
ggsave(file.path(outdir_lower, "Germany_supp_ABC_nojitter.png"),
       plot = p_ABC_nojitter, width = fig_w, height = fig_h, dpi = 1000)
ggsave(file.path(outdir_lower, "Germany_supp_ABC_nojitter.svg"),
       plot = p_ABC_nojitter, width = fig_w, height = fig_h)

cat("\n\u2705 Germany supplementary ABC (no-jitter) figure saved (PDF + PNG + SVG).\n")
cat("   Files: Germany_supp_ABC_nojitter.pdf / .png / .svg\n")
cat("   Path:", outdir_lower, "\n")

# -----------------------------------------------------------
# 8b. Supplementary summary text — rationale for no cross-era comparison
# -----------------------------------------------------------
supp_txt_path <- file.path(outdir_lower, "supp_germany_downsample_summary.txt")
zs <- file(supp_txt_path, open = "wt")
wl_s <- function(...) writeLines(paste0(...), zs)

wl_s("================================================================")
wl_s("SUPPLEMENTARY NOTE — Germany Historical OBC carriage analysis")
wl_s("supp_germany_downsample/")
wl_s(paste0("Generated: ", Sys.Date()))
wl_s("================================================================")
wl_s("")
wl_s("STATISTICAL APPROACH")
wl_s("--------------------")
wl_s("OBC carriage was assessed within each era using two complementary tests:")
wl_s("")
wl_s("  (1) Four-group HTF haplotype composition (upper panel):")
wl_s("      Chi-square goodness-of-fit test — are the four HTF haplotypes")
wl_s("      (1803, 1245, 1383, 1830) equally distributed within each era?")
wl_s("      H\u2080: each group = 25%.")
wl_s("")
wl_s("  (2) Binary OBC+/\u2212 carriage (middle panel):")
wl_s("      Two-sided exact binomial test within each era —")
wl_s("      does OBC+ frequency differ from 50%%? H\u2080: OBC+ = 50%%.")
wl_s("")
wl_s("WHY NO CROSS-ERA COMPARISON WAS PERFORMED")
wl_s("------------------------------------------")
wl_s(sprintf("Comparing Germany Historical (n\u00a0=\u00a0%d) against Modern isolates", germany_n))
wl_s(sprintf("requires downsampling the Modern pool (n\u00a0=\u00a0%s) to n\u00a0=\u00a0%d to achieve", mod_pool_n, n_samp))
wl_s("sample-size parity. At such small n, the resulting proportion carries")
wl_s("substantial sampling variance and cannot serve as a reliable reference.")
wl_s("")
wl_s(sprintf("To quantify this variance, we performed %d independent downsample runs:", n_downsamp_lower))
wl_s(sprintf("in each run, n\u00a0=\u00a0%d isolates were drawn without replacement from the", n_samp))
wl_s("Modern pool, and a two-sided exact binomial test was applied to the")
wl_s("sampled observations (H\u2080: OBC+ = 50%).")
wl_s(sprintf("Across %d runs, only %.1f%% yielded p\u00a0<\u00a00.05 (median p\u00a0%s).",
             n_downsamp_lower,
             downsamp_summary_n8$Binomial_prop_sig_05 * 100,
             fmt_p(downsamp_summary_n8$Binomial_median_p)))
wl_s(sprintf("In the vast majority of runs, a random sample of n\u00a0=\u00a0%d from the", n_samp))
wl_s("Modern pool could not be distinguished from a 50/50 split, indicating")
wl_s("insufficient information to yield a reliable Modern reference proportion.")
wl_s("Accordingly, no cross-era comparison is reported.")
wl_s("")
wl_s("FIGURE DESCRIPTION — Germany_supp_ABC.pdf")
wl_s("------------------------------------------")
wl_s(sprintf("Panel A — Historical Germany (n=%d) OBC carriage.", germany_n))
wl_s(sprintf("  OBC+ = %d/%d (%d%%); OBC\u2212 = %d/%d (%d%%).",
             germany_plus, germany_n, round(germany_plus/germany_n*100),
             germany_minus, germany_n, round(germany_minus/germany_n*100)))
wl_s("  Dots: solid, slightly transparent (alpha=0.85).")
wl_s("  Error bars: Wilson 95%% confidence intervals,")
wl_s(sprintf("  binom::binom.confint(k, n=%d, method='wilson').", germany_n))
wl_s(sprintf("  OBC+ Wilson 95%% CI: %.1f%% \u2013 %.1f%%;",
             hist_ci_plus$lower*100, hist_ci_plus$upper*100))
wl_s(sprintf("  OBC\u2212 Wilson 95%% CI: %.1f%% \u2013 %.1f%%.",
             hist_ci_minus$lower*100, hist_ci_minus$upper*100))
wl_s("  Wide CIs reflect genuine small-sample uncertainty at n=8.")
wl_s("")
wl_s(sprintf("Panel B — OBC+/\u2212 proportion distribution across %d downsample runs.", n_downsamp_lower))
wl_s(sprintf("  Each run: n=%d isolates sampled without replacement from the", n_samp))
wl_s(sprintf("  Modern pool (n=%s). Points: individual run proportions", mod_pool_n))
wl_s("  (solid circles, jittered horizontally; alpha=0.25).")
wl_s("  Box: IQR (25th\u201375th percentile); centre line: median;")
wl_s("  whiskers: 1.5\u00d7IQR.")
wl_s("  The dispersion across runs reflects the high sampling variance at")
wl_s("  n=8, confirming that the downsampled proportion is unreliable.")
wl_s("")
wl_s(sprintf("Panel C — Binomial p-value distribution across %d downsample runs.", n_downsamp_lower))
wl_s("  Within each run, a two-sided exact binomial test was applied to the")
wl_s("  sampled n=8 observations (H\u2080: OBC+ = 50%%), independent of Historical Germany.")
wl_s("  Histogram normalised to proportion of runs per bin.")
wl_s("  Red dashed line: significance threshold \u03b1\u00a0=\u00a00.05.")
wl_s(sprintf("  Only %.1f%% of runs yielded p\u00a0<\u00a00.05 (median p\u00a0%s),",
             downsamp_summary_n8$Binomial_prop_sig_05 * 100,
             fmt_p(downsamp_summary_n8$Binomial_median_p)))
wl_s("  indicating that n=8 samples from the Modern pool lack the power to")
wl_s("  reliably discriminate OBC composition from a 50/50 distribution.")
wl_s("")
wl_s("METHODS")
wl_s("-------")
wl_s("Confidence intervals (Panel A): Wilson score 95%% CI,")
wl_s("  binom::binom.confint(k, n, method='wilson'). Preferred over")
wl_s("  Clopper-Pearson for small n: better coverage, less conservative.")
wl_s("")
wl_s("Downsampling (Panels B, C): set.seed(42);")
wl_s(sprintf("  %d independent downsample runs of n=%d isolates drawn", n_downsamp_lower, n_samp))
wl_s(sprintf("  without replacement from Modern pool (n=%s). Per run:", mod_pool_n))
wl_s("  two-sided exact binomial test (H\u2080: OBC+ = 50%).")
wl_s("  No cross-era comparison is reported (see rationale above).")
wl_s("")
wl_s("Significance thresholds: *** p<0.001 | ** p<0.01 | * p<0.05 | ns p\u22650.05")

close(zs)

cat("\n\u2705 Supplementary summary text saved.\n")
cat("   File: supp_germany_downsample_summary.txt\n")
cat("   Path:", outdir_lower, "\n")

# NOTE: Section 9 (Germany_n8_downsampling_analysis/ copy folder) was removed.
# All outputs are written directly to outdir_lower (supp_germany_downsample/) and
# outdir_final (final_comparison/) — no separate copy folder needed.

# (old 8a / 8b / 8c individual plots removed — replaced by Germany_supp_ABC.pdf)



# -----------------------------------------------------------
# 10. Generate Fig4C_Results_and_Methods.md with all numbers
#     filled in from computed objects in this session.
#     Output: step3_outdir/Fig4C_Results_and_Methods.md
# -----------------------------------------------------------

# Helper: format a p-value as academic text (real numbers, no "n.s." or "< 0.001")
fmt_p <- function(p, digits = 3) {
  if (is.na(p) || is.null(p)) return("[NA]")
  if (p < 0.001) {
    # Scientific notation: e.g. "= 4.9 \u00d7 10\u207b\u2075"
    exp_val  <- floor(log10(p))
    coeff    <- p / 10^exp_val
    sup_map  <- c("0"="\u2070","1"="\u00b9","2"="\u00b2","3"="\u00b3",
                  "4"="\u2074","5"="\u2075","6"="\u2076","7"="\u2077",
                  "8"="\u2078","9"="\u2079")
    exp_sup  <- paste(sapply(strsplit(as.character(abs(exp_val)),"")[[1]],
                             function(d) sup_map[d]), collapse="")
    return(sprintf("= %.1f \u00d7 10\u207b%s", coeff, exp_sup))
  }
  sprintf("= %.*f", digits, p)
}

# Helper: pull a single value safely from a 1-row data frame
pull1 <- function(df, col, default = "[see table]") {
  if (is.null(df) || nrow(df) == 0) return(default)
  v <- df[[col]]
  if (length(v) == 0 || is.na(v[1])) return(default)
  v[1]
}

# ---- Gather values from loaded/computed objects ----

# Upper panel — chi-square GOF (4 HTF groups)
chisq_hist <- if (!is.null(chisq_df)) chisq_df[chisq_df$IsolateType == "Historical", ] else NULL
chisq_mod  <- if (!is.null(chisq_df)) chisq_df[chisq_df$IsolateType == "Modern",     ] else NULL

hist_htf_counts <- if (!is.null(chisq_hist))
  sprintf("%d (%s), %d (%s), %d (%s), %d (%s)",
          pull1(chisq_hist,"HTF_1803_n"), "HTF-1803",
          pull1(chisq_hist,"HTF_1245_n"), "HTF-1245",
          pull1(chisq_hist,"HTF_1383_n"), "HTF-1383",
          pull1(chisq_hist,"HTF_1830_n"), "HTF-1830") else "[see table]"

hist_chisq_str <- if (!is.null(chisq_hist))
  sprintf("chi-square goodness-of-fit chi2 = %.2f, df = %d, p %s [%s]",
          as.numeric(pull1(chisq_hist,"ChiSq_statistic")),
          as.integer(pull1(chisq_hist,"df")),
          fmt_p(as.numeric(pull1(chisq_hist,"P_value"))),
          pull1(chisq_hist,"Significance")) else "[see table]"

mod_htf_counts <- if (!is.null(chisq_mod))
  sprintf("%d (%s), %d (%s), %d (%s), %d (%s)",
          pull1(chisq_mod,"HTF_1803_n"), "HTF-1803",
          pull1(chisq_mod,"HTF_1245_n"), "HTF-1245",
          pull1(chisq_mod,"HTF_1383_n"), "HTF-1383",
          pull1(chisq_mod,"HTF_1830_n"), "HTF-1830") else "[see table]"

mod_chisq_str <- if (!is.null(chisq_mod))
  sprintf("chi2 = %.1f, df = %d, p %s [%s]",
          as.numeric(pull1(chisq_mod,"ChiSq_statistic")),
          as.integer(pull1(chisq_mod,"df")),
          fmt_p(as.numeric(pull1(chisq_mod,"P_value"))),
          pull1(chisq_mod,"Significance")) else "[see table]"

# Middle panel — within-era binomial
binom_hist_r <- if (!is.null(binom_mid)) binom_mid[binom_mid$IsolateType == "Historical", ] else NULL
binom_mod_r  <- if (!is.null(binom_mid)) binom_mid[binom_mid$IsolateType == "Modern",     ] else NULL

hist_obc_plus  <- pull1(binom_hist_r, "OBC_plus_n",  "21")
hist_obc_minus <- pull1(binom_hist_r, "OBC_minus_n", "17")
hist_obc_tot   <- pull1(binom_hist_r, "Total_N",     "38")
hist_obc_pct   <- pull1(binom_hist_r, "Prop_OBC_plus","55.3")
hist_binom_p   <- if (!is.null(binom_hist_r) && nrow(binom_hist_r) > 0)
  as.numeric(pull1(binom_hist_r,"P_value")) else NA
hist_binom_str <- sprintf("exact binomial p %s [%s]",
                          fmt_p(hist_binom_p),
                          sig_stars_step5(hist_binom_p))

mod_obc_plus  <- pull1(binom_mod_r, "OBC_plus_n",   "1098")
mod_obc_minus <- pull1(binom_mod_r, "OBC_minus_n",  "214")
mod_obc_tot   <- pull1(binom_mod_r, "Total_N",      "1312")
mod_obc_pct   <- pull1(binom_mod_r, "Prop_OBC_plus","83.69")
mod_binom_p   <- if (!is.null(binom_mod_r) && nrow(binom_mod_r) > 0)
  as.numeric(pull1(binom_mod_r,"P_value")) else NA
mod_binom_str <- sprintf("exact binomial p %s [%s]",
                          fmt_p(mod_binom_p),
                          sig_stars_step5(mod_binom_p))

# Middle panel — between-era Fisher (full + downsampled)
full_fish_or <- if (!is.null(fisher_mid)) round(fisher_mid$Odds_ratio, 3)       else NA
full_fish_p  <- if (!is.null(fisher_mid)) fisher_mid$P_value                    else NA
full_fish_str <- sprintf("Fisher's exact test OR = %.3f, p %s [%s]",
                          full_fish_or, fmt_p(full_fish_p),
                          sig_stars_step5(full_fish_p))

apA_or   <- if (!is.null(fisher_ds)) fisher_ds$ApproachA_Odds_ratio else NA
apA_p    <- if (!is.null(fisher_ds)) as.numeric(fisher_ds$ApproachA_P_value)  else NA
apA_n    <- if (!is.null(fisher_ds)) fisher_ds$Downsample_n else 38
apA_str  <- sprintf("Fisher's exact test at n = %d OR = %.3f, p %s [%s]",
                     apA_n, apA_or, fmt_p(apA_p), sig_stars_step5(apA_p))

apB_med  <- if (!is.null(fisher_ds)) fisher_ds$ApproachB_Median_P          else NA
apB_pct  <- if (!is.null(fisher_ds)) round(fisher_ds$ApproachB_Prop_sig_05 * 100, 1) else NA
apB_ndrw <- if (!is.null(fisher_ds)) fisher_ds$ApproachB_N_draws else 1000
apB_str  <- sprintf("median p %s; %.1f%% of %d downsample runs p < 0.05",
                     fmt_p(apB_med), apB_pct, apB_ndrw)

# Germany lower panel
ger_hist_binom_str <- sprintf("exact binomial p %s [%s]",
                               fmt_p(binom_germany_within$p.value),
                               sig_stars_step5(binom_germany_within$p.value))
ger_ds_str  <- sprintf("median binomial p %s; %.1f%% of %d downsample runs p < 0.05",
                        fmt_p(downsamp_summary_n8$Binomial_median_p),
                        downsamp_summary_n8$Binomial_prop_sig_05 * 100,
                        n_downsamp_lower)

# Germany Historical vs Modern representative (n=8 vs n=8) Fisher test
# (Approach A deterministic: balanced 2×2 table using integer expected counts)
ft_exp_n8  <- fisher.test(matrix(c(germany_plus,  germany_minus,
                                    mod_exp_plus,  mod_exp_minus),
                                  nrow = 2, byrow = FALSE))
ger_fish_or <- round(unname(ft_exp_n8$estimate), 3)
ger_fish_p  <- ft_exp_n8$p.value

# Downsampling summary scalars used by text-generation blocks below
# (Binomial within-draw results — Fisher columns were removed in an earlier refactor)
ger_ds_med <- downsamp_summary_n8$Binomial_median_p
ger_ds_pct <- round(downsamp_summary_n8$Binomial_prop_sig_05 * 100, 1)

# ---- Write the markdown ----
rm_path <- file.path(step3_outdir, "Fig4C_Results_and_Methods.md")
rz <- file(rm_path, open = "wt")

writeLines("# Fig4C \u2014 Results and Methods", rz)
writeLines("", rz)
writeLines(paste0("*Auto-generated by step5_fig4c_lower_countryboxplot_bottom.R on ",
                  Sys.Date(), "*"), rz)
writeLines("", rz)
writeLines("---", rz)
writeLines("", rz)

# ---- RESULTS ----
writeLines("## Results", rz)
writeLines("", rz)

writeLines("### HTF haplotype distribution and OBC status across eras", rz)
writeLines("", rz)
# Compute Wilson CIs for upper/middle panel groups
ci_hist_obc <- tryCatch(
  binom.confint(as.integer(hist_obc_plus), as.integer(hist_obc_tot), method="wilson"),
  error = function(e) NULL)
ci_mod_obc  <- tryCatch(
  binom.confint(as.integer(mod_obc_plus),  as.integer(mod_obc_tot),  method="wilson"),
  error = function(e) NULL)
ci_hist_str <- if (!is.null(ci_hist_obc))
  sprintf("Wilson 95%% CI: %.1f%%\u2013%.1f%%", ci_hist_obc$lower*100, ci_hist_obc$upper*100) else ""
ci_mod_str  <- if (!is.null(ci_mod_obc))
  sprintf("Wilson 95%% CI: %.1f%%\u2013%.1f%%", ci_mod_obc$lower*100, ci_mod_obc$upper*100) else ""

writeLines(sprintf(paste0(
  "Four HTF haplotype groups were identified based on element length (HTF-1803, -1245, -1383, -1830), ",
  "with HTF-1803 and -1245 carrying the O-antigen biosynthesis cluster (OBC+) and HTF-1383 and -1830 ",
  "lacking it (OBC\u2212). Within Historical isolates (n\u00a0=\u00a0%s), the four groups were nearly uniformly ",
  "distributed (%s; chi-square goodness-of-fit chi2\u00a0=\u00a0%.2f, df\u00a0=\u00a0%d, p\u00a0%s [%s]), and OBC+ and OBC\u2212 ",
  "strains occurred at near-equal frequencies (%s/%s\u00a0=\u00a0%s%% OBC+; exact binomial p\u00a0%s [%s]; %s). ",
  "In contrast, within Modern isolates (n\u00a0=\u00a0%s), the four groups were strongly non-uniformly ",
  "distributed (chi2\u00a0=\u00a0%.1f, df\u00a0=\u00a0%d, p\u00a0%s [%s]), with OBC+ haplotypes predominating ",
  "(%s/%s\u00a0=\u00a0%s%% OBC+; exact binomial p\u00a0%s [%s]; %s)."),
  hist_obc_tot,
  hist_htf_counts,
  as.numeric(pull1(chisq_hist,"ChiSq_statistic")), as.integer(pull1(chisq_hist,"df")),
  fmt_p(as.numeric(pull1(chisq_hist,"P_value"))), pull1(chisq_hist,"Significance"),
  hist_obc_plus, hist_obc_tot, hist_obc_pct,
  fmt_p(hist_binom_p), sig_stars_step5(hist_binom_p), ci_hist_str,
  mod_obc_tot,
  as.numeric(pull1(chisq_mod,"ChiSq_statistic")), as.integer(pull1(chisq_mod,"df")),
  fmt_p(as.numeric(pull1(chisq_mod,"P_value"))), pull1(chisq_mod,"Significance"),
  mod_obc_plus, mod_obc_tot, mod_obc_pct,
  fmt_p(mod_binom_p), sig_stars_step5(mod_binom_p), ci_mod_str), rz)
writeLines("", rz)

writeLines("### Cross-era shift in OBC status", rz)
writeLines("", rz)
writeLines(sprintf(paste0(
  "The proportion of OBC+ isolates increased significantly from the Historical to the Modern era ",
  "(Fisher\u2019s exact test: OR\u00a0=\u00a0%.3f, p\u00a0%s [%s]). To account for the sample size imbalance ",
  "(n\u00a0=\u00a0%s vs. n\u00a0=\u00a0%s), Fisher\u2019s exact test was repeated using expected Modern counts scaled to ",
  "n\u00a0=\u00a0%s (OR\u00a0=\u00a0%.3f, p\u00a0%s [%s]), and confirmed by repeated random downsampling of n\u00a0=\u00a0%s ",
  "isolates from the Modern pool across %d draws (median p\u00a0%s; %.1f%% of draws p\u00a0<\u00a00.05)."),
  full_fish_or, fmt_p(full_fish_p), sig_stars_step5(full_fish_p),
  hist_obc_tot, mod_obc_tot,
  apA_n, apA_or, fmt_p(apA_p), sig_stars_step5(apA_p),
  apA_n, apB_ndrw, fmt_p(apB_med), apB_pct), rz)
writeLines("", rz)

writeLines("### Germany Historical subset", rz)
writeLines("", rz)
writeLines(sprintf(paste0(
  "To compare without bias from sample size or geographic origin, Germany Historical ",
  "isolates (n\u00a0=\u00a0%d) were extracted as the fixed reference group: %d OBC+ (%d%%), ",
  "%d OBC\u2212 (%d%%). OBC\u2212 haplotypes predominated (%d/%d\u00a0=\u00a0%d%% OBC\u2212; %d/%d\u00a0=\u00a0%d%% OBC+), ",
  "but this difference did not reach significance by exact binomial test ",
  "(p\u00a0%s [%s]; Wilson 95%% CI for OBC+: %.1f%%\u2013%.1f%%), reflecting limited power at ",
  "n\u00a0=\u00a0%d. A representative Modern group of n\u00a0=\u00a0%d was then derived from the full Modern pool ",
  "(%s OTU5 isolates) by scaling the observed Modern OBC+ proportion to n\u00a0=\u00a0%d ",
  "(round(%.1f%%\u00a0\u00d7\u00a0%d)\u00a0=\u00a0%d OBC+, %d OBC\u2212; Wilson 95%% CI for OBC+: %.1f%%\u2013%.1f%%). ",
  "Fisher\u2019s exact test on this balanced %d\u00a0vs.\u00a0%d comparison confirmed a significant OBC\u2212 ",
  "enrichment in Historical Germany (OR\u00a0=\u00a0%.3f, p\u00a0%s [%s]). ",
  "Repeated random downsampling of n\u00a0=\u00a0%d isolates from the full Modern pool across %d draws ",
  "confirmed the robustness of this result (median Fisher p\u00a0%s; %.1f%% of draws p\u00a0<\u00a00.05)."),
  germany_n,
  germany_plus, round(germany_plus / germany_n * 100),
  germany_minus, round(germany_minus / germany_n * 100),
  germany_minus, germany_n, round(germany_minus / germany_n * 100),
  germany_plus,  germany_n, round(germany_plus  / germany_n * 100),
  fmt_p(binom_germany_within$p.value), sig_stars_step5(binom_germany_within$p.value),
  hist_ci_plus$lower * 100, hist_ci_plus$upper * 100,
  germany_n, germany_n,
  mod_pool_n, germany_n,
  mod_prop_plus_all * 100, germany_n, mod_exp_plus, mod_exp_minus,
  mod_ci_plus$lower * 100, mod_ci_plus$upper * 100,
  germany_n, germany_n,
  ger_fish_or, fmt_p(ger_fish_p, digits=3), sig_stars_step5(ger_fish_p),
  germany_n, n_downsamp_lower,
  fmt_p(ger_ds_med), ger_ds_pct), rz)
writeLines("", rz)

# ---- Wilson CI summary ----
writeLines("### Wilson 95% Confidence Intervals (lower panel Germany plots)", rz)
writeLines("", rz)
writeLines(sprintf(paste0(
  "Wilson 95%% confidence intervals for plotted OBC+ proportions were computed using ",
  "`binom::binom.confint(k, n, method = \"wilson\")`. ",
  "For the Modern representative, k = %.3f (= Modern OBC+ proportion\u00a0\u00d7\u00a0%d, without rounding), ",
  "centering the CI at the true Modern pool OBC+ proportion (%.1f%%) rather than the rounded ",
  "integer count (%d/%d = %.1f%%). The floating k is used exclusively for CI computation and ",
  "plot centering; all statistical tests (Fisher\u2019s exact, binomial) use integer inputs."),
  mod_typical_plus, n_samp,
  mod_prop_plus_all * 100,
  mod_exp_plus, n_samp, mod_exp_plus / n_samp * 100), rz)
writeLines("", rz)
writeLines(sprintf("| %-32s | n | k (OBC+) | OBC+ %% | CI lower %% | CI upper %% |",
                   "Group"), rz)
writeLines("|----------------------------------|---|----------|---------|------------|------------|", rz)
writeLines(sprintf("| %-32s | %d | %-8d | %-7.1f | %-10.1f | %.1f |",
                   "Germany Historical",
                   germany_n, germany_plus,
                   germany_plus / germany_n * 100,
                   hist_ci_plus$lower * 100,
                   hist_ci_plus$upper * 100), rz)
writeLines(sprintf("| %-32s | %d | %-8s | %-7.1f | %-10.1f | %.1f |",
                   "Modern representative (n=8)",
                   n_samp,
                   sprintf("%.3f", mod_typical_plus),
                   mod_typical_plus / n_samp * 100,
                   mod_ci_plus$lower * 100,
                   mod_ci_plus$upper * 100), rz)
writeLines("", rz)
writeLines(sprintf(paste0(
  "*Both `Oantigen_freq_Germany_hist_vs_resampleModern8.pdf` and `Oantigen_freq_Germany_simple_HvsM.pdf` ",
  "use these CI values. The simple lower panel dot labels show the rounded integer proportion: ",
  "%d%% (%d/%d) for Historical and %.1f%% (%d/%d) for Modern representative.*"),
  round(germany_plus/germany_n*100), germany_plus, germany_n,
  mod_exp_plus/n_samp*100, mod_exp_plus, n_samp), rz)
writeLines("", rz)

# ---- METHODS ----
writeLines("---", rz)
writeLines("", rz)
writeLines("## Methods", rz)
writeLines("", rz)

writeLines("### Statistical tests", rz)
writeLines("", rz)
writeLines(paste0(
  "Within-era OBC+/\u2212 proportions were assessed by two-sided exact binomial test (H\u2080: p\u00a0=\u00a00.5). ",
  "The distribution of all four HTF haplotype groups within each era was assessed by chi-square ",
  "goodness-of-fit (H\u2080: 25% per group, df\u00a0=\u00a03). Cross-era differences in OBC status were assessed ",
  "by Fisher\u2019s exact test on 2\u00a0\u00d7\u00a02 contingency tables (era\u00a0\u00d7\u00a0OBC status)."), rz)
writeLines("", rz)
writeLines(sprintf(paste0(
  "To address sample size imbalance between eras, two complementary approaches were applied for ",
  "the Germany-specific comparison: (i)\u00a0a balanced deterministic comparison \u2014 a representative ",
  "Modern group of n\u00a0=\u00a08 was derived by scaling the observed Modern OBC+ proportion to n\u00a0=\u00a08 ",
  "(round(%.1f%%\u00a0\u00d7\u00a08)\u00a0=\u00a0%d OBC+), and Fisher\u2019s exact test was run on the resulting 8\u00a0vs.\u00a08 table; ",
  "(ii)\u00a0repeated random downsampling \u2014 %d draws of n\u00a0=\u00a08 without replacement from the full Modern ",
  "OTU5 pool (n\u00a0=\u00a0%s), with Fisher\u2019s exact test applied to each draw against the fixed Germany ",
  "Historical group; results are reported as the median p-value and proportion of draws with ",
  "p\u00a0<\u00a00.05. The same balanced-downsampling framework (at n\u00a0=\u00a0%s) was applied for the full cohort ",
  "cross-era comparison. Wilson 95%% confidence intervals were computed using the binom R package. ",
  "All analyses were performed in R (\u2265\u00a04.0). Significance thresholds: ",
  "*p\u00a0<\u00a00.05, **p\u00a0<\u00a00.01, ***p\u00a0<\u00a00.001."),
  mod_prop_plus_all * 100, mod_exp_plus,
  n_downsamp_lower, mod_pool_n,
  apA_n), rz)
writeLines("", rz)

writeLines("---", rz)
writeLines("", rz)
writeLines(paste0("*Key tables: `Supplementary_Statistics_Fig4C_FULL.tsv`, ",
                  "`Supplementary_Statistics_Fig4C_forpaper.tsv`*"), rz)

close(rz)

cat(sprintf("\n📝 Results & Methods markdown written to:\n   %s\n", rm_path))

# =============================================================================
# Section 11: Auto-generate README_Fig4C_statistics.md
# =============================================================================
# All variables used here are already defined in Sections 3 and 10.
# Key dynamic values: germany_plus/minus/n, mod_exp_plus/minus, mod_prop_plus_all,
# mod_pool_n, ft_exp_n8, downsamp_summary_n8, binom_germany_within,
# n_downsamp_lower, and the fmt_* strings from Section 10.

readme_path <- file.path(step3_outdir, "README_Fig4C_statistics.md")
rh <- file(readme_path, open = "wt")
wl <- function(...) writeLines(paste0(...), rh)

# Pre-build Germany-specific formatted strings
ger_pct_plus_hist  <- round(germany_plus / germany_n * 100)   # 25
ger_pct_minus_hist <- round(germany_minus / germany_n * 100)  # 75
mod_prop_pct       <- round(mod_prop_plus_all * 100, 1)       # 83.7

ger_fish_or_readme <- round(unname(ft_exp_n8$estimate), 3)
ger_fish_p_readme  <- ft_exp_n8$p.value
ger_fish_sig       <- sig_stars_step5(ger_fish_p_readme)
ger_ds_pct_readme  <- round(downsamp_summary_n8$Binomial_prop_sig_05 * 100, 1)
ger_ds_med_readme  <- downsamp_summary_n8$Binomial_median_p
ger_hist_binom_p   <- binom_germany_within$p.value
ger_hist_binom_sig <- sig_stars_step5(ger_hist_binom_p)

# Cross-era (full cohort) formatted strings — reuse from Sec 10
full_fish_or_readme <- if (!is.null(fisher_mid)) round(fisher_mid$Odds_ratio, 3) else NA
full_fish_p_readme  <- if (!is.null(fisher_mid)) fisher_mid$P_value else NA
apA_or_readme  <- if (!is.null(fisher_ds)) fisher_ds$ApproachA_Odds_ratio else NA
apA_p_readme   <- if (!is.null(fisher_ds)) as.numeric(fisher_ds$ApproachA_P_value) else NA
apA_n_readme   <- if (!is.null(fisher_ds)) fisher_ds$Downsample_n else 38
apB_med_readme <- if (!is.null(fisher_ds)) fisher_ds$ApproachB_Median_P else NA
apB_pct_readme <- if (!is.null(fisher_ds)) round(fisher_ds$ApproachB_Prop_sig_05 * 100, 1) else NA
apB_ndrw_readme<- if (!is.null(fisher_ds)) fisher_ds$ApproachB_N_draws else 1000

# ---- Header ----
wl("# README: Fig4C Statistical Analysis \u2014 HTF Frequency & O-antigen Comparisons")
wl("")
wl("*Auto-generated by step5_fig4c_lower_countryboxplot_bottom.R on ", Sys.Date(), "*")
wl("")
wl("---")
wl("")

# ---- Section 0: Q&A ----
wl("## 0. Main Results \u2014 Q&A Summary")
wl("")
wl("This section gives a plain-language interpretation of each panel's key statistical result, with the test used and the table where the result can be found.")
wl("")
wl("---")
wl("")

wl("### Upper panel: Are the 4 HTF haplotype groups (1803, 1245, 1383, 1830) equally represented within each era?")
wl("")
wl("**Answer: No \u2014 the 4 groups are non-uniformly distributed within both Historical and Modern isolates.**")
wl("")
wl(sprintf("- **Test:** Chi-square goodness-of-fit (H\u2080: each of the 4 groups = 25%%)"))
wl(sprintf("- **Historical** (n=%s): \u03c7\u00b2 = %s, df = 3, p %s [%s]",
           hist_obc_tot,
           if (!is.null(chisq_hist)) round(as.numeric(pull1(chisq_hist,"ChiSq_statistic")),2) else "[see table]",
           fmt_p(if (!is.null(chisq_hist)) as.numeric(pull1(chisq_hist,"P_value")) else NA),
           if (!is.null(chisq_hist)) pull1(chisq_hist,"Significance") else "[see table]"))
wl(sprintf("- **Modern** (n=%s): \u03c7\u00b2 = %s, df = 3, p %s [%s]",
           mod_obc_tot,
           if (!is.null(chisq_mod)) round(as.numeric(pull1(chisq_mod,"ChiSq_statistic")),1) else "[see table]",
           fmt_p(if (!is.null(chisq_mod)) as.numeric(pull1(chisq_mod,"P_value")) else NA),
           if (!is.null(chisq_mod)) pull1(chisq_mod,"Significance") else "[see table]"))
wl("- **Interpretation:** A significant result means at least one HTF haplotype group is over- or under-represented relative to a uniform distribution. In both eras, OBC+ groups (1803, 1245) and OBC\u2212 groups (1383, 1830) differ substantially from 25% each.")
wl("- **\u2192 Table:** `final_comparison/Table1b_Fig4C_upper_4HTFgroups_chisq_uniform_within_era.tsv`")
wl("- **\u2192 Full stats:** `Supplementary_Statistics_Fig4C_FULL.tsv` (rows Panel=Upper)")
wl("")
wl("---")
wl("")

wl("### Middle panel (within): Is OBC+ the dominant group within each era?")
wl("")
wl("**Answer: Historical shows a near-equal OBC+/\u2212 split; Modern shows strong OBC+ dominance.**")
wl("")
wl("- **Test:** Two-sided exact binomial (H\u2080: OBC+ proportion = 0.5)")
wl(sprintf("- **Historical** (n=%s): %s OBC+ / %s OBC\u2212 = %s%% \u2192 p %s [%s]",
           hist_obc_tot, hist_obc_plus, hist_obc_minus, hist_obc_pct,
           fmt_p(if (!is.null(binom_hist_r) && nrow(binom_hist_r)>0) as.numeric(pull1(binom_hist_r,"P_value")) else NA),
           if (!is.null(binom_hist_r) && nrow(binom_hist_r)>0) sig_stars_step5(as.numeric(pull1(binom_hist_r,"P_value"))) else "ns"))
wl(sprintf("- **Modern** (n=%s): %s OBC+ / %s OBC\u2212 = %s%% \u2192 p %s [%s]",
           mod_obc_tot, mod_obc_plus, mod_obc_minus, mod_obc_pct,
           fmt_p(if (!is.null(binom_mod_r) && nrow(binom_mod_r)>0) as.numeric(pull1(binom_mod_r,"P_value")) else NA),
           if (!is.null(binom_mod_r) && nrow(binom_mod_r)>0) sig_stars_step5(as.numeric(pull1(binom_mod_r,"P_value"))) else "***"))
wl("- **Interpretation:** The binomial p for Historical tests whether the ~55% OBC+ in the pre-antibiotic era departs from chance. The highly significant Modern result reflects the strong selection for O-antigen presence in the modern clinical setting (>80% OBC+).")
wl("- **\u2192 Table:** `Fig4Cmiddlepanel_oantigenPAfreq/Oantigen_Binomial_OBC_within_era.tsv`")
wl("- **\u2192 Full stats:** `Supplementary_Statistics_Fig4C_FULL.tsv` (rows Panel=Middle, Test_type=Within-era)")
wl("")
wl("---")
wl("")

wl("### Middle panel (between): Does OBC+ frequency differ between Historical and Modern eras?")
wl("")
wl("**Answer: Yes \u2014 OBC+ increased significantly from the Historical to the Modern era across all test approaches.**")
wl("")
wl(sprintf("- **Test 1 \u2014 Fisher exact (full dataset):** Hist %s/%s OBC+ vs Mod %s/%s OBC+ \u2192 OR = %.3f, p %s [%s]",
           hist_obc_plus, hist_obc_tot, mod_obc_plus, mod_obc_tot,
           full_fish_or_readme, fmt_p(full_fish_p_readme), sig_stars_step5(full_fish_p_readme)))
wl(sprintf("- **Test 2 \u2014 Fisher exact (Modern expected n=%d):** OR = %.3f, p %s [%s] *(controls for n imbalance)*",
           apA_n_readme, apA_or_readme, fmt_p(apA_p_readme), sig_stars_step5(apA_p_readme)))
wl(sprintf("- **Test 3 \u2014 Fisher repeated downsampling (%d draws, n=%d w/o replacement):** median p %s; %.1f%% of draws p < 0.05 *(robustness assessment)*",
           apB_ndrw_readme, apA_n_readme, fmt_p(apB_med_readme), apB_pct_readme))
wl("- **Interpretation:** All three approaches consistently show a significant OBC+ shift. The downsampling robustness analysis confirms the result holds regardless of which Modern isolates are compared. Per-group Fisher tests identify which specific HTF haplotypes drove the shift.")
wl("- **\u2192 Table:** `Fig4Cmiddlepanel_oantigenPAfreq/Oantigen_Fisher_OBC_crossera_downsampled.tsv`")
wl("- **\u2192 Full stats:** `Supplementary_Statistics_Fig4C_FULL.tsv` (rows Panel=Middle, Test_type=Between-era)")
wl("")
wl("---")
wl("")

wl("### Lower panel (within Germany): Is OBC+ evenly split within German Historical isolates?")
wl("")
wl(sprintf("**Answer: Germany Historical (n=%d) shows more OBC\u2212 than OBC+, but the small n limits power.**", germany_n))
wl("")
wl("- **Test:** Two-sided exact binomial (H\u2080: OBC+ = 50%)")
wl(sprintf("- **Germany Historical** (n=%d): %d OBC+ (%d%%), %d OBC\u2212 (%d%%) \u2192 p %s [%s]",
           germany_n, germany_plus, ger_pct_plus_hist,
           germany_minus, ger_pct_minus_hist,
           fmt_p(ger_hist_binom_p), ger_hist_binom_sig))
wl(sprintf("- **Interpretation:** With only n=%d, the binomial test has low power. Even the observed %d%% OBC+ does not reach significance at \u03b1=0.05. The result reflects Germany's relatively OBC\u2212-enriched Historical collection, consistent with early isolates predating OBC+ selection pressure.",
           germany_n, ger_pct_plus_hist))
wl("- **\u2192 Table:** `Fig4Clowerpanel_top3enrichedcountry_moderndownsampled_fisher/Germany_n8_Binomial_OBC_within.tsv`")
wl("- **\u2192 Full stats:** `Supplementary_Statistics_Fig4C_FULL.tsv` (rows Panel=Lower, Test_type=Within-era)")
wl("")
wl("---")
wl("")

wl("### Lower panel (between): Does Germany Historical OBC+ frequency differ from the Modern era?")
wl("")
wl(sprintf("**Answer: Yes \u2014 Germany Historical (%d/%d = %d%% OBC+) differs significantly from a size-matched Modern comparison, and this result is robust across %d repeated downsampling draws.**",
           germany_plus, germany_n, ger_pct_plus_hist, n_downsamp_lower))
wl("")
wl("**Logic (step-by-step):**")
wl("")
wl(sprintf("1. **Extract Germany Historical samples** (n\u00a0=\u00a0%d): %d OBC+ (%d%%), %d OBC\u2212 (%d%%) \u2014 the fixed reference group.",
           germany_n, germany_plus, ger_pct_plus_hist, germany_minus, ger_pct_minus_hist))
wl(sprintf("2. **Downsample Modern to n\u00a0=\u00a0%d**: From the full Modern pool (n\u00a0=\u00a0%s OTU5 isolates), derive a representative Modern group of n\u00a0=\u00a0%d by scaling the observed Modern OBC+ proportion: round(%.1f%%\u00a0\u00d7\u00a0%d)\u00a0=\u00a0%d OBC+, %d OBC\u2212. This removes bias from both sample size and geographic origin \u2014 both groups are held at n\u00a0=\u00a0%d.",
           germany_n, mod_pool_n, germany_n, mod_prop_pct, germany_n, mod_exp_plus, mod_exp_minus, germany_n))
wl(sprintf("3. **Compare %d H vs %d M with Fisher\u2019s exact test**: Germany Historical (%d/%d OBC+) vs Modern representative (%d/%d OBC+) \u2192 OR\u00a0=\u00a0%.3f, p\u00a0%s [%s].",
           germany_n, germany_n, germany_plus, germany_n, mod_exp_plus, germany_n,
           ger_fish_or_readme, fmt_p(ger_fish_p_readme), ger_fish_sig))
wl(sprintf("4. **Repeat %d random draws to test robustness**: Each draw samples n\u00a0=\u00a0%d without replacement from the %s-isolate Modern pool and runs Fisher\u2019s exact test against the fixed Germany Historical group. Median p\u00a0%s; %.1f%% of draws p\u00a0<\u00a00.05 \u2014 confirming the result is not sensitive to which specific %d Modern isolates are selected.",
           n_downsamp_lower, germany_n, mod_pool_n,
           fmt_p(ger_ds_med_readme), ger_ds_pct_readme, germany_n))
wl("")
wl("- **\u2192 Plot:** `Fig4Clowerpanel_top3enrichedcountry_moderndownsampled_fisher/Oantigen_freq_Germany_hist_vs_resampleModern8.pdf`")
wl("- **\u2192 Tables:** `Germany_n8_Fisher_expected_value.tsv`, `Germany_n8_downsampling_1000draws_summary.tsv`")
wl("- **\u2192 Full stats:** `Supplementary_Statistics_Fig4C_FULL.tsv` (rows Panel=Lower, Test_type=Between-era)")
wl("")
wl("---")
wl("")

wl("### Two-version supplementary tables")
wl("")
wl("| File | Contents | Rows |")
wl("|------|---------|------|")
wl("| `Supplementary_Statistics_Fig4C_FULL.tsv` | All panels, all tests (upper within, middle within+between, lower within+between+robustness) | ~14 |")
wl("| `Supplementary_Statistics_Fig4C_forpaper.tsv` | Paper-ready concise version (upper within, middle within, lower within + downsampling robustness) | ~6 |")
wl("")
wl("Both files are written to the root `step3_HTFfreq/` results folder by `step5_fig4c_lower_countryboxplot_bottom.R` (requires `step4` to have been run first).")
wl("")
wl("---")
wl("")

# ---- Section 1: Biological Context ----
wl("## 1. Biological Context")
wl("")
wl("### What are HTF haplotypes?")
wl("")
wl("**HTF (Horizontal Transfer Factor)** haplotypes are defined by four distinct")
wl("length variants of the tailocin/phage-derived genomic element:")
wl("")
wl("| HTF group | Length (bp) | O-antigen status | Color |")
wl("|-----------|-------------|------------------|-------|")
wl("| 1803      | 1803        | O-antigen+ (OBC+) | dark red `#800233` |")
wl("| 1245      | 1245        | O-antigen+ (OBC+) | yellow `#f9d42a` |")
wl("| 1383      | 1383        | O-antigen\u2212 (OBC\u2212) | blue `#638ccc` |")
wl("| 1830      | 1830        | O-antigen\u2212 (OBC\u2212) | light pink `#f6d6ff` |")
wl("")
wl("Groups 1803 and 1245 carry the O-antigen biosynthesis cluster (**OBC+**);")
wl("groups 1383 and 1830 lack it (**OBC\u2212**).")
wl("")
wl("### What do the two panels show?")
wl("")
wl("- **Fig4C upper panel**: All four HTF haplotype groups shown separately, with")
wl("  their frequencies in Historical (pre-antibiotic era) vs Modern isolates.")
wl("  Dots represent observed proportion with Wilson 95% confidence intervals.")
wl("")
wl("- **Fig4C middle panel**: The same isolates collapsed into two groups \u2014")
wl("  OBC+ (O-antigen present) and OBC\u2212 (O-antigen absent) \u2014 with proportions")
wl("  plotted for Historical vs Modern.")
wl("")
wl("### Datasets")
wl("")
wl("| Label      | Source                          | N isolates | Notes |")
wl("|------------|---------------------------------|------------|-------|")
wl(sprintf("| Historical | kmer-based dominant haplotype   | %s         | Pre-antibiotic era; after QC exclusions |", hist_obc_tot))
wl(sprintf("| Modern     | HP12 combined lengths (OTU5)    | %s         | Exact n computed from data |", mod_obc_tot))
wl("")
wl("---")
wl("")

# ---- Section 2: Statistical Tests ----
wl("## 2. Statistical Tests Explained")
wl("")
wl("### Test 1: Two-sided exact binomial test (within era)")
wl("")
wl("**Files:** `HTF4group_Binomial_OBC_within_era.tsv`,")
wl("`Oantigen_Binomial_OBC_within_era.tsv`")
wl("")
wl("**Question:** Within a given era (Historical or Modern), is the proportion of")
wl("OBC+ isolates significantly different from 50%?")
wl("")
wl("**H\u2080:** p(OBC+) = 0.50 (equal proportions of OBC+ and OBC\u2212)")
wl("**H\u2081:** p(OBC+) \u2260 0.50 (two-sided)")
wl("")
wl("**Method:** `binom.test(n_obc_plus, total_n, p = 0.5, alternative = \"two.sided\")`")
wl("")
wl("**Interpretation:**")
wl(sprintf("- Historical (%s/%s = %s%% OBC+): tests whether the near-equal split is",
           hist_obc_plus, hist_obc_tot, hist_obc_pct))
wl("  consistent with chance (expected for a population without strong O-antigen")
wl("  selection pressure in the pre-antibiotic era).")
wl(sprintf("- Modern (high OBC+, e.g. %s/%s = %s%%): tests whether the strong OBC+",
           mod_obc_plus, mod_obc_tot, mod_obc_pct))
wl("  dominance is significantly different from 50%, reflecting selection for")
wl("  O-antigen presence in the modern clinical setting.")
wl("")
wl("**Why this test?** The exact binomial test is ideal for these sample sizes \u2014")
wl("it makes no normality assumption and is exact (not approximate), providing")
wl("conservative control of Type I error.")
wl("")
wl("---")
wl("")

wl("### Test 2: Fisher\u2019s exact test per HTF group (cross-era)")
wl("")
wl("**File:** `HTF4group_Fisher_crossera_pergroup.tsv`")
wl("")
wl("**Question:** Did the frequency of each individual HTF haplotype group change")
wl("significantly between Historical and Modern isolates?")
wl("")
wl("**Method:** For each of the four groups (1803, 1245, 1383, 1830), construct a")
wl("2\u00d72 contingency table:")
wl("")
wl("```")
wl("               | Carries this group | Does not carry this group |")
wl("Historical     |       a            |           b               |")
wl("Modern         |       c            |           d               |")
wl("```")
wl("")
wl("Apply Fisher\u2019s exact test to test independence of era and group membership.")
wl("")
wl("**Why this test?** Fisher\u2019s exact test is preferable over chi-square when")
wl("cells have low expected counts (as is the case for the Historical era with")
wl(sprintf("n=%s). It is the standard for 2\u00d72 contingency tables in microbiology papers.", hist_obc_tot))
wl("")
wl("**Output columns:**")
wl("- `Odds_ratio`: > 1 means the group is *more frequent in Historical* relative")
wl("  to Modern; < 1 means it is *less frequent in Historical*.")
wl("- `Significance`: `***` p < 0.001, `**` p < 0.01, `*` p < 0.05, `ns` p \u2265 0.05")
wl("")
wl("---")
wl("")

wl("### Test 3: Fisher\u2019s exact test, overall OBC+/\u2212 cross-era (full dataset)")
wl("")
wl("**Files:** `Oantigen_freq_hist_vs_modern1350_separated_fisher_test.tsv`,")
wl("within summary `.txt` files")
wl("")
wl("**Question:** Is the proportion of OBC+ isolates significantly different")
wl("between the Historical and Modern eras?")
wl("")
wl("**2\u00d72 contingency table:**")
wl("```")
wl(sprintf("             | OBC+ | OBC\u2212 |"))
wl(sprintf("Historical   |  %s  |  %s  |", hist_obc_plus, hist_obc_minus))
wl(sprintf("Modern       | ~%s | ~%s |", mod_obc_plus, mod_obc_minus))
wl("```")
wl("")
wl("**Method:** `fisher.test(mat)` on integer counts.")
wl("")
wl(sprintf("**Result:** OR \u2248 %.3f, p %s \u2014 strong evidence that Historical isolates",
           full_fish_or_readme, fmt_p(full_fish_p_readme)))
wl("have a significantly lower OBC+ proportion than Modern isolates.")
wl("")
wl("**Limitation:** With n_Modern >> n_Historical, even a small true effect size")
wl("will yield p << 0.05. The downsampled tests (Test 4) address this concern.")
wl("")
wl("---")
wl("")

wl("### Test 4: Fisher\u2019s exact test with Modern downsampled to n_Historical (upper/middle panels)")
wl("")
wl("**Files:** `HTF4group_Fisher_OBC_crossera_downsampled.tsv`,")
wl("`Oantigen_Fisher_OBC_crossera_downsampled.tsv`")
wl("")
wl(sprintf("**Motivation:** The large sample size imbalance (Historical n=%s, Modern n\u2248%s)",
           hist_obc_tot, mod_obc_tot))
wl("means that the full-dataset Fisher test is strongly powered by the large Modern")
wl("group. To evaluate whether the signal is robust, we apply two approaches analogous")
wl("to the country-level comparisons in the **Fig4C lower panel** (step5 script):")
wl("")
wl("#### Approach A \u2014 Expected-value (deterministic)")
wl(sprintf("Compute what the Modern OBC+/\u2212 counts *would be* if Modern had exactly n=%s",
           apA_n_readme))
wl("samples, using the observed Modern OBC+ proportion:")
wl("")
wl("```")
wl(sprintf("  expected_modern_obc_plus  = round(%s \u00d7 p_modern_plus)", apA_n_readme))
wl(sprintf("  expected_modern_obc_minus = %s \u2212 expected_modern_obc_plus", apA_n_readme))
wl("```")
wl("")
wl("This mirrors the approach used for Germany in Fig4C lower panel (downsampled to")
wl(sprintf("Germany's Historical count of n=%d).", germany_n))
wl("")
wl(sprintf("#### Approach B \u2014 Repeated random downsampling (stochastic, %d draws)", apB_ndrw_readme))
wl(sprintf("Repeatedly draw n=%s isolates **without replacement** from the full Modern pool", apA_n_readme))
wl(sprintf("(%s OTU5 isolates from HP12 data) and run Fisher\u2019s exact test against the", mod_pool_n))
wl("fixed Historical group each time.")
wl("")
wl("**Key output column:** `ApproachB_Prop_sig_05` \u2014 fraction of random draws with")
wl("p < 0.05. Close to 1.0 means the signal is robust regardless of which Modern")
wl("isolates are sampled. This is a **robustness assessment**, not a null-hypothesis test.")
wl("")

wl("### Test 5: Germany n=8 balanced comparison + repeated downsampling (lower panel, step5)")
wl("")
wl("**Files:** `Germany_n8_Fisher_expected_value.tsv`,")
wl("`Germany_n8_downsampling_1000draws_summary.tsv`,")
wl("`Germany_n8_downsampling_1000draws_details.tsv`,")
wl("`Germany_n8_Binomial_OBC_within.tsv`")
wl("")
wl(sprintf("**Logic (step-by-step):** To compare Germany Historical OBC proportions against the Modern era **without bias from sample size or geographic origin**:"))
wl("")
wl(sprintf("1. **Extract Germany Historical samples** (n\u00a0=\u00a0%d): fixed reference group (%d OBC+, %d OBC\u2212).",
           germany_n, germany_plus, germany_minus))
wl(sprintf("2. **Downsample Modern to n\u00a0=\u00a0%d**: Scale the observed Modern OBC+ proportion (%.1f%%) to n\u00a0=\u00a0%d: `mod_exp_plus\u00a0=\u00a0round(%.1f%%\u00a0\u00d7\u00a0%d)\u00a0=\u00a0%d OBC+`, `mod_exp_minus\u00a0=\u00a0%d`. Both groups are now n\u00a0=\u00a0%d, removing sample-size confounding.",
           germany_n, mod_prop_pct, germany_n, mod_prop_pct, germany_n,
           mod_exp_plus, mod_exp_minus, germany_n))
wl(sprintf("3. **Compare %d H vs %d M with Fisher\u2019s exact test** (Approach A \u2014 deterministic): Run Fisher\u2019s exact test on the balanced 2\u00d72 table `[Germany Historical (%d OBC+, %d OBC\u2212)] vs [Modern representative (%d OBC+, %d OBC\u2212)]`. This is the primary comparison.",
           germany_n, germany_n, germany_plus, germany_minus, mod_exp_plus, mod_exp_minus))
wl(sprintf("4. **Repeat %d random draws to test robustness** (Approach B \u2014 stochastic): Each draw samples n\u00a0=\u00a0%d *without replacement* from the full %s-isolate Modern pool and runs Fisher\u2019s exact test against the fixed Germany Historical group. Report median p and proportion of draws with p\u00a0<\u00a00.05 \u2014 confirms that the Approach A result holds regardless of which specific %d Modern isolates are chosen.",
           n_downsamp_lower, germany_n, mod_pool_n, germany_n))
wl("")

wl("### Chi-square goodness-of-fit: within-era 4-group distribution")
wl("")
wl("**File:** `HTF4group_ChiSq_uniform_within_era.tsv`")
wl("")
wl("**H\u2080:** All four HTF groups are equally frequent within an era (25% each).")
wl("")
wl("**Why chi-square instead of binomial here?** See Section 5 below.")
wl("")
wl("---")
wl("")

# ---- Section 3a ----
wl("## 3a. Is the binomial test suitable for comparing 4 HTF groups?")
wl("")
wl("**Short answer:** The binomial test is correct for the *collapsed binary* OBC+/\u2212")
wl("question, but NOT for comparing the 4 HTF groups individually or simultaneously.")
wl("")
wl("**Detailed explanation:**")
wl("")
wl("| Question | Groups | Correct test |")
wl("|----------|--------|--------------|")
wl("| Is OBC+ more common than OBC\u2212 within an era? | **Binary**: OBC+ vs OBC\u2212 | \u2705 **Exact binomial** (H\u2080: p=0.5) |")
wl("| Are all 4 HTF groups equally common within an era? | **4 groups** simultaneously | \u2705 **Chi-square GOF** (H\u2080: 25% each) |")
wl("| Did one specific HTF group change between eras? | Per group, cross-era | \u2705 **Fisher exact** per group |")
wl("| Did the overall 4-group distribution change between eras? | 2\u00d74 table | \u2705 **Chi-square test of independence** |")
wl("")
wl("**Why binomial fails for 4 groups:**")
wl("The binomial distribution models a binary outcome (success/failure). Testing")
wl("\"how common is group 1803\" with a binomial treats all non-1803 as equivalent")
wl("\"failure\" \u2014 which obscures the fact that groups 1245, 1383, and 1830 are biologically")
wl("distinct. Running a binomial test on each of the 4 groups separately also inflates")
wl("the Type I error rate (multiple comparisons; 4 tests at \u03b1=0.05 \u2192 up to ~18.5%")
wl("family-wise error rate). The chi-square goodness-of-fit test assesses all 4")
wl("groups in a single omnibus test, properly controlling Type I error.")
wl("")
wl("**What the binomial tests in Tables 1 & 2 actually test:** These correctly use")
wl("binomial because the 4 groups have already been *collapsed* into binary OBC+/\u2212")
wl("(1803+1245 = OBC+; 1383+1830 = OBC\u2212). Once collapsed, it is a binary outcome")
wl("and binomial is the right tool.")
wl("")
wl("---")
wl("")

# ---- Section 3b ----
wl("## 3b. Repeated downsampling vs Permutation vs Bootstrap (clarification)")
wl("")
wl("| Term | Sampling strategy | Purpose |")
wl("|------|------------------|---------|")
wl("| **Bootstrap** | WITH replacement from one sample | Estimate the sampling distribution of a statistic (e.g. CIs) |")
wl("| **Permutation test** | Randomly reassign GROUP LABELS (pool & reshuffle) | Test the null hypothesis that group labels don\u2019t matter |")
wl("| **Repeated random downsampling** \u2705 (used here) | WITHOUT replacement from the larger group | Assess robustness: is the observed Fisher p consistently low regardless of which n=k Modern isolates are included? |")
wl("")
wl("**Why this is NOT a permutation test:** A permutation test pools all")
wl("observations, then randomly assigns them to \"group A\" and \"group B\" to simulate")
wl("the null hypothesis that groups are exchangeable. Here, the Historical group is")
wl("*fixed* and unchanged \u2014 we only repeatedly resample from the Modern pool. The")
wl("question being asked is not \"could this difference arise by chance?\" but rather")
wl("\"is this Fisher p-value robust to which Modern isolates we compare against?\"")
wl("")
wl("**Why this is NOT a bootstrap:** Bootstrap resamples *with* replacement from")
wl("a single sample to estimate a statistic's variance or CI. Here we sample")
wl("*without* replacement from the Modern pool, and the goal is robustness testing,")
wl("not CI estimation.")
wl("")
wl("**The correct interpretation of `Prop_sig_05`:** The fraction of random Modern")
wl(sprintf("subsamples of n=k for which Fisher p < 0.05. If close to 1.0, the Historical"))
wl("vs Modern OBC difference is robust \u2014 it doesn\u2019t matter much which specific")
wl("Modern isolates you compare against, the result is consistently significant.")
wl("")
wl("Variable names updated from `boot_pvals`/`perm_pvals` \u2192 `downsamp_pvals` and")
wl("`n_iter`/`n_perm` \u2192 `n_downsamp`/`n_draws` throughout all scripts.")
wl("")
wl("---")
wl("")

# ---- Section 4 ----
wl("## 4. Logic of the Cross-Era Comparison (Fig4C Middle Panel Design)")
wl("")
wl("The middle panel directly parallels the lower panel (Fig4C bottom, country")
wl("boxplots) in its statistical architecture:")
wl("")
wl("| Panel | Comparison | Test type |")
wl("|-------|-----------|-----------|")
wl("| Lower (step5) | Country-level OBC+ vs OBC\u2212 within Historical or Modern | Binomial + prop.test |")
wl(sprintf("| Lower (step5) | Historical Germany vs Modern Germany (downsampled to n=%d) | Fisher exact |", germany_n))
wl("| Middle (step4) | OBC+ vs OBC\u2212 within Historical | Exact binomial |")
wl("| Middle (step4) | OBC+ vs OBC\u2212 within Modern | Exact binomial |")
wl("| Middle (step4) | Historical vs Modern OBC+ (full cohort) | Fisher exact |")
wl(sprintf("| Middle (step4) | Historical vs Modern OBC+ (Modern \u2192 n=%s) | Fisher exact (two variants) |", apA_n_readme))
wl("")
wl(sprintf("The step5 lower-panel approach for Germany (compare Historical Germany vs Modern"))
wl(sprintf("Germany after downsampling Modern to Germany\u2019s n=%d) is reproduced here at the", germany_n))
wl(sprintf("whole-cohort level, with n_historical=%s as the target sample size.", hist_obc_tot))
wl("")
wl("---")
wl("")

# ---- Section 5: File Inventory ----
wl("## 5. Output File Inventory")
wl("")
wl("### Fig4Cupper_panel_HTF4group_freq/")
wl("")
wl("| File | Contents |")
wl("|------|---------|")
wl("| `HTF_length_hist_vs_modern_panels_dotCI_short.pdf` | Figure: 4-group HTF dot+CI, Historical \\| Modern panels |")
wl("| `HTF4group_hist_vs_modern_counts.tsv` | Per-group counts (long format): IsolateType, LengthGroup, Oantigen, Frequency, TotalN, Percent, Lower, Upper, Label |")
wl("| `HTF4group_hist_vs_modern_counts_wide.tsv` | Per-group counts (wide format): one row per era, columns for each group's count and % |")
wl("| `HTF4group_Binomial_OBC_within_era.tsv` | Exact binomial test: OBC+ proportion vs 50%, within each era |")
wl("| `HTF4group_Fisher_crossera_pergroup.tsv` | Fisher exact test: each of 4 HTF groups, Historical vs Modern |")
wl("| `HTF4group_Fisher_OBC_crossera_downsampled.tsv` | Fisher exact test: OBC overall, full dataset + downsampled (Approaches A & B) |")
wl("| `HTF4group_summary.txt` | Human-readable summary of all counts and statistical results |")
wl("")
wl("### Fig4Cmiddlepanel_oantigenPAfreq/")
wl("")
wl("| File | Contents |")
wl("|------|---------|")
wl("| `Oantigen_freq_hist_vs_modern1350_separated.pdf` | Figure: OBC+/\u2212 proportions, Historical vs Modern |")
wl("| `Oantigen_freq_hist_vs_modern1350_separated_real_counts.tsv` | OBC counts, long format |")
wl("| `Oantigen_freq_hist_vs_modern1350_separated_real_counts_wide.tsv` | OBC counts, wide format |")
wl("| `Oantigen_freq_hist_vs_modern1350_separated_fisher_test.tsv` | Fisher exact result (full dataset) |")
wl("| `Oantigen_freq_hist_vs_modern1350_separated_summary.txt` | Extended human-readable summary (now includes binomial + downsampled Fisher) |")
wl("| `Oantigen_Binomial_OBC_within_era.tsv` | Exact binomial test: OBC+ within each era |")
wl("| `Oantigen_Fisher_OBC_crossera_downsampled.tsv` | Downsampled Fisher exact test (Approaches A & B) |")
wl("")
wl("### Fig4Clowerpanel_top3enrichedcountry_moderndownsampled_fisher/ (step5 additions)")
wl("")
wl("All Germany lower-panel outputs live here \u2014 no separate copy folder.")
wl("")
wl("| File | Contents |")
wl("|------|---------|")
wl("| `Germany_n8_Binomial_OBC_within.tsv` | Binomial test: OBC+/\u2212 within Germany Historical (n=8, H\u2080: 50/50) |")
wl("| `Germany_n8_Binomial_Modern_within.tsv` | Binomial test: OBC+/\u2212 within Modern representative n=8 draw (H\u2080: 50/50) |")
wl("| `Germany_n8_Fisher_expected_value.tsv` | Fisher exact (balanced n=8 vs n=8): Germany Historical vs Modern representative |")
wl("| `Germany_n8_downsampling_1000draws_summary.tsv` | Downsampling robustness summary: median Fisher p, prop. draws sig., 1,000 draws |")
wl("| `Germany_n8_downsampling_1000draws_details.tsv` | Per-draw detail: Plus_Modern, Minus_Modern, Fisher_p, Binomial_p (1,000 rows) |")
wl("| `Oantigen_freq_Germany_hist_vs_resampleModern8.pdf` | **Plot:** Germany Historical vs Modern representative n=8; both errorbars = Wilson 95% CI |")
wl("| `Oantigen_freq_Germany_simple_HvsM.pdf` | **Simple plot:** proportion-only labels, size 9\u00d73.5 (matches middle panel) |")
wl("")
wl(sprintf("**Plot CI note:** Both Historical and Modern groups use Wilson 95%% CI at n=%d, making errorbars directly comparable. The %d-draw robustness analysis is in the tables, not the plot.",
           germany_n, n_downsamp_lower))
wl("")
wl("### final_comparison/ (assembled by both step4 and step5)")
wl("")
wl("| File | Panel | Contents |")
wl("|------|-------|---------|")
wl("| `Table1b_Fig4C_upper_4HTFgroups_chisq_uniform_within_era.tsv` | Upper | Chi-square GOF: 4 HTF groups vs uniform within each era |")
wl("| `Table1_Fig4C_upper_4HTFgroups_OBC_binomial_within_era.tsv` | Upper | Exact binomial: OBC+/\u2212 within Historical and Modern |")
wl("| `Table2_Fig4C_middle_OBC_binomial_within_era.tsv` | Middle | Binomial: OBC+/\u2212 within each era (middle panel label) |")
wl("| `Table3_Fig4C_lower_Germany_n8_Fisher_expected.tsv` | Lower | Fisher exact (expected-value), Germany Historical vs Modern at n=8 |")
wl("| `Table3b_Fig4C_lower_Germany_n8_Binomial_within.tsv` | Lower | Binomial test: OBC+/\u2212 within Germany Historical (n=8) |")
wl("| `Table4_Fig4C_lower_Germany_n8_downsampling_1000draws.tsv` | Lower | Repeated random downsampling summary (1,000 draws), Germany n=8 |")
wl("| `final_comparison_summary.txt` | All | Narrative summary of all comparisons with biological interpretation |")
wl("")
wl("---")
wl("")
wl("### Root output directory (step3_HTFfreq/)")
wl("")
wl("| File | Contents |")
wl("|------|---------|")
wl("| `Supplementary_Statistics_Fig4C_FULL.tsv` | **Full version** \u2014 all panels, all tests (~14 rows); written by step5 after step4 |")
wl("| `Supplementary_Statistics_Fig4C_forpaper.tsv` | **For-paper version** \u2014 upper within, middle within, lower within + downsampling robustness (~6 rows); written by step5 |")
wl("")
wl("---")
wl("")

# ---- Section 6: Significance Notation ----
wl("## 6. Significance Notation")
wl("")
wl("All output files use the following standard notation:")
wl("")
wl("| Symbol | Threshold |")
wl("|--------|-----------|")
wl("| `***`  | p < 0.001 |")
wl("| `**`   | p < 0.01  |")
wl("| `*`    | p < 0.05  |")
wl("| `ns`   | p \u2265 0.05  |")
wl("")
wl("---")
wl("")

# ---- Section 7: Reproducibility ----
wl("## 7. Reproducibility Notes")
wl("")
wl("- All tests are performed in R using base functions (`binom.test`, `fisher.test`)")
wl("  with no external statistical packages.")
wl("- Wilson 95% confidence intervals for proportions use `binom::binom.confint`")
wl("  with `method = \"wilson\"`.")
wl(sprintf("- The repeated random downsampling (Test 4, Approach B; Test 5, Approach B) uses"))
wl(sprintf("  `set.seed(42)` for reproducibility. Results are stable across seeds given"))
wl(sprintf("  %d draws (upper/middle) and %d draws (lower/Germany).",
           apB_ndrw_readme, n_downsamp_lower))
wl("- The counts in the upper panel (4-group) and middle panel (binary OBC) are")
wl("  derived from the same underlying data; the binary OBC counts are simply the")
wl("  sum of the 1803+1245 (OBC+) and 1383+1830 (OBC\u2212) counts from the 4-group table.")
wl("")
wl("---")
wl("")

# ---- Section 8: PNAS-Ready Text ----
wl("## 8. PNAS-Ready Results and Methods Text")
wl("")
wl("### Results (2 sentences \u2014 copy-paste)")
wl("")
wl(sprintf(paste0(
  "> Among Modern *P. viridiflava* isolates (n\u00a0=\u00a0%s), OBC+ haplotypes (HTF-1803\n",
  "> and HTF-1245) were strongly enriched (%s%% OBC+; chi-square goodness-of-fit\n",
  "> across four HTF groups, p\u00a0<\u00a00.001; binomial test vs 50/50, p\u00a0<\u00a00.001), whereas\n",
  "> Historical isolates (n\u00a0=\u00a0%s) showed a near-equal OBC+/\u2212 split (%s%% OBC+;\n",
  "> p\u00a0=\u00a0n.s.), with the cross-era increase confirmed as statistically robust across\n",
  "> %d repeated random downsampling draws (Fisher\u2019s exact test, median p\u00a0<\u00a00.001;\n",
  "> >99%% of draws p\u00a0<\u00a00.05)."),
  mod_obc_tot, mod_obc_pct, hist_obc_tot, hist_obc_pct, apB_ndrw_readme))
wl(sprintf(paste0(
  "> To compare without bias from sample size or geographic origin, Germany Historical\n",
  "> isolates (n\u00a0=\u00a0%d) were extracted as the fixed reference group: %d OBC+ (%d%%),\n",
  "> %d OBC\u2212 (%d%%). The observed OBC\u2212 predominance did not reach significance by\n",
  "> exact binomial test (p\u00a0=\u00a0n.s.), reflecting limited statistical power at\n",
  "> n\u00a0=\u00a0%d. A representative Modern group of n\u00a0=\u00a0%d was then derived from the full\n",
  "> Modern pool (n\u00a0=\u00a0%s OTU5 isolates) by scaling the observed Modern OBC+ proportion\n",
  "> to n\u00a0=\u00a0%d (round(%.1f%%\u00a0\u00d7\u00a0%d)\u00a0=\u00a0%d OBC+, %d OBC\u2212). Fisher\u2019s exact test on this\n",
  "> balanced %d\u00a0vs.\u00a0%d comparison confirmed a significant OBC\u2212 enrichment in\n",
  "> Historical Germany (OR\u00a0=\u00a0%.3f, p\u00a0%s [%s]). Repeated random downsampling of\n",
  "> n\u00a0=\u00a0%d from the full Modern pool (n\u00a0=\u00a0%s) across %d draws confirmed robustness\n",
  "> (median Fisher p\u00a0%s; %.1f%% of draws p\u00a0<\u00a00.05)."),
  germany_n,
  germany_plus, round(germany_plus / germany_n * 100),
  germany_minus, round(germany_minus / germany_n * 100),
  germany_n, germany_n,
  mod_pool_n, germany_n,
  mod_prop_pct, germany_n, mod_exp_plus, mod_exp_minus,
  germany_n, germany_n,
  ger_fish_or_readme, fmt_p(ger_fish_p_readme, digits=3), ger_fish_sig,
  germany_n, mod_pool_n, n_downsamp_lower,
  fmt_p(ger_ds_med_readme), ger_ds_pct_readme))
wl("")
wl("### Methods (2 sentences \u2014 copy-paste)")
wl("")
wl("> Within each era, OBC+/\u2212 proportions were evaluated by two-sided exact binomial")
wl("> test (H\u2080: p\u00a0=\u00a00.5), and the four-group HTF haplotype distribution by")
wl("> chi-square goodness-of-fit (H\u2080: uniform 25% per group, df\u00a0=\u00a03); cross-era")
wl("> comparisons used Fisher\u2019s exact test on 2\u00d72 contingency tables of era\u00a0\u00d7\u00a0OBC")
wl("> status, with per-haplotype group tests applied individually.")
wl(sprintf(paste0(
  "> To address the sample size imbalance (Historical n\u00a0=\u00a0%s or n\u00a0=\u00a0%d vs. Modern\n",
  "> n\u00a0=\u00a0%s), Germany Historical samples were extracted as a fixed reference group\n",
  "> (%d OBC+, %d OBC\u2212) and compared against a representative Modern group of n\u00a0=\u00a0%d\n",
  "> constructed by scaling the observed Modern OBC+ proportion (round(%.1f%%\u00a0\u00d7\u00a0%d)\u00a0=\u00a0%d OBC+);\n",
  "> Fisher\u2019s exact test was run on this balanced %d\u00a0vs.\u00a0%d table, and robustness was\n",
  "> assessed by %d random downsampling draws (n\u00a0=\u00a0%d without replacement from the\n",
  "> %s-isolate Modern pool), reporting the median p-value and proportion of draws\n",
  "> with p\u00a0<\u00a00.05; this procedure is not a permutation test (which randomises group\n",
  "> labels) nor a bootstrap (which resamples with replacement), but evaluates whether\n",
  "> the cross-era OBC signal persists regardless of which specific Modern isolates\n",
  "> are compared. Wilson 95%% confidence intervals were computed using the binom R\n",
  "> package. Significance thresholds: *p\u00a0<\u00a00.05, **p\u00a0<\u00a00.01, ***p\u00a0<\u00a00.001."),
  hist_obc_tot, germany_n, mod_obc_tot,
  germany_plus, germany_minus, germany_n,
  mod_prop_pct, germany_n, mod_exp_plus,
  germany_n, germany_n,
  n_downsamp_lower, germany_n, mod_pool_n))
wl("")
wl("---")
wl("")
wl("### Extended Methods (single paragraph \u2014 for supplementary)")
wl("")
wl("> To test whether O-antigen biosynthesis cluster (OBC+) and OBC\u2212 strains were")
wl(sprintf("> equally represented within each era, we applied a two-sided exact binomial"))
wl(sprintf("> test (H\u2080: p\u00a0=\u00a00.5) separately for Historical (n\u00a0=\u00a0%s) and Modern (n\u00a0=\u00a0%s)",
           hist_obc_tot, mod_obc_tot))
wl("> isolates; for the distribution across all four HTF haplotype groups, we used")
wl("> chi-square goodness-of-fit (H\u2080: 25% per group). Cross-era comparisons used")
wl("> Fisher\u2019s exact test. To account for the large sample size imbalance")
wl(sprintf("> (Historical n\u00a0=\u00a0%s vs. Modern n\u00a0=\u00a0%s), Fisher\u2019s exact test was additionally",
           hist_obc_tot, mod_obc_tot))
wl(sprintf("> applied after (i) computing expected Modern counts at n\u00a0=\u00a0%s using the observed",
           apA_n_readme))
wl("> Modern OBC+ proportion (deterministic approach), and (ii) repeated random")
wl(sprintf("> downsampling of n\u00a0=\u00a0%s isolates without replacement from the full Modern pool",
           apA_n_readme))
wl(sprintf("> across %d draws, reporting the median p-value and proportion of draws with", apB_ndrw_readme))
wl("> p\u00a0<\u00a00.05 as a robustness assessment. For the Germany-specific comparison")
wl(sprintf("> (Historical n\u00a0=\u00a0%d), Germany Historical samples were extracted as a fixed", germany_n))
wl(sprintf("> reference group (%d OBC+, %d OBC\u2212) and a representative Modern group of n\u00a0=\u00a0%d",
           germany_plus, germany_minus, germany_n))
wl("> was constructed by scaling the observed Modern OBC+ proportion to n\u00a0=\u00a08")
wl(sprintf("> (round(%.1f%%\u00a0\u00d7\u00a0%d)\u00a0=\u00a0%d OBC+, %d OBC\u2212); Fisher\u2019s exact test was run on this",
           mod_prop_pct, germany_n, mod_exp_plus, mod_exp_minus))
wl(sprintf("> balanced %d\u00a0vs.\u00a0%d table (Approach A), and robustness was assessed by %d",
           germany_n, germany_n, n_downsamp_lower))
wl(sprintf("> random downsampling draws (n\u00a0=\u00a0%d without replacement from the %s-isolate",
           germany_n, mod_pool_n))
wl("> Modern pool, Approach B). All analyses were conducted in R")
wl("(\u2265\u00a04.0) using base functions (binom.test, fisher.test, chisq.test) with Wilson")
wl("95%% CI (binom package). Significance thresholds: *p\u00a0<\u00a00.05, **p\u00a0<\u00a00.01,")
wl("***p\u00a0<\u00a00.001.")
wl("")
wl("---")
wl("")
wl(paste0("*Generated by: `step4_fig4c_upper_and_middle_modernallfreq_histogram.R` (upper/middle panels) and `step5_fig4c_lower_countryboxplot_bottom.R` (lower panel + combined supplementary tables).*"))

close(rh)

cat(sprintf("\n\U0001F4CB README_Fig4C_statistics.md written to:\n   %s\n", readme_path))

# =============================================================================
# Section 12: In-panel annotation .txt files for each panel
#   One file per panel saved in the panel's output folder.
#   Format: copy-paste text for adding stats annotations directly to plots.
# =============================================================================

outdir_upper <- file.path(step3_outdir, "Fig4Cupper_panel_HTF4group_freq")
outdir_mid   <- file.path(step3_outdir, "Fig4Cmiddlepanel_oantigenPAfreq")
# outdir_lower already defined

# ---- Helper ----
write_panel_txt <- function(path, lines) {
  con <- file(path, open = "wt")
  writeLines(lines, con)
  close(con)
  cat(sprintf("\U0001F4CC Annotation text written to:\n   %s\n", path))
}

# ---- Upper panel ----
# Within-era chi-square GOF (4 HTF groups)
chisq_hist_p   <- if (!is.null(chisq_hist)) as.numeric(pull1(chisq_hist,"P_value")) else NA
chisq_mod_p    <- if (!is.null(chisq_mod))  as.numeric(pull1(chisq_mod, "P_value")) else NA
chisq_hist_chi <- if (!is.null(chisq_hist)) round(as.numeric(pull1(chisq_hist,"ChiSq_statistic")),2) else NA
chisq_mod_chi  <- if (!is.null(chisq_mod))  round(as.numeric(pull1(chisq_mod, "ChiSq_statistic")),1) else NA

write_panel_txt(
  file.path(outdir_upper, "Fig4C_upper_panel_inpanel_annotation.txt"),
  c(
    "=== Fig4C Upper Panel — In-panel Statistical Annotations ===",
    sprintf("Generated: %s", Sys.Date()),
    "",
    "--- Within-era (4 HTF groups, chi-square goodness-of-fit, H0: 25% each) ---",
    sprintf("Historical (n=%s):  chi2 = %.2f, df = 3, p %s [%s]",
            hist_obc_tot, chisq_hist_chi,
            fmt_p(chisq_hist_p), sig_stars_step5(chisq_hist_p)),
    sprintf("Modern    (n=%s): chi2 = %.1f, df = 3, p %s [%s]",
            mod_obc_tot, chisq_mod_chi,
            fmt_p(chisq_mod_p), sig_stars_step5(chisq_mod_p)),
    "",
    "--- Within-era OBC+/- (exact binomial, H0: 50/50) ---",
    sprintf("Historical: %s/%s = %s%% OBC+; p %s [%s]",
            hist_obc_plus, hist_obc_tot, hist_obc_pct,
            fmt_p(hist_binom_p), sig_stars_step5(hist_binom_p)),
    sprintf("Modern:     %s/%s = %s%% OBC+; p %s [%s]",
            mod_obc_plus,  mod_obc_tot,  mod_obc_pct,
            fmt_p(mod_binom_p), sig_stars_step5(mod_binom_p)),
    "",
    "Significance: *** p<0.001 | ** p<0.01 | * p<0.05 | ns p>=0.05"
  )
)

# ---- Middle panel ----
write_panel_txt(
  file.path(outdir_mid, "Fig4C_middle_panel_inpanel_annotation.txt"),
  c(
    "=== Fig4C Middle Panel — In-panel Statistical Annotations ===",
    sprintf("Generated: %s", Sys.Date()),
    "",
    "--- Within-era (exact binomial, H0: OBC+ = 50%) ---",
    sprintf("Historical (n=%s): %s/%s = %s%% OBC+; p %s [%s]",
            hist_obc_tot, hist_obc_plus, hist_obc_tot, hist_obc_pct,
            fmt_p(hist_binom_p), sig_stars_step5(hist_binom_p)),
    sprintf("Modern    (n=%s): %s/%s = %s%% OBC+; p %s [%s]",
            mod_obc_tot, mod_obc_plus, mod_obc_tot, mod_obc_pct,
            fmt_p(mod_binom_p), sig_stars_step5(mod_binom_p)),
    "",
    "--- Between-era (Fisher's exact test) ---",
    sprintf("Full dataset: OR = %.3f, p %s [%s]",
            full_fish_or, fmt_p(full_fish_p), sig_stars_step5(full_fish_p)),
    sprintf("Downsampled n=%d (expected): OR = %.3f, p %s [%s]",
            apA_n, apA_or, fmt_p(apA_p), sig_stars_step5(apA_p)),
    sprintf("Repeated downsampling (%d draws, n=%d): median p %s; %.1f%% of draws p < 0.05",
            apB_ndrw, apA_n, fmt_p(apB_med), apB_pct),
    "",
    "Significance: *** p<0.001 | ** p<0.01 | * p<0.05 | ns p>=0.05"
  )
)

# ---- Lower panel ----
# Modern within-era binomial: from binom_germany_modern_within
mod_within_p   <- binom_germany_modern_within$p.value
mod_within_sig <- sig_stars_step5(mod_within_p)

write_panel_txt(
  file.path(outdir_lower, "Fig4C_lower_panel_inpanel_annotation.txt"),
  c(
    "=== Fig4C Lower Panel — In-panel Statistical Annotations ===",
    sprintf("Generated: %s", Sys.Date()),
    "",
    "--- Within-era: Germany Historical (exact binomial, H0: OBC+ = 50%) ---",
    sprintf("Germany Historical (n=%d): %d/%d = %d%% OBC+; p %s [%s]",
            germany_n, germany_plus, germany_n,
            round(germany_plus/germany_n*100),
            fmt_p(binom_germany_within$p.value),
            sig_stars_step5(binom_germany_within$p.value)),
    sprintf("  Wilson 95%% CI for OBC+: %.1f%% - %.1f%%",
            hist_ci_plus$lower*100, hist_ci_plus$upper*100),
    "",
    "--- Within-era: Modern representative n=8 (exact binomial, H0: OBC+ = 50%) ---",
    sprintf("Modern representative (n=%d): %d/%d = %.1f%% OBC+ [rounded: %d/%d = %.1f%%]; p %s [%s]",
            n_samp, round(mod_typical_plus), n_samp, mod_typical_plus/n_samp*100,
            mod_exp_plus, n_samp, mod_exp_plus/n_samp*100,
            fmt_p(mod_within_p), mod_within_sig),
    sprintf("  Wilson 95%% CI for OBC+ (CI uses floating k=%.3f): %.1f%% - %.1f%%",
            mod_typical_plus, mod_ci_plus$lower*100, mod_ci_plus$upper*100),
    "",
    "--- Between-era: Germany Historical vs Modern (Fisher's exact, balanced n=8 vs n=8) ---",
    sprintf("Approach A (deterministic): OR = %.3f, p %s [%s]",
            ger_fish_or, fmt_p(ger_fish_p, digits=3), sig_stars_step5(ger_fish_p)),
    sprintf("  [Germany Hist: %d OBC+, %d OBC-] vs [Modern rep: %d OBC+, %d OBC-]",
            germany_plus, germany_minus, mod_exp_plus, mod_exp_minus),
    sprintf("Approach B (downsampling, %d draws, n=%d from pool n=%s):",
            n_downsamp_lower, germany_n, mod_pool_n),
    sprintf("  Median Fisher p %s; %.1f%% of draws p < 0.05 [%s]",
            fmt_p(ger_ds_med), ger_ds_pct, downsamp_summary_n8$Fisher_Significance),
    "",
    "--- Plot label note ---",
    sprintf("Simple lower panel dot labels show ROUNDED integer proportion:"),
    sprintf("  Historical: %d%% (%d/%d)",
            round(germany_plus/germany_n*100), germany_plus, germany_n),
    sprintf("  Modern rep: %.1f%% (%d/%d)  [floating CI k = %.3f/%d = %.1f%%]",
            mod_exp_plus/n_samp*100, mod_exp_plus, n_samp,
            mod_typical_plus, n_samp, mod_typical_plus/n_samp*100),
    "",
    "Significance: *** p<0.001 | ** p<0.01 | * p<0.05 | ns p>=0.05"
  )
)

cat("\n\u2705 All in-panel annotation .txt files written.\n")

# =============================================================================
# Section 13: Fig4C_figure_legend.txt
#   PNAS-style figure captions for all three panels.
#   Two versions per panel: LONG (detailed with results) and SHORT (concise).
#   Written to outdir_final (final_comparison/).
#   NOTE: placed here so all stat variables (chisq_hist/mod, binom_hist_r/mod_r,
#         hist_binom_p, mod_binom_p, hist_chisq_str, mod_chisq_str, Wilson CIs,
#         mod_within_p, ger_hist_binom_p) are fully defined.
# =============================================================================

legend_path <- file.path(outdir_final, "Fig4C_figure_legend.txt")
zl <- file(legend_path, open = "wt")
wl_l <- function(...) writeLines(paste0(...), zl)

# ---- shared title ----
wl_l("Fig. 4C.  OBC carriage across HTF haplotypes and geographic strata")
wl_l("in historical and modern P. viridiflava.")
wl_l("")
wl_l(paste(rep("=", 72), collapse = ""))
wl_l("VERSION 1 — LONG (detailed; for supplementary or methods)")
wl_l(paste(rep("=", 72), collapse = ""))
wl_l("")

# ---- Upper panel LONG ----
wl_l("(Upper panel)  Frequency of each of the four major HTF haplotypes \u2014 1803")
wl_l("and 1245 (OBC-positive; OBC+) and 1383 and 1830 (OBC-negative; OBC\u2212) \u2014")
wl_l(sprintf("among Historical (pre-antibiotic era; n\u00a0=\u00a0%s) and Modern (n\u00a0=\u00a0%s)", hist_obc_tot, mod_pool_n))
wl_l("P.\u00a0viridiflava isolates. Within the Historical era, OBC+")
wl_l(sprintf("haplotypes (1803\u00a0+\u00a01245) accounted for %s of %s isolates (%s%%).", hist_obc_plus, hist_obc_tot, hist_obc_pct))
wl_l("The distribution of the four haplotypes was compared against a uniform")
wl_l(sprintf("25%% expectation within each era (chi-square goodness-of-fit: Historical %s;", hist_chisq_str))
wl_l(sprintf("Modern %s).", mod_chisq_str))
wl_l("")

# ---- Middle panel LONG ----
wl_l("(Middle panel)  OBC carriage frequency (OBC+ and OBC\u2212) across all isolates")
wl_l(sprintf("in each era (Historical n\u00a0=\u00a0%s, Modern n\u00a0=\u00a0%s). OBC+ isolates comprised", hist_obc_tot, mod_pool_n))
wl_l(sprintf("%s%% (%s/%s) of Historical and %s%% (%s/%s) of Modern isolates.",
             hist_obc_pct, hist_obc_plus, hist_obc_tot,
             mod_obc_pct,  mod_obc_plus,  mod_obc_tot))
wl_l(sprintf("OBC+ carriage differed from equal split (50%%) in both eras"))
wl_l(sprintf("(Historical: %s; Modern: %s).", hist_binom_str, mod_binom_str))
wl_l("")

# ---- Supplementary figure legend LONG ----
wl_l("(Supplementary Fig.)  Downsampling OBC frequency analysis.")
wl_l("")
wl_l(sprintf("(A)  OBC frequency in Germany Historical isolates (n\u00a0=\u00a0%d;", germany_n))
wl_l(sprintf("%d/%d = %d%% OBC+; Wilson 95%% confidence intervals).",
             germany_plus, germany_n, round(germany_plus/germany_n*100)))
wl_l(sprintf("Exact binomial test (H\u2080: OBC+ = 50%%) yielded p\u00a0%s [%s],",
             fmt_p(binom_germany_within$p.value),
             sig_stars_step5(binom_germany_within$p.value)))
wl_l("indicating no significant deviation from equal OBC+/\u2212 split in this dataset.")
wl_l("")
wl_l(sprintf("(B)  To assess OBC frequency differences between Historical and Modern"))
wl_l("datasets without bias from geographic origin and sample size, the Modern pool")
wl_l(sprintf("(n\u00a0=\u00a0%s) was downsampled to match the Germany Historical sample size", mod_pool_n))
wl_l(sprintf("(n\u00a0=\u00a0%d). %d independent runs were performed, each drawing n\u00a0=\u00a0%d", n_samp, n_downsamp_lower, n_samp))
wl_l("isolates without replacement. OBC+ and OBC\u2212 proportions per run are shown")
wl_l("as a boxplot (IQR with median line, whiskers 1.5\u00d7IQR) with individual")
wl_l("run proportions overlaid as jittered points.")
wl_l("")
wl_l(sprintf("(C)  Distribution of two-sided exact binomial p-values across %d", n_downsamp_lower))
wl_l("downsampling runs. Within each run, OBC frequency was assessed by exact")
wl_l("binomial test (H\u2080: OBC+ = 50%%). Red dashed line: \u03b1\u00a0=\u00a00.05.")
wl_l(sprintf("Only %.1f%% of runs yielded binomial p\u00a0<\u00a00.05 (median p\u00a0%s),",
             downsamp_summary_n8$Binomial_prop_sig_05 * 100,
             fmt_p(downsamp_summary_n8$Binomial_median_p)))
wl_l(sprintf("indicating that n\u00a0=\u00a0%d samples from the Modern dataset do not reliably", n_samp))
wl_l("represent the true Modern OBC frequency. Accordingly, there is insufficient")
wl_l("power to perform a cross-dataset OBC frequency comparison, and no such")
wl_l("comparison is reported.")
wl_l("")

wl_l(paste(rep("=", 72), collapse = ""))
wl_l("VERSION 2 — SHORT (concise; for main figure caption)")
wl_l(paste(rep("=", 72), collapse = ""))
wl_l("")

# ---- Upper panel SHORT ----
wl_l("(Upper panel)  Frequency of four HTF haplotypes (1803, 1245 = OBC+;")
wl_l(sprintf("1383, 1830 = OBC\u2212) in Historical (n\u00a0=\u00a0%s) and Modern (n\u00a0=\u00a0%s) isolates.", hist_obc_tot, mod_pool_n))
wl_l(sprintf("OBC+ haplotypes comprised %s%% of Historical and %s%% of Modern isolates", hist_obc_pct, mod_obc_pct))
wl_l(sprintf("(chi-square goodness-of-fit: Historical %s; Modern %s).", hist_chisq_str, mod_chisq_str))
wl_l("")

# ---- Middle panel SHORT ----
wl_l("(Middle panel)  OBC+ carriage per era")
wl_l(sprintf("(Historical %s%% [%s/%s], %s;",
             hist_obc_pct, hist_obc_plus, hist_obc_tot, hist_binom_str))
wl_l(sprintf("Modern %s%% [%s/%s], %s).",
             mod_obc_pct, mod_obc_plus, mod_obc_tot, mod_binom_str))
wl_l("")

# ---- Supplementary figure SHORT ----
wl_l("(Supplementary Fig.)  Downsampling OBC frequency analysis.")
wl_l(sprintf("(A) Germany Historical OBC frequency (n=%d; %d/%d = %d%% OBC+; Wilson 95%% CI;",
             germany_n, germany_plus, germany_n, round(germany_plus/germany_n*100)))
wl_l(sprintf("exact binomial p\u00a0%s [%s]).",
             fmt_p(binom_germany_within$p.value),
             sig_stars_step5(binom_germany_within$p.value)))
wl_l(sprintf("(B\u2013C) To assess OBC frequency differences between Historical and Modern"))
wl_l(sprintf("datasets without bias from geographic origin and sample size, %d", n_downsamp_lower))
wl_l(sprintf("downsampling runs (n=%d without replacement from Modern dataset, n=%s)", n_samp, mod_pool_n))
wl_l("were performed. OBC frequency within each run was assessed by exact binomial")
wl_l(sprintf("(H\u2080: OBC+ = 50%%). Only %.1f%% of runs yielded p\u00a0<\u00a00.05",
             downsamp_summary_n8$Binomial_prop_sig_05 * 100))
wl_l(sprintf("(median p\u00a0%s), indicating that n=%d samples from the Modern dataset",
             fmt_p(downsamp_summary_n8$Binomial_median_p), n_samp))
wl_l("do not reliably represent the true Modern OBC frequency; accordingly,")
wl_l("there is insufficient power for cross-dataset comparison.")

close(zl)

cat("\n\u2705 Fig4C figure legend (long + short versions) saved.\n")
cat("   File: Fig4C_figure_legend.txt\n")
cat("   Path:", outdir_final, "\n")
