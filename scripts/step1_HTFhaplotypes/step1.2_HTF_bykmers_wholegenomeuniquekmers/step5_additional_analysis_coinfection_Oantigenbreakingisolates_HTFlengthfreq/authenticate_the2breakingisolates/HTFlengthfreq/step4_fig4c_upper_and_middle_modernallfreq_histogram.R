#!/usr/bin/env Rscript
# =============================================================================
# Fig4C Upper & Middle panels: HTF haplotype + O-antigen frequency histograms
# with comprehensive statistical tests for PNAS supplementary material
#
# Outputs:
#   Fig4Cupper_panel_HTF4group_freq/
#     ├── HTF_length_hist_vs_modern_panels_dotCI_short.pdf  (plot)
#     ├── HTF4group_hist_vs_modern_counts.tsv               (per-group counts, long format)
#     ├── HTF4group_hist_vs_modern_counts_wide.tsv          (per-group counts, wide format)
#     ├── HTF4group_Binomial_OBC_within_era.tsv             (binomial test OBC+/- within era)
#     ├── HTF4group_Fisher_crossera_pergroup.tsv            (Fisher per-HTF-group, H vs M)
#     ├── HTF4group_Fisher_OBC_crossera_downsampled.tsv     (Fisher OBC overall + downsampled)
#     └── HTF4group_summary.txt                            (human-readable summary)
#
#   Fig4Cmiddlepanel_oantigenPAfreq/
#     ├── Oantigen_freq_hist_vs_modern1350_separated.pdf    (plot, existing)
#     ├── Oantigen_freq_hist_vs_modern1350_separated_real_counts.tsv     (existing)
#     ├── Oantigen_freq_hist_vs_modern1350_separated_real_counts_wide.tsv (existing)
#     ├── Oantigen_freq_hist_vs_modern1350_separated_fisher_test.tsv     (existing)
#     ├── Oantigen_freq_hist_vs_modern1350_separated_summary.txt         (existing)
#     ├── Oantigen_Binomial_OBC_within_era.tsv             (NEW: binomial within-era)
#     └── Oantigen_Fisher_OBC_crossera_downsampled.tsv     (NEW: Fisher + downsampled)
#
# Statistical tests applied:
#   1. Binomial test (two-sided, H0: p = 0.5):
#        Within each era (Historical, Modern), test whether O-antigen+ proportion
#        significantly departs from a 50/50 split.
#   2. Fisher's exact test (per HTF group, cross-era):
#        For each of the 4 HTF haplotype groups, compare Historical vs Modern
#        proportions using a 2×2 contingency table.
#   3. Fisher's exact test (O-antigen overall, cross-era, full dataset):
#        Compare total O-antigen+/- counts between Historical and Modern (already in
#        the existing middle panel output; reproduced here for completeness).
#   4. Fisher's exact test (O-antigen overall, cross-era, with repeated downsampling):
#        To make the cross-era comparison fair given the large sample size imbalance
#        (Historical n=38 vs Modern n~1312), Modern is randomly subsampled to n=38
#        across 1,000 repeated draws. Fisher's exact test is run on each draw.
#        The distribution of p-values assesses whether the result is robust to
#        which specific Modern isolates happen to be included. This is a repeated
#        random downsampling robustness analysis — NOT a permutation test (which
#        randomises labels) and NOT a bootstrap (which resamples with replacement).
# =============================================================================

# ---- Libraries ----
library(ggplot2)
library(readr)
library(dplyr)
library(data.table)
library(binom)
library(patchwork)
library(tidyr)
library(scales)

# =============================================================================
# INPUT / OUTPUT PATHS  — absolute paths, verified against directory listings
# =============================================================================
base_dir    <- "/Users/jiajuncui/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57"

# ---- Script-side inputs (raw data) ----
scripts_htf <- file.path(base_dir,
  "scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers",
  "step5_additional_analysis_coinfection_Oantigenbreakingisolates_HTFlengthfreq/HTFlengthfreq")

infile      <- file.path(scripts_htf, "HP12_combined_lengths_sorted.txt")
otu5_file   <- file.path(scripts_htf, "1355OTU5andp8A2.txt")

# ---- Results-side inputs (kmer proportion tables) ----
results_kmers <- file.path(base_dir,
  "results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers")

modern_path <- file.path(results_kmers, "newsummarymodern57_tailocin_kmer_propnorm.tsv")
hist_path   <- file.path(results_kmers, "newsummaryhistorical49_tailocin_kmer_propnorm_05sd.tsv")

# ---- Base output directory (step3_HTFfreq results folder) ----
base_outdir <- file.path(base_dir,
  "results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step4_additionalanalysis/step3_HTFfreq")

# Subfolders
outdir_upper  <- file.path(base_outdir, "Fig4Cupper_panel_HTF4group_freq")
outdir_middle <- file.path(base_outdir, "Fig4Cmiddlepanel_oantigenPAfreq")

dir.create(outdir_upper,  showWarnings = FALSE, recursive = TRUE)
dir.create(outdir_middle, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# COLOURS AND FACTOR LEVELS
# =============================================================================
length_colors  <- c("1830" = "#f6d6ff",
                    "1383" = "#638ccc",
                    "1803" = "#800233",
                    "1245" = "#f9d42a")
length_groups  <- rev(c("1803","1245","1383","1830"))   # display order
oantigen_colors <- c("O-antigen-" = "#9e9ac8", "O-antigen+" = "#c94f7c")

# Helper: significance stars
sig_stars <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

# =============================================================================
# SECTION A: MODERN DATA (OTU5 strains, HP12 combined lengths)
# =============================================================================
df_raw       <- read_tsv(infile, show_col_types = FALSE)
otu5_strains <- fread(otu5_file, header = FALSE)$V1

df_raw <- df_raw %>% filter(Sample %in% otu5_strains)

# Summarise per-length counts (Modern, all OTU5)
df_summary <- df_raw %>%
  filter(!is.na(Length)) %>%
  count(Length) %>%
  mutate(
    prop       = n / sum(n) * 100,
    prop_label = paste0(round(prop, 1), "%"),
    Length     = factor(as.character(Length), levels = rev(length_groups))
  ) %>%
  filter(!is.na(Length))

total_modern <- sum(df_summary$n)   # actual OTU5 count (e.g. 1312)

# Annotate O-antigen status for Modern summary
df_summary <- df_summary %>%
  mutate(Oantigen = case_when(
    Length %in% c("1803","1245") ~ "O-antigen+",
    Length %in% c("1383","1830") ~ "O-antigen-",
    TRUE ~ "Unknown"
  ))

# =============================================================================
# SECTION B: HISTORICAL + MODERN COMBINED (kmer-based dominant haplotype tables)
# =============================================================================
modern_kmer   <- read.table(modern_path,   header = TRUE, sep = "\t")
historical_kmer <- read.table(hist_path,   header = TRUE, sep = "\t")

# Remove QC-flagged isolates
modern_kmer   <- modern_kmer   %>% filter(!Isolate %in% c("p12.F2","p13.C7","p6.A10","p9.C4"))
historical_kmer <- historical_kmer %>% filter(Isolate != "64.GBR_1933b_S36")

modern_kmer$IsolateType   <- "Modern"
historical_kmer$IsolateType <- "Historical"
df_all <- bind_rows(modern_kmer, historical_kmer)

# Exclude low-quality isolates
exclude <- c("HB0828","HB0863","PL0066","PL0108","PL0203","PL0258",
             "64.GBR_1933b_S36","PL0065","PL0026","PL0027","PL0053")
df_all  <- df_all %>% filter(!Isolate %in% exclude)

# Assign length groups
df_all <- df_all %>%
  mutate(LengthGroup = case_when(
    Reference == "HTF_p7.G11"  ~ "1830",
    Reference == "HTF_p25.A12" ~ "1383",
    Reference %in% c("HTF_p5.D5","HTF_p23.B8","HTF_p26.D6","HTF_p25.C2") ~ "1803",
    Reference == "HTF_p21.F9"  ~ "1245",
    TRUE ~ "Unknown"
  )) %>%
  filter(LengthGroup != "Unknown")

# Keep dominant haplotype per isolate
dominant_df <- df_all %>%
  group_by(Isolate) %>%
  filter(PropNorm == max(PropNorm, na.rm = TRUE)) %>%
  ungroup()

# Per-era total n
total_counts <- dominant_df %>%
  distinct(Isolate, IsolateType) %>%
  count(IsolateType, name = "TotalN")

# Per-group per-era frequency table with Wilson 95% CI
freq_table_length <- dominant_df %>%
  group_by(LengthGroup, IsolateType) %>%
  summarise(Frequency = n(), .groups = "drop") %>%
  left_join(total_counts, by = "IsolateType") %>%
  mutate(
    Proportion = Frequency / TotalN,
    ci     = binom.confint(Frequency, TotalN, method = "wilson"),
    Lower  = ci$lower,
    Upper  = ci$upper
  ) %>%
  select(LengthGroup, IsolateType, Frequency, TotalN, Proportion, Lower, Upper)

# Add O-antigen annotation
freq_table_length <- freq_table_length %>%
  mutate(
    Oantigen    = case_when(
      LengthGroup %in% c("1803","1245") ~ "O-antigen+",
      LengthGroup %in% c("1383","1830") ~ "O-antigen-"
    ),
    LengthGroup = factor(LengthGroup, levels = length_groups)
  )

# =============================================================================
# SECTION C: BUILD PLOT DATAFRAMES
# =============================================================================

# ---- C1. Historical (4-group dot-CI plot) ----
df_hist4 <- freq_table_length %>%
  filter(IsolateType == "Historical") %>%
  mutate(
    Proportion = Proportion * 100,
    Lower      = Lower * 100,
    Upper      = Upper * 100
  )

# ---- C2. Modern (4-group dot-CI plot, from HP12/OTU5 summary) ----
df_mod4 <- df_summary %>%
  mutate(LengthGroup = factor(as.character(Length), levels = length_groups)) %>%
  group_by(Oantigen, LengthGroup) %>%
  summarise(Frequency = sum(n), .groups = "drop") %>%
  mutate(
    IsolateType = "Modern",
    TotalN      = total_modern,
    ci          = binom.confint(Frequency, TotalN, method = "wilson"),
    Proportion  = ci$mean * 100,
    Lower       = ci$lower * 100,
    Upper       = ci$upper * 100
  ) %>%
  select(IsolateType, Oantigen, LengthGroup, Frequency, TotalN, Proportion, Lower, Upper)

# Combined 4-group dataframe
df_combined4 <- bind_rows(df_hist4, df_mod4) %>%
  mutate(IsolateType = factor(IsolateType, levels = c("Historical","Modern")))

# ---- C3. Binary O-antigen (Historical vs Modern, for middle panel) ----
df_hist_obc <- freq_table_length %>%
  filter(IsolateType == "Historical") %>%
  group_by(IsolateType, Oantigen) %>%
  summarise(Frequency = sum(Frequency), TotalN = unique(TotalN), .groups = "drop") %>%
  mutate(
    Proportion = Frequency / TotalN,
    ci     = binom.confint(Frequency, TotalN, method = "wilson"),
    Lower  = ci$lower, Upper = ci$upper
  ) %>%
  select(IsolateType, Oantigen, Frequency, TotalN, Proportion, Lower, Upper)

df_mod_obc <- df_summary %>%
  group_by(Oantigen) %>%
  summarise(Frequency = sum(n), TotalN = total_modern, .groups = "drop") %>%
  mutate(
    Proportion = Frequency / TotalN,
    ci     = binom.confint(Frequency, TotalN, method = "wilson"),
    Lower  = ci$lower, Upper = ci$upper,
    IsolateType = "Modern"
  ) %>%
  select(IsolateType, Oantigen, Frequency, TotalN, Proportion, Lower, Upper)

df_combined <- bind_rows(df_hist_obc, df_mod_obc) %>%
  mutate(
    IsolateType = factor(IsolateType, levels = c("Historical","Modern")),
    Oantigen    = factor(Oantigen, levels = c("O-antigen-","O-antigen+")),
    Label       = paste0(Frequency,"/",TotalN," (",round(Proportion*100,1),"%)")
  )

# Convenience scalars for statistical sections
hist_plus  <- df_combined$Frequency[df_combined$IsolateType=="Historical" & df_combined$Oantigen=="O-antigen+"]
hist_minus <- df_combined$Frequency[df_combined$IsolateType=="Historical" & df_combined$Oantigen=="O-antigen-"]
hist_total <- hist_plus + hist_minus
mod_plus   <- df_combined$Frequency[df_combined$IsolateType=="Modern"     & df_combined$Oantigen=="O-antigen+"]
mod_minus  <- df_combined$Frequency[df_combined$IsolateType=="Modern"     & df_combined$Oantigen=="O-antigen-"]
mod_total  <- mod_plus + mod_minus

# =============================================================================
# SECTION D: PLOTS
# =============================================================================

# ---- D1. Upper panel: 4-group Historical vs Modern dot+CI ----
plot_panel_4grp <- function(df, panel_title) {
  dodge <- position_dodge(width = 0.8)
  ggplot(df, aes(x = Oantigen, y = Proportion, color = LengthGroup)) +
    geom_point(position = dodge, size = 4) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper, group = LengthGroup),
                  width = 0.2, position = dodge, linewidth = 0.7,
                  show.legend = FALSE) +
    geom_text(aes(label = paste0(round(Proportion,1),"%"), group = LengthGroup),
              position = dodge, vjust = -2, size = 5, fontface = "bold",
              color = "black") +
    scale_color_manual(values = length_colors,
                       guide = guide_legend(override.aes = list(size=5))) +
    scale_y_continuous(labels = function(x) paste0(x,"%"), limits = c(0,50)) +
    labs(title = panel_title, x = "", y = "Proportion (%)") +
    theme_classic(base_size = 14) +
    theme(plot.title    = element_text(size=15, hjust=0.5, face="bold"),
          legend.position = "top",
          legend.title  = element_blank(),
          legend.key    = element_blank(),
          axis.text     = element_text(size=14, face="bold"),
          axis.title    = element_text(size=14, face="bold"))
}

p_hist4 <- plot_panel_4grp(filter(df_combined4, IsolateType=="Historical"), "Historical")
p_mod4  <- plot_panel_4grp(filter(df_combined4, IsolateType=="Modern"),     "Modern") +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y  = element_blank())

p_upper <- (p_hist4 + p_mod4) + plot_layout(guides = "collect")
ggsave(file.path(outdir_upper, "HTF_length_hist_vs_modern_panels_dotCI_short.pdf"),
       plot = p_upper, width = 14, height = 4, dpi = 600)

# ---- D2. Middle panel: binary O-antigen Historical vs Modern dot+CI ----
dodge <- position_dodge(width = 0.8)
p_middle <- ggplot(df_combined,
                   aes(x = IsolateType, y = Proportion, color = Oantigen)) +
  geom_point(position = dodge, size = 3, show.legend = TRUE) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.1,
                position = dodge, show.legend = FALSE) +
  geom_text(aes(label = paste0(round(Proportion*100,1),"%"), group = Oantigen),
            position = dodge, vjust = -1.0, size = 4, fontface = "bold",
            color = "black", show.legend = FALSE) +
  scale_color_manual(values = oantigen_colors) +
  scale_y_continuous(labels = percent_format(accuracy=1), limits = c(0,1)) +
  labs(x = "", y = "Proportion (%)") +
  theme_classic(base_size = 12) +
  theme(axis.text    = element_text(size=11, face="bold"),
        axis.title   = element_text(size=11, face="bold"),
        legend.title = element_blank(),
        legend.position = "top",
        legend.text  = element_text(size=11, face="bold"))

ggsave(file.path(outdir_middle, "Oantigen_freq_hist_vs_modern1350_separated.pdf"),
       plot = p_middle, width = 9, height = 3.5, dpi = 600)

# =============================================================================
# SECTION E: STATISTICAL TESTS
# =============================================================================

# ---- E1. COUNT TABLES (upper panel) ----
# Long format: one row per (IsolateType × LengthGroup)
upper_counts_long <- df_combined4 %>%
  mutate(
    Percent = round(Proportion, 1),
    Label   = paste0(Frequency, "/", TotalN, " (", Percent, "%)")
  ) %>%
  arrange(IsolateType, Oantigen, LengthGroup) %>%
  select(IsolateType, LengthGroup, Oantigen, Frequency, TotalN, Percent, Lower, Upper, Label)

write.table(upper_counts_long,
            file.path(outdir_upper, "HTF4group_hist_vs_modern_counts.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# Wide format: one row per IsolateType, columns for each group
upper_counts_wide <- upper_counts_long %>%
  select(IsolateType, LengthGroup, Frequency, Percent, TotalN) %>%
  pivot_wider(names_from = LengthGroup,
              values_from = c(Frequency, Percent),
              names_glue = "{LengthGroup}_{.value}") %>%
  arrange(IsolateType)

write.table(upper_counts_wide,
            file.path(outdir_upper, "HTF4group_hist_vs_modern_counts_wide.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# ---- E2. COUNT TABLES (middle panel: existing format) ----
oantigen_count_table <- df_combined %>%
  mutate(Percent = round(Proportion * 100, 1),
         Label   = paste0(Frequency,"/",TotalN," (",Percent,"%)")) %>%
  arrange(IsolateType, Oantigen) %>%
  select(IsolateType, Oantigen, Frequency, TotalN, Percent, Lower, Upper, Label)

write.table(oantigen_count_table,
            file.path(outdir_middle, "Oantigen_freq_hist_vs_modern1350_separated_real_counts.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

oantigen_count_wide <- df_combined %>%
  select(IsolateType, Oantigen, Frequency, TotalN, Proportion) %>%
  mutate(Percent = round(Proportion * 100, 1)) %>%
  select(-Proportion) %>%
  pivot_wider(names_from = Oantigen, values_from = c(Frequency, Percent)) %>%
  arrange(IsolateType)

write.table(oantigen_count_wide,
            file.path(outdir_middle, "Oantigen_freq_hist_vs_modern1350_separated_real_counts_wide.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# ---- E3. BINOMIAL TEST: O-antigen+/- within each era (H0: p_plus = 0.5) ----
#
# Rationale: Within Historical or Modern isolates, test whether O-antigen+
# and O-antigen- strains are equally frequent (H0: 50/50 split).
# A significant result means one group is enriched within that era.
# Method: Two-sided exact binomial test.
#
compute_binomial_tests <- function(plus_h, minus_h, plus_m, minus_m) {
  eras <- list(
    list(era = "Historical", plus = plus_h, minus = minus_h),
    list(era = "Modern",     plus = plus_m, minus = minus_m)
  )
  lapply(eras, function(e) {
    total <- e$plus + e$minus
    bt <- binom.test(e$plus, total, p = 0.5, alternative = "two.sided")
    data.frame(
      IsolateType      = e$era,
      OBC_plus_n       = e$plus,
      OBC_minus_n      = e$minus,
      Total_N          = total,
      Prop_OBC_plus    = round(e$plus / total * 100, 2),
      Test             = "Two-sided exact binomial (H0: p_plus = 0.5)",
      P_value          = bt$p.value,
      CI_lower_plus    = round(bt$conf.int[1] * 100, 2),
      CI_upper_plus    = round(bt$conf.int[2] * 100, 2),
      Significance     = sig_stars(bt$p.value),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
}

binomial_tests <- compute_binomial_tests(hist_plus, hist_minus, mod_plus, mod_minus)

write.table(binomial_tests,
            file.path(outdir_upper,  "HTF4group_Binomial_OBC_within_era.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(binomial_tests,
            file.path(outdir_middle, "Oantigen_Binomial_OBC_within_era.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# ---- E4. FISHER EXACT TEST: per HTF group, Historical vs Modern ----
#
# Rationale: For each of the 4 HTF haplotype groups (1803, 1245, 1383, 1830),
# test whether its frequency differs between Historical and Modern isolates.
# Method: Fisher's exact test on a 2×2 table:
#         [isolates WITH group | isolates WITHOUT group]
#         [  Historical row    |     Modern row        ]
# This identifies which specific haplotypes underwent frequency shifts.
#
fisher_pergroup <- lapply(levels(df_combined4$LengthGroup), function(grp) {
  h <- filter(df_combined4, IsolateType == "Historical", LengthGroup == grp)
  m <- filter(df_combined4, IsolateType == "Modern",     LengthGroup == grp)
  h_freq  <- h$Frequency;   h_total <- h$TotalN
  m_freq  <- round(m$Frequency); m_total <- m$TotalN   # Modern counts rounded to integer

  mat <- matrix(c(h_freq,        h_total - h_freq,
                  m_freq,        m_total - m_freq),
                nrow = 2, byrow = TRUE,
                dimnames = list(c("Historical","Modern"),
                                c("WithGroup","WithoutGroup")))
  ft <- fisher.test(mat)
  data.frame(
    LengthGroup     = grp,
    Oantigen        = ifelse(grp %in% c("1803","1245"), "O-antigen+", "O-antigen-"),
    Hist_count      = h_freq,
    Hist_total      = h_total,
    Hist_pct        = round(h_freq / h_total * 100, 1),
    Mod_count       = m_freq,
    Mod_total       = m_total,
    Mod_pct         = round(m_freq / m_total * 100, 1),
    Odds_ratio      = round(unname(ft$estimate), 4),
    P_value         = ft$p.value,
    Significance    = sig_stars(ft$p.value),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

write.table(fisher_pergroup,
            file.path(outdir_upper, "HTF4group_Fisher_crossera_pergroup.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# ---- E5. FISHER EXACT TEST: O-antigen overall, full dataset ----
#
# Rationale: Test whether the total proportion of O-antigen+ isolates
# differs significantly between Historical and Modern.
# Method: Fisher's exact test on the 2×2 O-antigen+/- × Era table.
# This is the primary cross-era comparison (same test as in Fig4C middle panel).
#
mat_full <- matrix(c(hist_plus, hist_minus, mod_plus, mod_minus),
                   nrow = 2, byrow = TRUE,
                   dimnames = list(c("Historical","Modern"),
                                   c("O-antigen+","O-antigen-")))
ft_full <- fisher.test(mat_full)

fisher_full_result <- data.frame(
  Test             = "Fisher exact (full dataset)",
  Hist_OBC_plus    = hist_plus,
  Hist_OBC_minus   = hist_minus,
  Hist_total       = hist_total,
  Mod_OBC_plus     = mod_plus,
  Mod_OBC_minus    = mod_minus,
  Mod_total        = mod_total,
  Odds_ratio       = round(unname(ft_full$estimate), 4),
  P_value          = ft_full$p.value,
  Significance     = sig_stars(ft_full$p.value),
  stringsAsFactors = FALSE
)

# Write to middle panel folder (existing format)
write.table(
  data.frame(
    Test = "Fisher_exact",
    Historical_Oantigen_plus  = hist_plus,
    Historical_Oantigen_minus = hist_minus,
    Modern_Oantigen_plus      = mod_plus,
    Modern_Oantigen_minus     = mod_minus,
    Odds_ratio = unname(ft_full$estimate),
    P_value    = ft_full$p.value
  ),
  file.path(outdir_middle, "Oantigen_freq_hist_vs_modern1350_separated_fisher_test.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)

# ---- E6. FISHER EXACT TEST: cross-era with repeated random downsampling ----
#
# Rationale: The sample size imbalance (Historical n=38, Modern n~1312) means
# Fisher's exact test on the full dataset is dominated by the Modern group.
# To evaluate robustness, we apply two approaches:
#
#   Approach A – Expected-value (deterministic): compute what Modern O-antigen+/-
#   counts would be if Modern had exactly n_hist=38 samples, using the observed
#   Modern proportion. Run Fisher's exact test on these expected counts.
#   This matches the approach used for country-level comparisons in Fig4C lower panel.
#
#   Approach B – Repeated random downsampling (stochastic, 1,000 draws):
#   Randomly draw n_hist=38 isolates WITHOUT REPLACEMENT from the Modern pool
#   and run Fisher's exact test on each draw. Repeat 1,000 times.
#   The goal is to assess ROBUSTNESS: does the Fisher test consistently yield
#   p < 0.05 regardless of which 38 Modern isolates are included?
#   If Prop_sig_05 ≈ 1.0, the signal is robust to the specific Modern sample.
#
#   TERMINOLOGY NOTE: This is called "repeated random downsampling" (or
#   "resampling robustness analysis"), NOT a permutation test (which randomises
#   group labels to test a null hypothesis) and NOT a bootstrap (which resamples
#   WITH replacement to estimate a statistic's sampling distribution).
#
n_hist    <- hist_total   # 38
mod_prop_plus <- mod_plus / mod_total

# Approach A: expected counts
exp_mod_plus  <- round(n_hist * mod_prop_plus)
exp_mod_minus <- n_hist - exp_mod_plus

mat_exp <- matrix(c(hist_plus, hist_minus, exp_mod_plus, exp_mod_minus),
                  nrow = 2, byrow = TRUE,
                  dimnames = list(c("Historical","Modern_expected_n38"),
                                  c("O-antigen+","O-antigen-")))
ft_exp <- fisher.test(mat_exp)

# Approach B: repeated random downsampling (without replacement each draw)
set.seed(42)
n_downsamp   <- 1000
modern_pool  <- c(rep("O-antigen+", mod_plus), rep("O-antigen-", mod_minus))

downsamp_pvals <- numeric(n_downsamp)
for (i in seq_len(n_downsamp)) {
  samp    <- sample(modern_pool, n_hist, replace = FALSE)  # random subsample, no replacement
  s_plus  <- sum(samp == "O-antigen+")
  s_minus <- n_hist - s_plus
  mat_b   <- matrix(c(hist_plus, hist_minus, s_plus, s_minus), nrow = 2, byrow = TRUE)
  downsamp_pvals[i] <- fisher.test(mat_b)$p.value
}

# Compile downsampled results table
fisher_downsamp <- data.frame(
  # Shared inputs
  Hist_OBC_plus             = hist_plus,
  Hist_OBC_minus            = hist_minus,
  Hist_total                = hist_total,
  Mod_OBC_plus_full         = mod_plus,
  Mod_OBC_minus_full        = mod_minus,
  Mod_total_full            = mod_total,
  Downsample_n              = n_hist,
  # Approach A
  ApproachA_Test            = "Fisher exact (expected counts, Mod downsampled to n=38)",
  ApproachA_Mod_exp_plus    = exp_mod_plus,
  ApproachA_Mod_exp_minus   = exp_mod_minus,
  ApproachA_Odds_ratio      = round(unname(ft_exp$estimate), 4),
  ApproachA_P_value         = ft_exp$p.value,
  ApproachA_Significance    = sig_stars(ft_exp$p.value),
  # Approach B: repeated random downsampling
  ApproachB_Test            = "Fisher exact (repeated random downsampling, Mod n=38 w/o replacement, 1,000 draws)",
  ApproachB_N_draws         = n_downsamp,
  ApproachB_Median_P        = round(median(downsamp_pvals), 6),
  ApproachB_Mean_P          = round(mean(downsamp_pvals), 6),
  ApproachB_Prop_sig_05     = round(mean(downsamp_pvals < 0.05), 4),
  ApproachB_Significance    = sig_stars(median(downsamp_pvals)),
  stringsAsFactors          = FALSE
)

# ---- E6b. CHI-SQUARE GOODNESS-OF-FIT: 4 HTF groups within each era ----
#
# Rationale: The binomial test above addresses only the binary OBC+/- split.
# For comparing the distribution ACROSS ALL FOUR HTF haplotype groups within
# each era (is any group disproportionately represented?), a chi-square goodness
# of fit test is more appropriate than 4 separate binomial tests.
#
# H0: All four groups are equally frequent within each era (p = 0.25 each).
# H1: At least one group departs from uniform distribution.
#
# NOTE on binomial appropriateness: the binomial test (E3) is valid for the
# COLLAPSED binary question (OBC+ vs OBC-). For the 4-group question, use
# chi-square or multinomial. Both are provided here for completeness.
#
chi_tests <- lapply(c("Historical","Modern"), function(era) {
  d <- df_combined4 %>%
    filter(IsolateType == era) %>%
    arrange(LengthGroup)
  # Integer counts per group
  obs <- round(d$Frequency)
  total_n <- unique(d$TotalN)
  # Chi-square goodness of fit vs. uniform (p=0.25 each)
  chisq_uniform <- chisq.test(obs, p = rep(0.25, length(obs)))
  data.frame(
    IsolateType       = era,
    HTF_1803_n        = obs[d$LengthGroup == "1803"],
    HTF_1245_n        = obs[d$LengthGroup == "1245"],
    HTF_1383_n        = obs[d$LengthGroup == "1383"],
    HTF_1830_n        = obs[d$LengthGroup == "1830"],
    Total_N           = total_n,
    Test              = "Chi-square goodness-of-fit (H0: uniform 25% each)",
    ChiSq_statistic   = round(chisq_uniform$statistic, 3),
    df                = chisq_uniform$parameter,
    P_value           = chisq_uniform$p.value,
    Significance      = sig_stars(chisq_uniform$p.value),
    Interpretation    = "H0: each of 4 HTF groups equally likely within era",
    stringsAsFactors  = FALSE
  )
}) %>% bind_rows()

write.table(chi_tests,
            file.path(outdir_upper, "HTF4group_ChiSq_uniform_within_era.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

write.table(fisher_downsamp,
            file.path(outdir_upper,  "HTF4group_Fisher_OBC_crossera_downsampled.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(fisher_downsamp,
            file.path(outdir_middle, "Oantigen_Fisher_OBC_crossera_downsampled.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# =============================================================================
# SECTION F: COMPREHENSIVE SUPPLEMENTARY STATISTICS TABLE
# (Single table combining all tests, ready for PNAS supplementary material)
# =============================================================================
supp_table <- bind_rows(
  # Row 1: Binomial – Historical OBC
  data.frame(
    Panel            = "Upper & Middle",
    Comparison       = "Within Historical: O-antigen+ vs O-antigen- (H0: p=0.5)",
    Test             = "Two-sided exact binomial",
    Group            = "All HTF groups aggregated",
    Historical_n     = paste0(hist_plus, "+, ", hist_minus, "-  (total ", hist_total, ")"),
    Modern_n         = NA_character_,
    Statistic_type   = "p-value",
    Statistic_value  = signif(binomial_tests$P_value[binomial_tests$IsolateType=="Historical"], 4),
    Significance     = binomial_tests$Significance[binomial_tests$IsolateType=="Historical"],
    Notes            = "Testing departure from 50/50 O-antigen split within Historical era",
    stringsAsFactors = FALSE
  ),
  # Row 2: Binomial – Modern OBC
  data.frame(
    Panel            = "Upper & Middle",
    Comparison       = "Within Modern: O-antigen+ vs O-antigen- (H0: p=0.5)",
    Test             = "Two-sided exact binomial",
    Group            = "All HTF groups aggregated",
    Historical_n     = NA_character_,
    Modern_n         = paste0(mod_plus, "+, ", mod_minus, "-  (total ", mod_total, ")"),
    Statistic_type   = "p-value",
    Statistic_value  = signif(binomial_tests$P_value[binomial_tests$IsolateType=="Modern"], 4),
    Significance     = binomial_tests$Significance[binomial_tests$IsolateType=="Modern"],
    Notes            = "Testing departure from 50/50 O-antigen split within Modern era",
    stringsAsFactors = FALSE
  ),
  # Rows 3-6: Fisher per group
  fisher_pergroup %>%
    mutate(
      Panel           = "Upper",
      Comparison      = paste0("Cross-era: Historical vs Modern (HTF group ", LengthGroup, ")"),
      Test            = "Fisher's exact test",
      Group           = paste0("HTF-", LengthGroup, " (", Oantigen, ")"),
      Historical_n    = paste0(Hist_count, "/", Hist_total, " (", Hist_pct, "%)"),
      Modern_n        = paste0(Mod_count,  "/", Mod_total,  " (", Mod_pct,  "%)"),
      Statistic_type  = "OR; p-value",
      Statistic_value = signif(P_value, 4),
      Notes           = paste0("Odds ratio = ", Odds_ratio)
    ) %>%
    select(Panel, Comparison, Test, Group, Historical_n, Modern_n,
           Statistic_type, Statistic_value, Significance, Notes),
  # Row 7: Fisher overall full dataset
  data.frame(
    Panel            = "Middle",
    Comparison       = "Cross-era: Historical vs Modern, O-antigen+/- (full dataset)",
    Test             = "Fisher's exact test",
    Group            = "O-antigen+/- totals",
    Historical_n     = paste0(hist_plus, "/", hist_total, " O-antigen+ (", round(hist_plus/hist_total*100,1), "%)"),
    Modern_n         = paste0(mod_plus,  "/", mod_total,  " O-antigen+ (", round(mod_plus/mod_total*100,1),   "%)"),
    Statistic_type   = "OR; p-value",
    Statistic_value  = signif(ft_full$p.value, 4),
    Significance     = sig_stars(ft_full$p.value),
    Notes            = paste0("Odds ratio = ", round(unname(ft_full$estimate),4),
                              "; full cohort, unadjusted for n imbalance"),
    stringsAsFactors = FALSE
  ),
  # Row 8: Fisher downsampled Approach A
  data.frame(
    Panel            = "Middle",
    Comparison       = "Cross-era: Historical vs Modern (Modern expected n=38)",
    Test             = "Fisher's exact test (expected-value downsampling)",
    Group            = "O-antigen+/- totals",
    Historical_n     = paste0(hist_plus, "/", hist_total, " O-antigen+"),
    Modern_n         = paste0(exp_mod_plus, "/", n_hist, " O-antigen+ (expected at Modern proportion)"),
    Statistic_type   = "OR; p-value",
    Statistic_value  = signif(ft_exp$p.value, 4),
    Significance     = sig_stars(ft_exp$p.value),
    Notes            = paste0("Modern downsampled to n=38 using observed proportion ",
                              round(mod_prop_plus*100,1),"%; OR=",
                              round(unname(ft_exp$estimate),4)),
    stringsAsFactors = FALSE
  ),
  # Row 9: Fisher with repeated random downsampling (Approach B)
  data.frame(
    Panel            = "Middle",
    Comparison       = "Cross-era: Historical vs Modern (repeated random downsampling, Mod n=38 w/o replacement, 1,000 draws)",
    Test             = "Fisher's exact test (repeated random downsampling robustness)",
    Group            = "O-antigen+/- totals",
    Historical_n     = paste0(hist_plus, "/", hist_total, " O-antigen+"),
    Modern_n         = paste0(mod_plus, "/", mod_total, " Modern pool (randomly subsampled to n=38 per draw)"),
    Statistic_type   = "Median p-value (10k random draws)",
    Statistic_value  = round(median(downsamp_pvals), 6),
    Significance     = sig_stars(median(downsamp_pvals)),
    Notes            = paste0("Prop. draws p<0.05 = ",
                              round(mean(downsamp_pvals < 0.05)*100, 1),
                              "%; repeated downsampling robustness (not permutation, not bootstrap)"),
    stringsAsFactors = FALSE
  ),
  # Row 10: Chi-square within Historical
  data.frame(
    Panel            = "Upper",
    Comparison       = "Within Historical: 4 HTF groups vs uniform distribution (H0: 25% each)",
    Test             = "Chi-square goodness-of-fit",
    Group            = "4 HTF groups (1803, 1245, 1383, 1830)",
    Historical_n     = paste0(chi_tests$Total_N[chi_tests$IsolateType=="Historical"], " total"),
    Modern_n         = NA_character_,
    Statistic_type   = "Chi-sq statistic; p-value",
    Statistic_value  = signif(chi_tests$P_value[chi_tests$IsolateType=="Historical"], 4),
    Significance     = chi_tests$Significance[chi_tests$IsolateType=="Historical"],
    Notes            = paste0("Chi-sq = ", chi_tests$ChiSq_statistic[chi_tests$IsolateType=="Historical"],
                              ", df=", chi_tests$df[chi_tests$IsolateType=="Historical"],
                              "; use for 4-group comparison, not binomial"),
    stringsAsFactors = FALSE
  ),
  # Row 11: Chi-square within Modern
  data.frame(
    Panel            = "Upper",
    Comparison       = "Within Modern: 4 HTF groups vs uniform distribution (H0: 25% each)",
    Test             = "Chi-square goodness-of-fit",
    Group            = "4 HTF groups (1803, 1245, 1383, 1830)",
    Historical_n     = NA_character_,
    Modern_n         = paste0(chi_tests$Total_N[chi_tests$IsolateType=="Modern"], " total"),
    Statistic_type   = "Chi-sq statistic; p-value",
    Statistic_value  = signif(chi_tests$P_value[chi_tests$IsolateType=="Modern"], 4),
    Significance     = chi_tests$Significance[chi_tests$IsolateType=="Modern"],
    Notes            = paste0("Chi-sq = ", chi_tests$ChiSq_statistic[chi_tests$IsolateType=="Modern"],
                              ", df=", chi_tests$df[chi_tests$IsolateType=="Modern"],
                              "; use for 4-group comparison, not binomial"),
    stringsAsFactors = FALSE
  )
)

# NOTE: Two-version combined supplementary tables (FULL + forpaper) are
# assembled by step5_fig4c_lower_countryboxplot_bottom.R, which reads the
# individual table files written above and appends the lower panel (Germany n=8) rows.
# Run step4 first, then step5 to produce:
#   Supplementary_Statistics_Fig4C_FULL.tsv      (all panels, all tests)
#   Supplementary_Statistics_Fig4C_forpaper.tsv  (upper within, middle within, lower key tests)

# =============================================================================
# SECTION G: HUMAN-READABLE SUMMARY TEXT FILES
# =============================================================================

# ---- G1. Upper panel summary ----
zz <- file(file.path(outdir_upper, "HTF4group_summary.txt"), open = "wt")
writeLines("=== HTF 4-group frequency summary: Historical vs Modern ===", zz)
writeLines(paste0("Date: ", Sys.Date()), zz)
writeLines("", zz)

writeLines("--- Per-group counts used in the plot ---", zz)
capture.output(print(as.data.frame(upper_counts_long), row.names=FALSE), file=zz)
writeLines("", zz)

writeLines("--- O-antigen totals aggregated from the 4 groups ---", zz)
writeLines(sprintf("Historical: O-antigen+ = %d/%d (%.1f%%), O-antigen- = %d/%d (%.1f%%)",
                   hist_plus, hist_total, hist_plus/hist_total*100,
                   hist_minus, hist_total, hist_minus/hist_total*100), zz)
writeLines(sprintf("Modern:     O-antigen+ = %d/%d (%.1f%%), O-antigen- = %d/%d (%.1f%%)",
                   mod_plus, mod_total, mod_plus/mod_total*100,
                   mod_minus, mod_total, mod_minus/mod_total*100), zz)
writeLines("", zz)

writeLines("--- Statistical tests ---", zz)
writeLines("", zz)

writeLines("1. Binomial test within Historical (H0: OBC+ proportion = 0.5):", zz)
bt_h <- binomial_tests[binomial_tests$IsolateType=="Historical",]
writeLines(sprintf("   %d OBC+ vs %d OBC- (total n=%d): p = %.6f [%s]",
                   bt_h$OBC_plus_n, bt_h$OBC_minus_n, bt_h$Total_N,
                   bt_h$P_value, bt_h$Significance), zz)
writeLines("", zz)

writeLines("2. Binomial test within Modern (H0: OBC+ proportion = 0.5):", zz)
bt_m <- binomial_tests[binomial_tests$IsolateType=="Modern",]
writeLines(sprintf("   %d OBC+ vs %d OBC- (total n=%d): p = %.2e [%s]",
                   bt_m$OBC_plus_n, bt_m$OBC_minus_n, bt_m$Total_N,
                   bt_m$P_value, bt_m$Significance), zz)
writeLines("", zz)

writeLines("3. Fisher exact test per HTF group (Historical vs Modern):", zz)
for (i in seq_len(nrow(fisher_pergroup))) {
  r <- fisher_pergroup[i,]
  writeLines(sprintf("   HTF-%s (%s): Hist %d/%d (%.1f%%) vs Mod %d/%d (%.1f%%) | OR=%.3f | p=%.4e [%s]",
                     r$LengthGroup, r$Oantigen,
                     r$Hist_count, r$Hist_total, r$Hist_pct,
                     r$Mod_count,  r$Mod_total,  r$Mod_pct,
                     r$Odds_ratio, r$P_value, r$Significance), zz)
}
writeLines("", zz)

writeLines("4. Fisher exact test: O-antigen overall, Historical vs Modern (full dataset):", zz)
writeLines(sprintf("   Hist %d/%d vs Mod %d/%d | OR=%.4f | p=%.4e [%s]",
                   hist_plus, hist_total, mod_plus, mod_total,
                   round(unname(ft_full$estimate),4), ft_full$p.value,
                   sig_stars(ft_full$p.value)), zz)
writeLines("", zz)

writeLines("5. Fisher exact test: O-antigen, Modern downsampled to n=38 (expected-value approach A):", zz)
writeLines(sprintf("   Hist %d/%d vs Mod_exp %d/%d | OR=%.4f | p=%.4e [%s]",
                   hist_plus, n_hist, exp_mod_plus, n_hist,
                   round(unname(ft_exp$estimate),4), ft_exp$p.value,
                   sig_stars(ft_exp$p.value)), zz)
writeLines("", zz)

writeLines("6. Fisher exact test: O-antigen, repeated random downsampling (1,000 draws, n=38 w/o replacement):", zz)
writeLines(sprintf("   Median p = %.6f | Mean p = %.6f | Prop. draws sig (p<0.05) = %.1f%% [%s]",
                   median(downsamp_pvals), mean(downsamp_pvals),
                   mean(downsamp_pvals<0.05)*100,
                   sig_stars(median(downsamp_pvals))), zz)
writeLines("", zz)
writeLines("7. Chi-square goodness-of-fit: 4 HTF groups within each era (H0: uniform 25% each):", zz)
for (era in c("Historical","Modern")) {
  r <- chi_tests[chi_tests$IsolateType==era,]
  writeLines(sprintf("   %s: Chi-sq=%.3f, df=%d, p=%.4e [%s]",
                     era, r$ChiSq_statistic, r$df, r$P_value, r$Significance), zz)
}
close(zz)

# ---- G2. Middle panel summary (same format as existing, extended) ----
zz2 <- file(file.path(outdir_middle, "Oantigen_freq_hist_vs_modern1350_separated_summary.txt"), open = "wt")
writeLines("=== O-antigen frequency summary: Historical vs Modern (OTU5 strains) ===", zz2)
writeLines(paste0("Date: ", Sys.Date()), zz2)
writeLines("", zz2)

writeLines("Counts used in the plot:", zz2)
capture.output(print(as.data.frame(oantigen_count_table), row.names=FALSE), file=zz2)
writeLines("", zz2)

writeLines("2x2 contingency table used for Fisher's exact test:", zz2)
capture.output(print(mat_full), file=zz2)
writeLines("", zz2)

writeLines(paste0("Fisher's exact test p-value (full dataset): ", signif(ft_full$p.value, 6)), zz2)
writeLines(paste0("Odds ratio: ", round(unname(ft_full$estimate), 4)), zz2)
writeLines("", zz2)

writeLines("Binomial test within Historical (H0: OBC+ = 0.5):", zz2)
writeLines(sprintf("  %d/%d OBC+ | p = %.6f [%s]",
                   hist_plus, hist_total,
                   binomial_tests$P_value[binomial_tests$IsolateType=="Historical"],
                   binomial_tests$Significance[binomial_tests$IsolateType=="Historical"]), zz2)
writeLines("", zz2)

writeLines("Binomial test within Modern (H0: OBC+ = 0.5):", zz2)
writeLines(sprintf("  %d/%d OBC+ | p = %.2e [%s]",
                   mod_plus, mod_total,
                   binomial_tests$P_value[binomial_tests$IsolateType=="Modern"],
                   binomial_tests$Significance[binomial_tests$IsolateType=="Modern"]), zz2)
writeLines("", zz2)

writeLines("Fisher's exact test (Modern downsampled to n=38, expected-value approach):", zz2)
writeLines(sprintf("  Hist %d/%d vs Mod_exp %d/%d | p = %.4e [%s]",
                   hist_plus, n_hist, exp_mod_plus, n_hist,
                   ft_exp$p.value, sig_stars(ft_exp$p.value)), zz2)
writeLines("", zz2)

writeLines("Fisher's exact test (repeated random downsampling, 1,000 draws, n=38 w/o replacement):", zz2)
writeLines(sprintf("  Median p = %.6f | Prop. draws sig (p<0.05) = %.1f%% [%s]",
                   median(downsamp_pvals), mean(downsamp_pvals<0.05)*100,
                   sig_stars(median(downsamp_pvals))), zz2)
writeLines("", zz2)

writeLines("Summary:", zz2)
writeLines(sprintf("Historical: O-antigen+ = %d/%d (%.1f%%), O-antigen- = %d/%d (%.1f%%).",
                   hist_plus, hist_total, hist_plus/hist_total*100,
                   hist_minus, hist_total, hist_minus/hist_total*100), zz2)
writeLines(sprintf("Modern:     O-antigen+ = %d/%d (%.1f%%), O-antigen- = %d/%d (%.1f%%).",
                   mod_plus, mod_total, mod_plus/mod_total*100,
                   mod_minus, mod_total, mod_minus/mod_total*100), zz2)
close(zz2)

# =============================================================================
cat("\n✅ All outputs written.\n")
cat("   Upper panel folder:", outdir_upper, "\n")
cat("   Middle panel folder:", outdir_middle, "\n")
cat("   Supplementary tables (FULL + forpaper): written by step5 after running both scripts\n")
cat("\nKey results:\n")
cat(sprintf("  Binomial (Historical OBC+/- = 50/50): p = %.4f [%s]\n",
            binomial_tests$P_value[binomial_tests$IsolateType=="Historical"],
            binomial_tests$Significance[binomial_tests$IsolateType=="Historical"]))
cat(sprintf("  Binomial (Modern OBC+/- = 50/50):     p = %.2e [%s]\n",
            binomial_tests$P_value[binomial_tests$IsolateType=="Modern"],
            binomial_tests$Significance[binomial_tests$IsolateType=="Modern"]))
cat(sprintf("  Fisher (cross-era, full):    p = %.2e [%s]\n",
            ft_full$p.value, sig_stars(ft_full$p.value)))
cat(sprintf("  Fisher (cross-era, n38 exp): p = %.2e [%s]\n",
            ft_exp$p.value, sig_stars(ft_exp$p.value)))
cat(sprintf("  Fisher (cross-era, downsampling robustness median p): %.4e [%s]\n",
            median(downsamp_pvals), sig_stars(median(downsamp_pvals))))

# =============================================================================
# SECTION H: WRITE TO final_comparison FOLDER
# Upper panel (4HTF binomial) and middle panel (OBC binomial) tables are
# contributed here. Lower panel (Germany n=8 Fisher + permutation) is
# contributed by step5_fig4c_lower_countryboxplot_bottom.R.
# =============================================================================
outdir_final <- file.path(base_outdir, "final_comparison")
dir.create(outdir_final, showWarnings = FALSE, recursive = TRUE)

# Table 1 of 4: Upper panel — 4 HTF groups, binomial within each era
write.table(binomial_tests,
            file.path(outdir_final, "Table1_Fig4C_upper_4HTFgroups_OBC_binomial_within_era.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# Table 2 of 4: Middle panel — OBC binary, binomial within each era
write.table(binomial_tests,   # same underlying test; identical data
            file.path(outdir_final, "Table2_Fig4C_middle_OBC_binomial_within_era.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# Also copy the chi-square 4-group test as an additional upper-panel table
write.table(chi_tests,
            file.path(outdir_final, "Table1b_Fig4C_upper_4HTFgroups_chisq_uniform_within_era.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

cat("  final_comparison tables 1 & 2 written to:", outdir_final, "\n")
cat("  (Tables 3 & 4 — lower panel Germany n=8 — are written by step5.R)\n")
