
#any pair with low divergence <10-5 and NA recombine site, was set to recombination prop=0 


suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(scales)
})

setwd('/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step5_recophy/report/recophy_downloads/')

#/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step5_recophy/report/recophy_downloads/mixmodel.tsv
pairwise_file <- "mixmodel.tsv"

# ---- Load ----
df <- read.table(pairwise_file)
colnames(df) <- c("strain1","strain2","divergence","clonal_fraction",
                  "SNP_rate_clonal","SNP_rate_recomb","SNPs_clonal","SNPs_recomb")

# ---- Remove redundant reference, the same as p25.c2 ----
df <- df %>%
  filter(strain1 != "ref" & strain2 != "ref")


# ==========================================================
# Step 1: Identify near-identical (clonal) pairs first
# ==========================================================
df <- df %>%
  mutate(across(c(divergence, SNPs_clonal, SNPs_recomb,
                  SNP_rate_clonal, SNP_rate_recomb),
                ~ suppressWarnings(as.numeric(.))))

nearclonal_df <- df %>%
  filter(!is.na(divergence) & divergence <= 1e-5)

if (nrow(nearclonal_df) > 0) {
  write_tsv(nearclonal_df, "pairs_nearclonal_bydivergence.tsv")
  cat("ℹ️  Found", nrow(nearclonal_df),
      "near-identical pairs (divergence ≤ 1e-5). Saved to pairs_nearclonal_bydivergence.tsv\n")
} else {
  cat("✅  No pairs with divergence ≤ 1e-5 detected.\n")
}

# ==========================================================
# Step 2: Derived values
# ==========================================================
df <- df %>%
  mutate(
    SNPs_clonal = ifelse(is.na(SNPs_clonal), 0, SNPs_clonal),
    SNPs_recomb = ifelse(is.na(SNPs_recomb), 0, SNPs_recomb),
    recomb_fraction = ifelse(
      (SNPs_clonal + SNPs_recomb) == 0, 0,
      SNPs_recomb / (SNPs_clonal + SNPs_recomb)
    ),
    log_ratio_recomb = ifelse(
      SNPs_clonal <= 0 & SNPs_recomb <= 0, 0,
      ifelse(SNPs_clonal <= 0, 0, log10(SNPs_recomb / SNPs_clonal))
    ),
    recomb_fraction = ifelse(divergence <= 1e-5, 0, recomb_fraction),
    log_ratio_recomb = ifelse(divergence <= 1e-5, 0, log_ratio_recomb)
  )

# ==========================================================
# Step 3: Summary statistics
# ==========================================================
summary_stats <- df %>%
  summarise(
    n_pairs = n(),
    mean_clonal = mean(clonal_fraction, na.rm=TRUE),
    sd_clonal = sd(clonal_fraction, na.rm=TRUE),
    median_clonal = median(clonal_fraction, na.rm=TRUE),
    mean_divergence = mean(divergence, na.rm=TRUE),
    mean_recomb_frac = mean(recomb_fraction, na.rm=TRUE),
    sd_recomb_frac = sd(recomb_fraction, na.rm=TRUE),
    median_recomb_frac = median(recomb_fraction, na.rm=TRUE)
  )
write_tsv(summary_stats, "recophy_pairwise_summary.txt")
# ==========================================================
# Step 3.5: Add recombination dominance summary
# ==========================================================

prop_high_recomb <- df %>%
  filter(!is.na(log_ratio_recomb)) %>%
  summarise(
    n_high = sum(log_ratio_recomb > 1, na.rm = TRUE),
    prop_high = mean(log_ratio_recomb > 1, na.rm = TRUE)
  )

write_tsv(prop_high_recomb, "prop_high_recomb_pairs.txt")
summary_text <- paste0(
  "Among ", summary_stats$n_pairs, " ATUE5 genome pairs, the mean clonal fraction was ",
  round(summary_stats$mean_clonal, 3), " ± ", round(summary_stats$sd_clonal, 3),
  " (median = ", round(summary_stats$median_clonal, 3), "). ",
  "On average, ", round(summary_stats$mean_recomb_frac*100, 1), "% ± ",
  round(summary_stats$sd_recomb_frac*100, 1),
  "% of SNPs were introduced by recombination (median = ",
  round(summary_stats$median_recomb_frac*100, 1), "%). ",
  "Notably, ", prop_high_recomb$n_high, " pairs (",
  round(prop_high_recomb$prop_high*100, 1), "%) showed log₁₀(Recombined/Clonal SNPs) > 1, ",
  "indicating recombination generated more than tenfold the number of clonal SNPs. ",
  "This highlights pervasive horizontal gene transfer in ATUE5, yet HTF–O-antigen linkage remains fully conserved."
)
writeLines(summary_text, "summary.words.txt")

# ==========================================================
# Step 4: Plot style
# ==========================================================
main_color <- "#4A90E2"  # uniform blue
  theme_nat <- theme_minimal(base_size = 18) +
    theme(
      panel.grid = element_blank(),  # remove ALL gridlines
      axis.line = element_line(colour = "black", linewidth = 0.6),
      axis.ticks = element_line(colour = "black"),
      axis.text = element_text(color = "black", size = 14),
      axis.title = element_text(size = 16),
      plot.title = element_text(face = "bold", size = 17, hjust = 0.5),
      legend.position = "top",
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 12),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
    )
  
  
  n_total <- nrow(df)
  
  # ==========================================================
  # Panel A – Histogram of clonal fraction
  # ==========================================================
  mean_clonal <- mean(df$clonal_fraction, na.rm=TRUE)
  
  pA <- ggplot(df, aes(x = clonal_fraction)) +
    geom_histogram(aes(y = after_stat(count / n_total)),
                   binwidth = 0.02, fill = main_color, color = "white", linewidth=0.4) +
    geom_vline(aes(xintercept = mean_clonal),
               color = "black", linetype="dashed", linewidth=0.8) +
    scale_y_continuous(labels = percent_format(accuracy=1)) +
    labs(x = "Clonal fraction",
         y = "Proportion of pairwise comparisons",
         title = "A. Recombination landscape",
         caption = "Dashed line = mean clonal fraction") +
    theme_nat
  
  ggsave("recophy_A_clonal_fraction_hist.pdf", pA, width = 6.5, height = 5, useDingbats=FALSE)
  
  # ==========================================================
  # Panel B – Divergence vs clonal fraction (with R² & p-value, no grid)
  # ==========================================================
  
  # ---- Fit linear regression ----
  fit <- lm(clonal_fraction ~ divergence, data = df)
  fit_sum <- summary(fit)
  r2_val <- fit_sum$r.squared
  p_val  <- fit_sum$coefficients[2, 4]
  
  # ---- Format statistics ----
  r2_text <- paste0("R² = ", sprintf("%.3f", r2_val))
  p_text  <- ifelse(p_val < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", p_val)))
  stats_label <- paste(r2_text, p_text, sep = ",  ")
  
  # ---- Label position (top-right corner) ----
  x_pos <- max(df$divergence, na.rm = TRUE) * 0.98
  y_pos <- max(df$clonal_fraction, na.rm = TRUE) * 0.98
  
  # ---- Plot ----
  pB <- ggplot(df, aes(x = divergence, y = clonal_fraction)) +
    geom_point(color = main_color, size = 2.5, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE,
                color = "black", linetype = "dashed", linewidth = 0.8) +
    annotate("text",
             x = x_pos, y = y_pos,
             label = stats_label,
             hjust = 1, vjust = 1,
             size = 5, fontface = "italic", color = "black") +
    labs(
      x = "Pairwise divergence",
      y = "Clonal fraction",
      title = "B. Divergence vs recombination",
      caption = "Dashed line = linear regression fit"
    ) +
    theme_minimal(base_size = 18) +
    theme(
      panel.grid = element_blank(),      # ← remove all gridlines
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black", size = 14),
      axis.title = element_text(size = 16),
      plot.title = element_text(face = "bold", size = 17, hjust = 0.5),
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
    )
  
  ggsave("recophy_B_divergence_vs_clonal_fraction.pdf",
         pB, width = 6.5, height = 5, useDingbats = FALSE)
  
  # ==========================================================
  # Panel C – Histogram of log10(Recombined/Clonal SNPs)
  # ==========================================================
  mean_logratio <- mean(df$log_ratio_recomb, na.rm=TRUE)
  pC <- ggplot(df, aes(x = log_ratio_recomb)) +
    geom_histogram(aes(y = after_stat(count / n_total)),
                   binwidth = 0.1, fill = main_color, color = "white", linewidth=0.4) +
    geom_vline(aes(xintercept = mean_logratio),
               color = "black", linetype = "dashed", linewidth = 0.8) +
    scale_y_continuous(labels = percent_format(accuracy=1)) +
    labs(x = expression(Log[10]~"(Recombined / Clonal SNPs)"),
         y = "Proportion of pairwise comparisons",
         title = "C. Ratio of recombined to clonal SNPs",
         caption = "Dashed line = mean log10 ratio") +
    theme_nat
  
  ggsave("recophy_C_logratio_recomb_vs_clonal.pdf", pC, width = 6.5, height = 5, useDingbats=FALSE)
  
  # ==========================================================
  # Console summary
  # ==========================================================
  cat("\nSummary for manuscript:\n", summary_text, "\n")
  cat("\nSaved files:\n",
      " - recophy_A_clonal_fraction_hist.pdf\n",
      " - recophy_B_divergence_vs_clonal_fraction.pdf\n",
      " - recophy_C_logratio_recomb_vs_clonal.pdf\n",
      " - summary.words.txt\n",
      " - recophy_pairwise_summary.txt\n\n")
  
  
  
  
  # ==========================================================
  # Panel B2 – Divergence vs clonal fraction colored by O-antigen status
  # ==========================================================
  
  # ---- Load O-antigen metadata ----
  
  meta_file <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step3_combine/combined_HTF_oantigen_m57_h43.txt"
  
  meta <- read_tsv(meta_file, show_col_types = FALSE) %>%
    select(sample, Oantigen = HTFgroup_Oantigen_PA) %>%
    mutate(
      # For samples starting with lowercase 'p', replace '.' with '_'
      sample = ifelse(grepl("^p", sample),
                      gsub("\\.", "_", sample),
                      sample)
    )
  
  # ---- Attach O-antigen info to both strains in each pair ----
  df_col <- df %>%
    left_join(meta, by = c("strain1" = "sample")) %>%
    rename(Oantigen1 = Oantigen) %>%
    left_join(meta, by = c("strain2" = "sample")) %>%
    rename(Oantigen2 = Oantigen) %>%
    mutate(
      same_Oantigen = case_when(
        is.na(Oantigen1) | is.na(Oantigen2) ~ NA_character_,
        Oantigen1 == Oantigen2 ~ "Same O-antigen P/A",
        TRUE ~ "Different O-antigen P/A"
      )
    )
  
  # ---- Fit regression again (same model) ----
  fit2 <- lm(clonal_fraction ~ divergence, data = df_col)
  fit2_sum <- summary(fit2)
  r2_val2 <- fit2_sum$r.squared
  p_val2  <- fit2_sum$coefficients[2, 4]
  r2_text2 <- paste0("R² = ", sprintf("%.3f", r2_val2))
  p_text2  <- ifelse(p_val2 < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", p_val2)))
  stats_label2 <- paste(r2_text2, p_text2, sep = ",  ")
  
  x_pos2 <- max(df_col$divergence, na.rm = TRUE) * 0.98
  y_pos2 <- max(df_col$clonal_fraction, na.rm = TRUE) * 0.98
  
  # ---- Color palette ----
  pal <- c("Same O-antigen P/A" = "#4A90E2", "Different O-antigen P/A" = "#E74C3C")
  #c("O-antigen-" = "#9e9ac8",
  #  "O-antigen+" = "#c94f7c")) +
  
  # ---- Plot ----
  pB2 <- ggplot(df_col, aes(x = divergence, y = clonal_fraction, color = same_Oantigen)) +
    geom_point(size = 2.5, alpha = 0.9) +
    geom_smooth(method = "lm", se = TRUE, color = "black",
                linetype = "dashed", linewidth = 0.8) +
    scale_color_manual(values = pal, na.value = "grey70") +
    annotate("text",
             x = x_pos2, y = y_pos2,
             label = stats_label2,
             hjust = 1, vjust = 1,
             size = 4, fontface = "italic", color = "black") +
    labs(
      x = "Pairwise divergence",
      y = "Clonal fraction",
      color = "O-antigen status",
      title = "B2. Divergence vs recombination by O-antigen",
      caption = "Dashed line = linear regression fit"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black", size = 13),
      axis.title = element_text(size = 13),
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      legend.position = "top",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
    )
  
  ggsave("recophy_B2_divergence_vs_clonal_fraction_byOantigen.pdf",
         pB2, width = 7, height = 6, useDingbats = FALSE)
  
  
  
  
  #option 1 transparent:
  # ==========================================================
  # Panel A2 – Semi-transparent overlap
  # ==========================================================
  mean_clonal_all <- mean(df_col$clonal_fraction, na.rm = TRUE)
  
  pA2 <- ggplot(df_col, aes(x = clonal_fraction, fill = same_Oantigen)) +
    geom_histogram(
      aes(y = after_stat(count / n_total)),
      position = "identity",
      binwidth = 0.02, color = "white", linewidth = 0.4, alpha = 0.55   # semi-transparent
    ) +
    geom_vline(aes(xintercept = mean_clonal_all),
               color = "black", linetype = "dashed", linewidth = 0.8) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_fill_manual(values = pal, na.value = "grey80") +
    labs(
      x = "Clonal fraction",
      y = "Proportion of pairwise comparisons",
      fill = "O-antigen status",
      title = "A2. Recombination landscape by O-antigen",
      caption = "Dashed line = mean clonal fraction"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black", size = 13),
      axis.title = element_text(size = 13),
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      legend.position = "top",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
    )
  
  ggsave("recophy_A2_clonal_fraction_hist_byOantigen_transparent.pdf",
         pA2, width = 7, height = 6, useDingbats = FALSE)
  # ==========================================================
  # Panel C2 – Semi-transparent overlap
  # ==========================================================
  mean_logratio_all <- mean(df_col$log_ratio_recomb, na.rm = TRUE)
  
  pC2 <- ggplot(df_col, aes(x = log_ratio_recomb, fill = same_Oantigen)) +
    geom_histogram(
      aes(y = after_stat(count / n_total)),
      position = "identity",
      binwidth = 0.1, color = "white", linewidth = 0.4, alpha = 0.55
    ) +
    geom_vline(aes(xintercept = mean_logratio_all),
               color = "black", linetype = "dashed", linewidth = 0.8) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_fill_manual(values = pal, na.value = "grey80") +
    labs(
      x = expression(Log[10]~"(Recombined / Clonal SNPs)"),
      y = "Proportion of pairwise comparisons",
      fill = "O-antigen status",
      title = "C2. Ratio of recombined to clonal SNPs by O-antigen",
      caption = "Dashed line = mean log10 ratio"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black", size = 13),
      axis.title = element_text(size = 13),
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      legend.position = "top",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
    )
  
  ggsave("recophy_C2_logratio_recomb_vs_clonal_byOantigen_transparent.pdf",
         pC2, width = 7, height = 6, useDingbats = FALSE)
  
  
  #option2 density curve:
  # ==========================================================
  # Panel A2 – Density curves
  # ==========================================================
  pA2_density <- ggplot(df_col, aes(x = clonal_fraction,
                                    color = same_Oantigen, fill = same_Oantigen)) +
    geom_density(alpha = 0.35, linewidth = 1.2) +
    geom_vline(aes(xintercept = mean_clonal_all),
               color = "black", linetype = "dashed", linewidth = 0.8) +
    scale_color_manual(values = pal, na.value = "grey70") +
    scale_fill_manual(values = pal, na.value = "grey70") +
    labs(
      x = "Clonal fraction",
      y = "Density",
      color = "O-antigen status",
      fill  = "O-antigen status",
      title = "A2. Clonal fraction distribution by O-antigen",
      caption = "Dashed line = mean clonal fraction"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black", size = 13),
      axis.title = element_text(size = 13),
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      legend.position = "top",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
    )
  
  ggsave("recophy_A2_clonal_fraction_density_byOantigen.pdf",
         pA2_density, width = 7, height = 6, useDingbats = FALSE)
  # ==========================================================
  # Panel C2 – Density curves
  # ==========================================================
  pC2_density <- ggplot(df_col, aes(x = log_ratio_recomb,
                                    color = same_Oantigen, fill = same_Oantigen)) +
    geom_density(alpha = 0.35, linewidth = 1.2) +
    geom_vline(aes(xintercept = mean_logratio_all),
               color = "black", linetype = "dashed", linewidth = 0.8) +
    scale_color_manual(values = pal, na.value = "grey70") +
    scale_fill_manual(values = pal, na.value = "grey70") +
    labs(
      x = expression(Log[10]~"(Recombined / Clonal SNPs)"),
      y = "Density",
      color = "O-antigen status",
      fill  = "O-antigen status",
      title = "C2. Recombined vs clonal SNPs by O-antigen (density)",
      caption = "Dashed line = mean log10 ratio"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black", size = 13),
      axis.title = element_text(size = 13),
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      legend.position = "top",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
    )
  
  ggsave("recophy_C2_logratio_recomb_vs_clonal_density_byOantigen.pdf",
         pC2_density, width = 7, height = 6, useDingbats = FALSE)
  