# ============================
# Supp Fig 1: Endogenous composition + Coverage/Depth
# ============================

# ---- Libraries ----
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(ggtext)   # for markdown/italics in axis/legend text

# ---- Paths ----
setwd("/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step1_h40metadata/visualization_suppfig1")

infile <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step1_h40metadata/allh46_metadata.txt"

# ---- Load ----
df <- read.table(infile, header = TRUE, sep = "\t")

# ---- Colors ----
col_at    <- "#4575b4"   # academic blue
  col_ps    <- "#d73027"   # academic red
    col_other <- "grey80"    # others grey
      
    # ============================
    # (A) Endogenous composition stacked barplot
    # ============================
    df_long <- df %>%
      mutate(Other_percent = 1 - (At_percent + Ps_percent)) %>%
      arrange(At_percent) %>%   # order by At proportion
      mutate(samplename = factor(samplename, levels = rev(samplename))) %>%
      select(samplename, At_percent, Ps_percent, Other_percent) %>%
      pivot_longer(cols = c("At_percent","Ps_percent","Other_percent"),
                   names_to = "Species", values_to = "Prop") %>%
      mutate(Species = factor(Species,
                              levels = c("Other_percent","Ps_percent","At_percent"),
                              labels = c("Other",
                                         "*Pseudomonas* sp.",
                                         "*Arabidopsis thaliana*")))
    
    p_a <- ggplot(df_long, aes(x = samplename, y = Prop, fill = Species)) +
      geom_bar(stat = "identity", width = 0.8) +
      scale_fill_manual(values = c("Other" = col_other,
                                   "*Pseudomonas* sp." = col_ps,
                                   "*Arabidopsis thaliana*" = col_at)) +
      labs(x = "Sample", y = "Proportion of endogenous reads", fill = "") +
      theme_classic(base_size = 12) +
      theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5),
            axis.ticks.x = element_blank(),
            legend.position = "top",
            legend.title = element_blank(),
            legend.text = element_markdown(size = 11, face = "bold"),  # italics in legend
            axis.title.x = element_text(size = 11, face = "bold"),  # <-- Sample not bold
            axis.title.y = element_text(size = 11, face = "bold"),
            axis.text.y = element_text(size = 11, face = "bold"))
    
    # ============================
    # (B) Covered proportion vs depth scatter (pseudo broken x-axis)
    # ============================
    df_scatter <- df %>%
      select(samplename, At_cov = At_covered, At_depth,
             Ps_cov = Ps_covered, Ps_depth) %>%
      pivot_longer(cols = everything(),
                   names_to = c("Species",".value"),
                   names_pattern = "(.*)_(.*)") %>%
      mutate(Species = recode(Species,
                              "At"="*Arabidopsis thaliana*",
                              "Ps"="*Pseudomonas* sp.")) %>%
      filter(!is.na(depth) & !is.na(cov))
    
    # squash >200 depth to 220
    df_scatter <- df_scatter %>%
      mutate(depth_plot = ifelse(depth > 200, 220, depth))
    
    p_b <- ggplot(df_scatter, aes(x = depth_plot, y = cov, color = Species)) +
      geom_point(size = 3, alpha = 0.9) +
      scale_color_manual(values = c("*Arabidopsis thaliana*" = col_at,
                                    "*Pseudomonas* sp." = col_ps),
                         labels = c("*Arabidopsis thaliana*", "*Pseudomonas* sp.")) +
      labs(x = "Mean depth", y = "Covered proportion", color = "") +
      scale_x_continuous(breaks = c(0,50,100,150,220),
                         labels = c("0","50","100","150","500")) +
      # add crossing slashes with segments
      
      coord_cartesian(ylim = c(0,1)) +
      theme_classic(base_size = 12) +
      theme(legend.position = "top",
            legend.title = element_blank(),
            legend.text = element_markdown(size = 11, face = "bold"), # italics in legend
            axis.title = element_text(size = 11, face = "bold"),
            axis.text = element_text(size = 11, face = "bold"))
    
    # ============================
    # Combine panels A + B
    # ============================
    p_final <- p_a + p_b + plot_layout(ncol = 2, widths = c(1.5,1)) +
      plot_annotation(tag_levels = 'A')
    
    # ---- Save ----
    ggsave("SuppFig1_Endogenous_AT_Ps_fixed.pdf",
           plot = p_final, width = 13, height = 5.5, dpi = 600)
    ggsave("SuppFig1a_Endogenous_AT_Ps_fixed.pdf",
           plot = p_a, width = 6, height = 5.5, dpi = 600)
    ggsave("SuppFig1b_Endogenous_AT_Ps_fixed.pdf",
           plot = p_b, width = 6, height = 5.5, dpi = 600)