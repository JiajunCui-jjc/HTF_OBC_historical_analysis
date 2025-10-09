# ============================
# Supp Fig 1: Endogenous composition + Coverage/Depth (final corrected)
# ============================

# ---- Libraries ----
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(ggtext)

# ---- Paths ----
setwd("/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step1_h40metadata/visualization_suppfig1")

infile <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step1_h40metadata/allh46_metadata.txt"
df <- read.table(infile, header = TRUE, sep = "\t")

# ---- Colors ----
col_at <- "#d5da6d"   # blue
  col_ps <- "#e86b7d"   # 
    col_other <- "grey85" # 
      
    # ============================
    # (A) Endogenous composition stacked barplot
    # ============================
    
    df_long <- df %>%
      mutate(Other_percent = 1 - (At_percent + Ps_percent)) %>%
      arrange(Ps_percent) %>%
      mutate(samplename = factor(samplename, levels = rev(samplename))) %>%
      select(samplename, At_percent, Ps_percent, Other_percent) %>%
      pivot_longer(cols = c("Other_percent","At_percent","Ps_percent"),
                   names_to = "Species", values_to = "Prop") %>%
      mutate(Species = factor(Species,
                              levels = c("Other_percent","At_percent","Ps_percent"),
                              labels = c("Other",
                                         "*Arabidopsis thaliana*",
                                         "*Pseudomonas* sp.")))
    
    # The order above ensures stacking: top = Other → At → Ps
    
    p_a <- ggplot(df_long, aes(x = samplename, y = Prop, fill = Species)) +
      geom_bar(stat = "identity", width = 0.8) +
      scale_fill_manual(values = c("Other" = col_other,
                                   "*Arabidopsis thaliana*" = col_at,
                                   "*Pseudomonas* sp." = col_ps)) +
      labs(x = "Sample", y = "Proportion of endogenous reads") +
      theme_classic(base_size = 12) +
      theme(
        axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_markdown(size = 11, face = "bold"),
        axis.title = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold")
      )
    
    # ============================
    # (B) Covered proportion vs mean depth (separate species panels)
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
      filter(!is.na(depth) & !is.na(cov)) %>%
      mutate(depth_plot = ifelse(depth > 200, 220, depth))
    #make the x more closer but the real number is 488
    # --- Panel B1: Arabidopsis ---
    p_b1 <- df_scatter %>%
      filter(Species == "*Arabidopsis thaliana*") %>%
      ggplot(aes(x = depth_plot, y = cov)) +
      geom_point(color = col_at, size = 3, alpha = 0.9) +
      scale_x_continuous(breaks = c(0,5,10,15, 20 ,25),
                         labels = c("0","5","10","15","20","25")) +
      coord_cartesian(xlim = c(0,25), ylim = c(0,1)) +
      labs(x = "Mean depth", y = "Covered proportion",
           title = expression(italic("Arabidopsis thaliana"))) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(size = 13, face = "bold.italic", hjust = 0.5),
        axis.title = element_text(size = 11, face = "bold"),
        axis.text = element_text(size = 11, face = "bold")
      )
    
    # --- Panel B2: Pseudomonas ---
    p_b2 <- df_scatter %>%
      filter(Species == "*Pseudomonas* sp.") %>%
      ggplot(aes(x = depth_plot, y = cov)) +
      geom_point(color = col_ps, size = 3, alpha = 0.9) +
      scale_x_continuous(breaks = c(0,50,100,150,200,232),
                         labels = c("0","50","100","150","200","500")) +
      coord_cartesian(xlim = c(0,232), ylim = c(0,1)) +
      labs(x = "Mean depth", y = "Covered proportion",
           title = expression(italic("Pseudomonas")~"sp.")) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(size = 13, face = "bold.italic", hjust = 0.5),
        axis.title = element_text(size = 11, face = "bold"),
        axis.text = element_text(size = 11, face = "bold")
      )
    
    # ============================
    # Combine all panels
    # ============================
    
    p_final <- (p_a / (p_b1 + p_b2)) +
      plot_annotation(tag_levels = 'A') +
      plot_layout(heights = c(1.3, 1))
    
    # ---- Save ----
    ggsave("SuppFig1_Endogenous_CoverageDepth_final.pdf",
           plot = p_final, width = 12, height = 8, dpi = 600)
    
    ggsave("SuppFig1a_Endogenous_sortedByPs.pdf",
           plot = p_a, width = 6.5, height = 5.5, dpi = 600)
    
    ggsave("SuppFig1b_At_Ps_separate.pdf",
           plot = p_b1 + p_b2 + plot_layout(ncol = 2),
           width = 10, height = 4.5, dpi = 600)