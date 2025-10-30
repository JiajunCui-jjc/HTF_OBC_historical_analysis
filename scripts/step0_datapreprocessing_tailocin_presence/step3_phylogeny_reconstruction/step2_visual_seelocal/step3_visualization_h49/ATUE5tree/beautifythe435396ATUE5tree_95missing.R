# ============================================================
# Supp Fig SX: Diversity within ATUE5 phylotype
# (1) Midpoint-rooted rectangular tree
# (2) Unrooted circular tree with dashed outer lines
# ============================================================

library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(stringr)
library(phangorn)

# ---- Paths ----
treefile <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/scripts/step0_datapreprocessing_tailocin_presence/step3_phylogeny_tree/results_h43_m53/tree/bialleliconly_h43_m53_96samples.min96.phy.treefile"

out_dir   <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/scripts/step0_datapreprocessing_tailocin_presence/step3_phylogeny_tree/step3_visualization_h49/ATUE5tree/"

atue5_list <- file.path("/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/scripts/step0_datapreprocessing_tailocin_presence/step3_phylogeny_tree/step3_visualization_h49/samples_m53_h43.txt")

out_pdf_midroot <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/scripts/step0_datapreprocessing_tailocin_presence/step3_phylogeny_tree/step3_visualization_h49/ATUE5tree/supfigh43_ATUE5_tree_midroot.pdf"
out_pdf_unrooted <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/scripts/step0_datapreprocessing_tailocin_presence/step3_phylogeny_tree/step3_visualization_h49/ATUE5tree/supfigh43_ATUE5_tree_unrooted_dashed.pdf"

dir.create(dirname(out_pdf_midroot), showWarnings = FALSE, recursive = TRUE)

# ---- Colors ----
col_mod  <- "#1b9e77"   # green – Modern ATUE5
  col_hist <- "#a6631b"   # brown – Historical ATUE5
    col_branch <- "black"
      
    # ---- Load tree and ATUE5 labels ----
    tr <- read.tree(treefile)
    atue5_labels <- read.table(atue5_list, stringsAsFactors = FALSE)$V1
    
    # ---- Keep only ATUE5 tips ----
    keep_tips <- tr$tip.label[tr$tip.label %in% atue5_labels]
    tr_sub <- drop.tip(tr, setdiff(tr$tip.label, keep_tips))
    tr_sub <- ladderize(tr_sub)
    
    # ============================================================
    # 1️⃣ MIDPOINT-ROOTED RECTANGULAR TREE
    # ============================================================
    
    tr_mid <- phangorn::midpoint(tr_sub)
    tr_mid <- ladderize(tr_mid)
    
    p_mid <- ggtree(tr_mid, layout = "rectangular", size = 0.35, color = col_branch)
    
    tip_df_mid <- p_mid$data %>%
      filter(isTip) %>%
      mutate(
        time_group = case_when(
          str_detect(label, "^p") ~ "Modern ATUE5",
          label == "DC3000" ~ "Modern ATUE5",
          TRUE ~ "Historical ATUE5"
        )
      )
    
    p_mid_final <- p_mid +
      geom_tippoint(
        data = tip_df_mid,
        aes(color = time_group, shape = time_group),
        size = 2.6, stroke = 0.25, alpha = 0.95
      ) +
      scale_color_manual(
        name = NULL,
        values = c("Historical ATUE5" = col_hist, "Modern ATUE5" = col_mod)
      ) +
      scale_shape_manual(
        name = NULL,
        values = c("Historical ATUE5" = 16, "Modern ATUE5" = 17)
      ) +
      guides(
        color = guide_legend(override.aes = list(size = 3.5, alpha = 1)),
        shape = guide_legend(override.aes = list(size = 3.5, alpha = 1))
      ) +
      theme_tree2() +
      theme(
        legend.position = c(0.85, 0.85),
        legend.text = element_text(size = 13, face = "bold"),
        plot.margin = margin(8, 10, 8, 10),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15)
      ) 
    ggsave(out_pdf_midroot, plot = p_mid_final,
           width = 10, height = 8, dpi = 600, useDingbats = FALSE)
    cat("✅ Midpoint-rooted tree saved to:", out_pdf_midroot, "\n")
    
    # ============================================================
    # 2️⃣ UNROOTED CIRCULAR TREE WITH DASHED SPOKES
    # ============================================================
    
    p_unroot <- ggtree(tr_sub, layout = "circular", size = 0.4, color = col_branch)
    
    tip_df_unroot <- p_unroot$data %>%
      filter(isTip) %>%
      mutate(
        time_group = case_when(
          str_detect(label, "^p") ~ "Modern ATUE5",
          label == "mDC3000" ~ "Modern ATUE5",
          TRUE ~ "Historical ATUE5"
        ),
        x_outer = x + 0.01,  # offset for dashed line
        y_outer = y
      )
    
    # ---- Add dashed outer lines (spokes) ----
    seg_layer <- geom_segment(
      data = tip_df_unroot,
      aes(x = x, y = y, xend = x_outer, yend = y_outer),
      linetype = "dotted", linewidth = 0.18, inherit.aes = FALSE, show.legend = FALSE
    )
    
    # ---- Add tip points ----
    dot_layer <- geom_point(
      data = tip_df_unroot,
      aes(x = x_outer, y = y_outer, color = time_group, shape = time_group),
      size = 3.0, stroke = 0.4, alpha = 0.95, inherit.aes = FALSE
    )
    
    p_unroot_final <- p_unroot +
      seg_layer + dot_layer +
      scale_color_manual(
        name = NULL,
        values = c("Historical ATUE5" = col_hist, "Modern ATUE5" = col_mod)
      ) +
      scale_shape_manual(
        name = NULL,
        values = c("Historical ATUE5" = 16, "Modern ATUE5" = 17)
      ) +
      guides(
        color = guide_legend(override.aes = list(size = 4, alpha = 1)),
        shape = guide_legend(override.aes = list(size = 4, alpha = 1))
      ) +
      theme_void() +
      theme(
        legend.position = c(0.75, 0.65),
        legend.text = element_text(size = 13, face = "bold"),
        plot.margin = margin(10, 10, 10, 10),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15)
      ) 
    ggsave(out_pdf_unrooted, plot = p_unroot_final,
           width = 10, height = 10, dpi = 600, useDingbats = FALSE)
    cat("✅ Unrooted dashed-line tree saved to:", out_pdf_unrooted, "\n")