# =========================
# Beautify ATUE5 tree plot (Conference-ready, branch + tip colors)
# =========================

# ---- Libraries ----
library(ape)
library(ggplot2)
library(ggtree)
library(dplyr)
library(readr)
library(stringr)
library(phangorn)
library(ggnewscale)  # <-- allows multiple color scales

# ---- Paths ----
out_dir   <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step3_h46_m85_tree/step3_visualization/"
treefile  <- file.path("/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/scripts/step0_datapreprocessing_tailocin_presence/step3_h46_m85_phylogeny_reconstruction/step3_visualization/bialleliconly_filtered.anew2024_all131.min131.phy.treefile")
atue5_file <- file.path("/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/scripts/step0_datapreprocessing_tailocin_presence/step3_h46_m85_phylogeny_reconstruction/step3_visualization/samples_m55_h40.txt")
outfile_pdf <- file.path(out_dir, "material2_treeofATUE5andnon.pdf")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Colors ----
col_mod_atue5 <- "#1b9e77"   # green (modern ATUE5)
  col_hist_atue5 <- "#a6631b"  # brown (historical ATUE5)
    col_non <- "grey80"          # non-ATUE5 (grey)
      branch_col <- c("ATUE5" = "black", "non-ATUE5" = "grey80")
      
      # ---- Read input ----
      tr <- ape::read.tree(treefile)
      atue5_labels <- read_table(atue5_file, col_names = FALSE)$X1
      
      # ---- Midpoint root + ladderize ----
      tr <- phangorn::midpoint(tr)
      tr <- ladderize(tr)
      
      # ---- Identify branch type ----
      has_atue5 <- function(node, tree, atue5_labels) {
        desc <- phangorn::Descendants(tree, node, type = "tips")[[1]]
        any(tree$tip.label[desc] %in% atue5_labels)
      }
      branch_df <- tibble(
        node = 1:(Nnode(tr) + length(tr$tip.label)),
        branch_type = sapply(1:(Nnode(tr) + length(tr$tip.label)), function(n)
          if (has_atue5(n, tr, atue5_labels)) "ATUE5" else "non-ATUE5")
      )
      
      # ---- Plot base tree with branch colors ----
      p <- ggtree(tr, layout = "circular", size = 0.5) %<+% branch_df +
        aes(color = branch_type) +
        scale_color_manual(values = branch_col, guide = "none")
      
      # ---- Add new color scale for tips ----
      p <- p + ggnewscale::new_scale_color()
      
      # ---- Tip metadata ----
      tip_df <- p$data %>%
        filter(isTip) %>%
        mutate(
          time_group = case_when(
            str_detect(label, "^p") ~ "Modern",
            label == "mDC3000" ~ "Modern",
            TRUE ~ "Historical"
          ),
          group_ATUE5 = ifelse(label %in% atue5_labels, "ATUE5", "non-ATUE5"),
          full_group = case_when(
            group_ATUE5 == "ATUE5" & time_group == "Modern" ~ "Modern ATUE5",
            group_ATUE5 == "ATUE5" & time_group == "Historical" ~ "Historical ATUE5",
            group_ATUE5 == "non-ATUE5" & time_group == "Modern" ~ "Modern non-ATUE5",
            TRUE ~ "Historical non-ATUE5"
          ),
          full_group = factor(
            full_group,
            levels = c("Historical ATUE5", "Modern ATUE5",
                       "Historical non-ATUE5", "Modern non-ATUE5")
          )
        )
      
      # ---- Prepare coords ----
      outer_offset <- 0.22
      tip_outer <- tip_df %>%
        mutate(
          x_outer = ifelse(group_ATUE5 == "ATUE5", x + outer_offset, x),
          y_outer = y
        )
      
      # ---- Dashed spokes for ATUE5 ----
      seg_layer <- geom_segment(
        data = tip_outer %>% filter(group_ATUE5 == "ATUE5"),
        aes(x = x, y = y, xend = x_outer, yend = y_outer),
        linetype = "dotted", linewidth = 0.15, inherit.aes = FALSE,
        show.legend = FALSE
      )
      
      # ---- Dot layer ----
      dot_layer <- geom_point(
        data = tip_outer,
        aes(
          x = x_outer, y = y_outer,
          color = full_group,
          shape = full_group,
          size = group_ATUE5
        ),
        stroke = 0.4, alpha = 0.95, inherit.aes = FALSE
      )
      
      # ---- Assemble ----
      p <- p +
        seg_layer +
        dot_layer +
        scale_size_manual(values = c("ATUE5" = 4.5, "non-ATUE5" = 2.8), guide = "none") +
        scale_color_manual(
          name = NULL,
          values = c(
            "Historical ATUE5" = col_hist_atue5,
            "Modern ATUE5" = col_mod_atue5,
            "Historical non-ATUE5" = col_non,
            "Modern non-ATUE5" = col_non
          )
        ) +
        scale_shape_manual(
          name = NULL,
          values = c(
            "Historical ATUE5" = 16,   # ●
            "Modern ATUE5" = 17,       # ▲
            "Historical non-ATUE5" = 16,
            "Modern non-ATUE5" = 17
          )
        ) +
        guides(
          color = guide_legend(
            override.aes = list(size = 4, alpha = 1),
            order = 1
          ),
          shape = guide_legend(
            override.aes = list(size = 4, alpha = 1),
            order = 1
          )
        ) +
        theme_void() +
        theme(
          legend.position = "right",
          legend.title = element_blank(),
          legend.text = element_text(size = 16),
          plot.margin = margin(4, 4, 4, 4)
        ) +
        xlim_tree(c(NA, max(tip_outer$x_outer) + 0.02))
      
      # ---- Rotate tree ----
      p <- rotate_tree(p, angle = 180)
      
      # ---- Save ----
      ggsave(outfile_pdf, plot = p,
             width = 12, height = 12, device = grDevices::pdf, useDingbats = FALSE)