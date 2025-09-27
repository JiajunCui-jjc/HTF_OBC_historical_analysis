# =========================
# Beautify ATUE5 tree plot (Conference-ready)
# =========================

# ---- Libraries ----
library(ape)
library(ggplot2)
library(ggtree)   # from Bioconductor
library(dplyr)
library(readr)
library(stringr)
library(phangorn)

# ---- Paths ----
out_dir   <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step3_h46_m85_tree/step3_visualization/"
treefile  <- file.path("/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/scripts/step0_datapreprocessing_tailocin_presence/step3_h46_m85_phylogeny_reconstruction/step3_visualization/bialleliconly_filtered.anew2024_all131.min131.phy.treefile")
atue5_file <- file.path("/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/scripts/step0_datapreprocessing_tailocin_presence/step3_h46_m85_phylogeny_reconstruction/step3_visualization/samples_m55_h40.txt")
outfile_pdf <- file.path(out_dir, "material2_treeofATUE5andnon.pdf")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Colors & shapes ----
col_mod <- "#1b9e77"       # Modern (teal green)
  col_hist <- "#a6631b"      # Historical (brown)
    col_nonatue5 <- "grey50"   # non-ATUE5 grey
      shape_all <- 16            # circle for all
      
      # ---- Read tree ----
      tr <- ape::read.tree(treefile)
      
      # ---- Midpoint root + ladderize ----
      tr <- phangorn::midpoint(tr)
      tr <- ladderize(tr)
      
      # ---- Base tree ----
      p0 <- ggtree(tr, layout = "circular", size = 0.85)
      
      # ---- Build tip metadata ----
      atue5_labels <- read_table(atue5_file, col_names = FALSE)$X1
      
      tip_df <- p0$data %>%
        filter(isTip) %>%
        mutate(
          era = case_when(
            str_detect(label, "^p") ~ "Modern ATUE5",
            label == "mDC3000"     ~ "Modern ATUE5",
            TRUE                   ~ "Historical ATUE5"
          ),
          era = factor(era, levels = c("Modern ATUE5", "Historical ATUE5")),
          group_ATUE5 = ifelse(label %in% atue5_labels, "ATUE5", "non-ATUE5"),
          group_ATUE5 = factor(group_ATUE5, levels = c("ATUE5", "non-ATUE5"))
        )
      
      # ---- Prepare coords ----
      outer_offset <- 0.15
      tip_outer <- tip_df %>%
        mutate(
          x_outer = ifelse(group_ATUE5 == "ATUE5", x + outer_offset, x),
          y_outer = y
        )
      
      # ---- Dashed spokes only for ATUE5 ----
      seg_layer <- geom_segment(
        data = tip_outer %>% filter(group_ATUE5 == "ATUE5"),
        aes(x = x, y = y, xend = x_outer, yend = y_outer),
        linetype = "dotted", linewidth = 0.3, inherit.aes = FALSE
      )
      
      # ---- Dot layer ----
      dot_layer <- geom_point(
        data = tip_outer,
        aes(
          x = x_outer, y = y_outer,
          color = ifelse(group_ATUE5 == "ATUE5", as.character(era), "non-ATUE5")
        ),
        shape = shape_all,
        size = 2.6, stroke = 0.3, alpha = 0.95, inherit.aes = FALSE
      )
      
      # ---- Assemble & style ----
      p <- p0 +
        seg_layer +
        dot_layer +
        scale_color_manual(
          values = c(
            "Modern ATUE5" = col_mod,
            "Historical ATUE5" = col_hist,
            "non-ATUE5" = col_nonatue5
          )
        ) +
        guides(
          color = guide_legend(override.aes = list(size = 4, alpha = 1), order = 1)
        ) +
        theme_void() +
        theme(
          legend.position = "right",
          legend.title = element_blank(),
          legend.text = element_text(size = 18),
          plot.margin = margin(4, 4, 4, 4)
        ) +
        xlim_tree(c(NA, max(tip_outer$x_outer) + 0.02))
      
      # ---- Flip vertically (180° rotation) ----
      p <- rotate_tree(p, angle = 180)
      
      # ---- Save ----
      ggsave(outfile_pdf, plot = p,
             width = 12, height = 12, device = grDevices::pdf, useDingbats = FALSE)
      
      # ---- Save sample table ----
      modern_hist_table <- tip_df %>%
        select(label, era, group_ATUE5) %>%
        arrange(era, group_ATUE5, label)
      
      write.table(modern_hist_table, file.path(out_dir, "Modern_vs_Historical_samples.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)