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

# ---- Paths ----
out_dir   <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_Aug_balgiumtalk_materials_jiajun/material2_treeofATUE5andnon/howto"
treefile  <- file.path(out_dir, "bialleliconly_filtered.anew2024_all131.min131.phy.treefile")
atue5_file <- file.path(out_dir, "samples_m55_h40.txt")
outfile_pdf <- file.path(out_dir, "material2_treeofATUE5andnon.pdf")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Colors & shapes ----
col_mod <- "#1b9e77"   # Modern (teal green)
  col_hist <- "#a6631b"  # Historical (brown)
    shape_atue5     <- 16  # circle
    shape_non_atue5 <- 17  # triangle
    
    # ---- Read tree ----
    tr <- ape::read.tree(treefile)
    
    # ---- Midpoint root ----
    # ---- Midpoint root ----
    library(phangorn)   # needed for midpoint()
    tr <- phangorn::midpoint(tr)     # midpoint root
    tr <- ladderize(tr)              # ladderize for aesthetics
    # ---- Prepare tip metadata ----
    tips <- tibble::tibble(label = tr$tip.label)
    
    # ---- Read ATUE5 list ----
    atue5_raw <- tryCatch(read_tsv(atue5_file, col_names = FALSE, show_col_types = FALSE),
                          error = function(e) tryCatch(read_csv(atue5_file, col_names = FALSE, show_col_types = FALSE),
                                                       error = function(e2) read_table(atue5_file, col_names = FALSE)))
    
    if (ncol(atue5_raw) == 1) {
      colnames(atue5_raw) <- c("label")
      atue5_raw$era <- NA_character_
    } else {
      colnames(atue5_raw)[1:2] <- c("label", "era")
    }
    
    # ---- Build metadata ----
    meta <- tips %>%
      left_join(atue5_raw %>% select(label), by = "label") %>%
      mutate(
        group_ATUE5 = if_else(label %in% atue5_raw$label, "ATUE5", "non-ATUE5")
      )
    
    # ---- Era assignment (custom rules) ----
    meta <- meta %>%
      mutate(
        era = case_when(
          str_detect(label, regex("^p", ignore_case = FALSE)) ~ "Modern",  # lowercase "p" start
          label == "mDC3000" ~ "Modern",
          TRUE ~ "Historical"
        ),
        era = factor(era, levels = c("Modern", "Historical"))
      )
    
    # ---- Shapes and colors ----
    meta <- meta %>%
      mutate(
        shape = if_else(group_ATUE5 == "ATUE5", shape_atue5, shape_non_atue5),
        color = if_else(era == "Modern", col_mod, col_hist),
        group_ATUE5 = factor(group_ATUE5, levels = c("ATUE5", "non-ATUE5")),
        era = factor(era, levels = c("Modern", "Historical"))
      )
    
    # ---- Plot (conference ready) ----
   # p <- ggtree(tr, layout = "rectangular", size = 0.45) %<+% meta +
  #unrooted
    # Radial (circular) unrooted tree
    # --- Base tree (no tip dots) ---
    p0 <- ggtree(tr, layout = "circular", size = 0.85)
    
    # --- Build tip metadata directly from labels (robust; no reliance on %<+%) ---
    tip_df <- p0$data %>%
      dplyr::filter(isTip) %>%
      dplyr::mutate(
        # era by your rules
        era = dplyr::case_when(
          stringr::str_detect(label, "^p") ~ "Modern",
          label == "mDC3000"              ~ "Modern",
          TRUE                            ~ "Historical"
        ),
        era = factor(era, levels = c("Modern", "Historical")),
        # ATUE5 set membership using your list
        group_ATUE5 = ifelse(label %in% atue5_raw$label, "ATUE5", "non-ATUE5"),
        group_ATUE5 = factor(group_ATUE5, levels = c("ATUE5", "non-ATUE5"))
      )
    
    # --- Move dots to an outer ring ---
    outer_offset <- 0.15  # try 0.15–0.40
    tip_outer <- tip_df %>%
      dplyr::mutate(x_outer = x + outer_offset)
    
    # --- (Optional) dotted spokes from tip to outer dot ---
    seg_layer <- geom_segment(
      data = tip_outer,
      aes(x = x, y = y, xend = x_outer, yend = y),
      linetype = "dotted", linewidth = 0.3, inherit.aes = FALSE
    )
    
    # --- Outer dots ONLY (no dots at original tips) ---
    dot_layer <- geom_point(
      data = tip_outer,
      aes(x = x_outer, y = y, color = era, shape = group_ATUE5),
      size = 2.6, stroke = 0.3, alpha = 0.95, inherit.aes = FALSE
    )
    
    # --- Assemble & style ---
    p <- p0 +
      seg_layer +               # remove this line if you want no spokes
      dot_layer +
      scale_color_manual(values = c(Modern = col_mod, Historical = col_hist)) +
      scale_shape_manual(values = c(ATUE5 = shape_atue5, `non-ATUE5` = shape_non_atue5)) +
      guides(
        color = guide_legend(override.aes = list(size = 4, alpha = 1), order = 1),
        shape = guide_legend(override.aes = list(size = 4, alpha = 1), order = 2)
      ) +
      theme_void() +
      theme(
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        plot.margin = margin(4, 4, 4, 4)
      ) +
      # ensure outer dots aren’t clipped
      xlim_tree(c(NA, max(tip_outer$x_outer) + 0.02))
    
    # --- Save bigger page so the tree looks larger ---
    ggsave(outfile_pdf, plot = p,
           width = 12, height = 12, device = grDevices::pdf, useDingbats = FALSE)
    # ---- Save Modern/Historical lists ----
    modern_hist_table <- meta %>%
      select(label, era, group_ATUE5) %>%
      arrange(era, group_ATUE5, label)
    
    write.table(modern_hist_table, file.path(out_dir, "Modern_vs_Historical_samples.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    