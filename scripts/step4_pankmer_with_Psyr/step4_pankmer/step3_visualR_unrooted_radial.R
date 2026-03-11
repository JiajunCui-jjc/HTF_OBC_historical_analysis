#!/usr/bin/env Rscript

setwd("/Users/cuijiajun/Desktop/others/tmphernan/2026_PNAS_revision/results/pankmer")

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(grid)    # unit()
  library(ggtext)  # element_markdown for Helvetica-consistent italics
})

# ==============================
# Output dirs
# ==============================
dir_S1 <- "final_figs/figS1"
dir_S2 <- "final_figs/figS2"
dir.create(dir_S1, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_S2, showWarnings = FALSE, recursive = TRUE)

# ==============================
# EXACT COLORS YOU GAVE
# ==============================
T1_COL <- "#1a9e78"
T2_COL <- "#7570b3"
NA_COL <- "#e6e6e6"

phy_order <- c(
  "P. viridiflava ATUE5",
  "P. syringae 1a",
  "P. syringae 1b",
  "P. syringae 2a",
  "P. syringae 2b",
  "P. syringae 2d",
  "P. syringae 3",
  "P. syringae 5",
  "P. syringae 6",
  "P. syringae 10a",
  "P. syringae 10b",
  "NA"
)

phy_scale_vals <- c(
  "P. viridiflava ATUE5" = "#2ba02b",
  "P. syringae 1a"  = "#464d86",
  "P. syringae 1b"  = "#7680ca",
  "P. syringae 2a"  = "#ca743e",
  "P. syringae 2b"  = "#c95573",
  "P. syringae 2d"  = "#d62628",
  "P. syringae 3"   = "#8c9f3d",
  "P. syringae 5"   = "#9467bd",
  "P. syringae 6"   = "#6cb990",
  "P. syringae 10a" = "#0f101f",
  "P. syringae 10b" = "#8e93b0",
  "NA"              = NA_COL
)

# ==============================
# Labels (markdown italics; keeps Helvetica via ggtext)
# ==============================
phy_labels_ital <- c(
  "P. viridiflava ATUE5" = "<i>P. viridiflava</i> ATUE5",
  "P. syringae 1a"  = "<i>P. syringae</i> 1a",
  "P. syringae 1b"  = "<i>P. syringae</i> 1b",
  "P. syringae 2a"  = "<i>P. syringae</i> 2a",
  "P. syringae 2b"  = "<i>P. syringae</i> 2b",
  "P. syringae 2d"  = "<i>P. syringae</i> 2d",
  "P. syringae 3"   = "<i>P. syringae</i> 3",
  "P. syringae 5"   = "<i>P. syringae</i> 5",
  "P. syringae 6"   = "<i>P. syringae</i> 6",
  "P. syringae 10a" = "<i>P. syringae</i> 10a",
  "P. syringae 10b" = "<i>P. syringae</i> 10b",
  "NA"              = "NA"
)

type_labels_ital <- c(
  "P. syringae T1" = "<i>P. syringae</i> T1",
  "P. syringae T2" = "<i>P. syringae</i> T2",
  "P. viridiflava ATUE5 T1" = "<i>P. viridiflava</i> ATUE5 T1",
  "P. viridiflava ATUE5 T2" = "<i>P. viridiflava</i> ATUE5 T2",
  "NA" = "NA"
)

# ==============================
# PNAS-like style (~7 pt) + Helvetica-consistent markdown legend
# ==============================
theme_pnas7 <- function() {
  theme(
    text = element_text(family = "Helvetica", size = 7),
    
    # IMPORTANT: legend markdown renderer to preserve Helvetica + italics
    legend.title = ggtext::element_markdown(family = "Helvetica", size = 7),
    legend.text  = ggtext::element_markdown(family = "Helvetica", size = 7),
    
    legend.position = "right",
    legend.key.height = unit(3.0, "mm"),
    legend.key.width  = unit(3.0, "mm"),
    legend.spacing.y  = unit(1.2, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(2, 2, 2, 2),
    
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank()
  )
}

# ==============================
# Metadata
# ==============================
meta <- read_tsv("metainfo/metainfo.ALL.tsv", show_col_types = FALSE) %>%
  mutate(
    SAMPLE = as.character(SAMPLE),
    Taxa = trimws(as.character(Taxa)),
    Phylogroup = trimws(as.character(Phylogroup)),
    Tailocin_type = toupper(trimws(as.character(Tailocin_type)))
  ) %>%
  mutate(
    Phylogroup = na_if(Phylogroup, "NA"),
    Phylogroup = na_if(Phylogroup, ""),
    Tailocin_type = na_if(Tailocin_type, "NA"),
    Tailocin_type = na_if(Tailocin_type, "")
  )

bad_type <- meta %>%
  filter(is.na(Tailocin_type) | !Tailocin_type %in% c("T1", "T2")) %>%
  distinct(SAMPLE, Tailocin_type)

if (nrow(bad_type) > 0) {
  stop(
    "Found Tailocin_type not in {T1,T2} or NA. Fix metainfo. Examples:\n",
    paste(capture.output(print(head(bad_type, 20))), collapse = "\n")
  )
}

# ==============================
# Helper: read tree + join meta
# ==============================
get_tree_data <- function(treefile) {
  tr <- read.tree(treefile)
  p0 <- ggtree(tr, layout = "unrooted")
  
  df <- p0$data %>% left_join(meta, by = c("label" = "SAMPLE"))
  
  bad_join <- df %>% filter(isTip) %>% filter(is.na(Taxa) | is.na(Tailocin_type))
  if (nrow(bad_join) > 0) {
    stop(
      "Some tips missing Taxa/Tailocin_type after join. Examples:\n",
      paste(capture.output(print(head(bad_join %>% select(label, Taxa, Tailocin_type), 20))), collapse = "\n")
    )
  }
  
  list(base_plot = p0, df = df)
}

# ==============================
# Tip category builders
# ==============================
build_tip_phylogroup <- function(df) {
  df %>%
    filter(isTip) %>%
    mutate(
      cat = case_when(
        Taxa == "P_viridiflava" ~ "P. viridiflava ATUE5",
        Taxa == "P_syringae" & !is.na(Phylogroup) ~ paste0("P. syringae ", Phylogroup),
        Taxa == "P_syringae" ~ "NA",
        TRUE ~ "NA"
      ),
      cat = ifelse(cat %in% names(phy_scale_vals), cat, "NA"),
      cat = factor(cat, levels = phy_order)
    )
}

build_tip_type <- function(df) {
  tip_df <- df %>%
    filter(isTip) %>%
    mutate(
      cat = case_when(
        Taxa == "P_syringae" & Tailocin_type == "T1" ~ "P. syringae T1",
        Taxa == "P_syringae" & Tailocin_type == "T2" ~ "P. syringae T2",
        Taxa == "P_viridiflava" & Tailocin_type == "T1" ~ "P. viridiflava ATUE5 T1",
        Taxa == "P_viridiflava" & Tailocin_type == "T2" ~ "P. viridiflava ATUE5 T2",
        TRUE ~ "NA"
      ),
      outline_col = ifelse(cat == "P. viridiflava ATUE5 T1", "black", "transparent")
    )
  
  type_order <- c(
    "P. syringae T1",
    "P. syringae T2",
    "P. viridiflava ATUE5 T1",
    "P. viridiflava ATUE5 T2",
    "NA"
  )
  tip_df$cat <- factor(tip_df$cat, levels = type_order)
  tip_df
}

# ==============================
# Plotters (no title, no scale)
# ==============================
plot_unrooted_phy <- function(treefile, point_size = 2.2, stroke_lw = 0.22) {
  obj <- get_tree_data(treefile)
  p0  <- obj$base_plot
  df  <- obj$df
  tip_df <- build_tip_phylogroup(df)
  
  p0 +
    theme_tree() +
    theme_pnas7() +
    geom_tippoint(
      data = tip_df,
      aes(fill = cat),
      shape = 21, size = point_size,
      colour = "transparent", stroke = stroke_lw
    ) +
    scale_fill_manual(
      values = phy_scale_vals,
      drop = FALSE,
      name = "Ring: phylogroup",
      breaks = phy_order,
      labels = phy_labels_ital
    ) +
    guides(
      fill = guide_legend(
        override.aes = list(shape = 21, colour = "transparent", stroke = 0)
      )
    )
}

plot_unrooted_type <- function(treefile, point_size = 2.2, stroke_lw = 0.5) {
  obj <- get_tree_data(treefile)
  p0  <- obj$base_plot
  df  <- obj$df
  tip_df <- build_tip_type(df)
  
  type_scale_vals <- c(
    "P. syringae T1" = T1_COL,
    "P. syringae T2" = T2_COL,
    "P. viridiflava ATUE5 T1" = T1_COL,
    "P. viridiflava ATUE5 T2" = T2_COL,
    "NA" = NA_COL
  )
  
  # legend: remove ATUE5 T2 and NA
  legend_keep <- c("P. syringae T1", "P. syringae T2", "P. viridiflava ATUE5 T1")
  
  # Legend border only for ATUE5 T1
  legend_border <- ifelse(legend_keep == "P. viridiflava ATUE5 T1", "black", "transparent")
  legend_stroke <- ifelse(legend_keep == "P. viridiflava ATUE5 T1", stroke_lw, 0)
  
  p0 +
    theme_tree() +
    theme_pnas7() +
    geom_tippoint(
      data = tip_df,
      aes(fill = cat, colour = outline_col),
      shape = 21, size = point_size,
      stroke = stroke_lw
    ) +
    scale_fill_manual(
      values = type_scale_vals,
      drop = FALSE,
      name = "Ring: tail fiber type",
      breaks = legend_keep,
      labels = type_labels_ital[legend_keep]
    ) +
    scale_colour_identity(guide = "none") +
    guides(
      fill = guide_legend(
        override.aes = list(
          shape = 21,
          colour = legend_border,
          stroke = legend_stroke
        )
      )
    )
}

# ==============================
# Trees
# ==============================
tree_without <- "NJ_trees/without_HTF.distance.jaccard.NJ.newick"
tree_only   <- "NJ_trees/only_HTF.distance.jaccard.NJ.newick"

# ==============================
# Make plots
# ==============================
pS1_without <- plot_unrooted_phy(tree_without)
pS1_only    <- plot_unrooted_phy(tree_only)

pS2_without <- plot_unrooted_type(tree_without)
pS2_only    <- plot_unrooted_type(tree_only)

# ==============================
# Save PDFs (size includes legend)
# Use cairo_pdf for Helvetica embedding + consistent rendering
# ==============================
W <- 5.6
H <- 3.8

ggsave(file.path(dir_S1, "withoutHTF.unrooted_radial.phylogroup.pdf"),
       pS1_without, width = W, height = H)
ggsave(file.path(dir_S1, "onlyHTF.unrooted_radial.phylogroup.pdf"),
       pS1_only, width = W, height = H)

ggsave(file.path(dir_S2, "withoutHTF.unrooted_radial.type.pdf"),
       pS2_without, width = W, height = H)
ggsave(file.path(dir_S2, "onlyHTF.unrooted_radial.type.pdf"),
       pS2_only, width = W, height = H)

message("DONE. Helvetica-consistent italics via ggtext; S2 legend dropped ATUE5 T2 + NA; no titles/scales.")
message("  ", normalizePath(dir_S1))
message("  ", normalizePath(dir_S2))
