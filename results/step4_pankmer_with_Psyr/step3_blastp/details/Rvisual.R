#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(stringr)
})

# -------------------------
# Input / output
# -------------------------
infile <- "/Users/cuijiajun/Desktop/others/tmphernan/2026_PNAS_revision/results/T1T2_assignATUE5/step3_blastp/blastp_best_hit.wide.tsv"
outpdf <- "/Users/cuijiajun/Desktop/others/tmphernan/2026_PNAS_revision/results/T1T2_assignATUE5/step3_blastp/HTF_alignment_coverage_vs_identity_byHTFlength_shapeT1.pdf"

# -------------------------
# HTF length color scheme
# -------------------------
length_colors <- c(
  "1830" = "#f6d6ff",
  "1383" = "#638ccc",
  "1803" = "#800233",
  "1245" = "#f9d42a"
)

shape_map <- c(
  "T1" = 16,   # filled
  "T2" = 1,    # open
  "T3" = 1
)

# -------------------------
# Read and parse
# -------------------------
df <- read_tsv(infile, show_col_types = FALSE) %>%
  separate(`aligned_len/ref_len`,
           into = c("aligned_len", "ref_len"),
           sep = "/",
           convert = TRUE) %>%
  mutate(
    coverage_pct = 100 * aligned_len / ref_len,
    best_identity = as.numeric(best_identity),
    HTF_length = str_extract(query, "(?<=trimHTF)[0-9]+"),
    HTF_length = factor(HTF_length, levels = names(length_colors)),
    best_ref_type = factor(best_ref_type, levels = c("T1","T2","T3"))
  ) %>%
  filter(!is.na(HTF_length))

# -------------------------
# Plot
# -------------------------
p <- ggplot(df, aes(
  x = coverage_pct,
  y = best_identity,
  color = HTF_length,
  shape = best_ref_type
)) +
  geom_point(size = 2.8, alpha = 0.85, stroke = 1) +
  scale_color_manual(
    values = length_colors,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = shape_map,
    name = "Best ref type"
  ) +
  labs(
    x = "Alignment coverage (%)",
    y = "Percent identity",
    color = "HTF length (bp)",
    title = "HTF best-hit alignment coverage vs identity"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# -------------------------
# Save
# -------------------------
ggsave(outpdf, p, width = 6.8, height = 5.2, useDingbats = FALSE)

message("Saved plot to: ", outpdf)