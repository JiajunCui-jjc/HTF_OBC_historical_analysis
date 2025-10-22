# === Libraries ===
library(data.table)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(stringr)

# === Working directory ===
setwd('/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/figures_tables/suppfig1_confidence_matrix_HTF/step2_Rvisual/')

# === Input files ===
merged_file <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/figures_tables/supptables/supptable2sum/merged_final_HTF_tailocin.txt"
modern53_file <- "modern53-infig.txt"

# === Load merged dataset ===
df_merge <- fread(merged_file)

# --- Extract numeric HTF lengths from haplotype columns ---
df_all <- df_merge %>%
  mutate(
    HTF_by_local = str_extract(HTF_haplotype_bylocalassembly, "\\((\\d+)\\)") %>% str_replace_all("[()]", ""),
    HTF_by_kmer  = str_extract(HTF_haplotype_bykmer, "\\((\\d+)\\)") %>% str_replace_all("[()]", "")
  ) %>%
  select(Strain = sample,
         `By Local Assembly` = HTF_by_local,
         `By Kmer` = HTF_by_kmer)

# === Define base HTF color mapping ===
base_colors <- c(
  "1830" = "#f6d6ff",
  "1383" = "#638ccc",
  "1803" = "#800233",
  "1245" = "#f9d42a"
)

# === Helper plotting function ===
plot_htf_heatmap <- function(df_mod, output_pdf, highlight_historical = FALSE) {
  
  # --- Reshape to long → wide ---
  df_long <- df_mod %>%
    pivot_longer(cols = -Strain, names_to = "method", values_to = "length")
  df_wide <- df_long %>%
    pivot_wider(names_from = Strain, values_from = length)
  
  # --- Assign colors ---
  all_lens <- unique(unlist(df_wide[, -1]))
  all_lens <- all_lens[!is.na(all_lens)]
  other_lens <- setdiff(all_lens, names(base_colors))
  if (length(other_lens) > 0) {
    other_names <- paste0("other_", other_lens)
    names(other_names) <- other_lens
    other_colors <- setNames(colorRampPalette(c("#bbbbbb", "#444444"))(length(other_lens)), other_names)
  } else {
    other_colors <- character(0)
  }
  all_colors <- c(base_colors, other_colors, "NA" = "white")
  
  # --- Convert to color codes (NA → white) ---
  df_display <- df_wide
  for (i in 2:ncol(df_display)) {
    df_display[[i]] <- sapply(df_display[[i]], function(val) {
      if (is.na(val)) "NA"
      else if (val %in% names(base_colors)) val
      else "NA"
    })
  }
  
  # --- Sorting (By Kmer present first) ---
  df_sort <- df_mod %>%
    mutate(
      kmer_present = !is.na(`By Kmer`),
      kmer_group = ifelse(`By Kmer` %in% names(base_colors),
                          `By Kmer`,
                          paste0("other_", `By Kmer`)),
      has_mapping = ifelse(is.na(`By Local Assembly`), 0, 1)
    ) %>%
    arrange(desc(kmer_present),
            factor(kmer_group, levels = c(names(base_colors),
                                          sort(setdiff(unique(kmer_group), names(base_colors))))),
            desc(has_mapping),
            Strain)
  
  df_sort <- df_sort %>% filter(Strain %in% colnames(df_display))
  strain_order <- df_sort$Strain
  
  # --- Label colors ---
  label_colors <- sapply(strain_order, function(s) {
    if (grepl("^p", s)) {
      "black"  # modern
    } else if (highlight_historical) {
      "#804111"  # historical
    } else {
      "#cccccc"
    }
  })
  names(label_colors) <- strain_order
  
  # --- Matrix for plotting ---
  mat <- as.matrix(df_display[, strain_order])
  rownames(mat) <- df_wide[[1]]
  
  # --- Draw heatmap ---
  ht <- Heatmap(
    mat,
    name = "HTF length",
    col = all_colors,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    column_names_side = "bottom",
    column_names_rot = 45,
    column_names_gp = gpar(col = label_colors),
    heatmap_legend_param = list(title = "HTF length"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.rect(x, y, width, height, gp = gpar(fill = fill, col = "grey"))
    }
  )
  
  pdf(output_pdf, width = length(strain_order) * 0.25 + 2, height = 4)
  draw(ht)
  dev.off()
  cat("✅ Saved:", output_pdf, "\n")
}

# === 🟢 Figure 1: All samples (97 total) ===
plot_htf_heatmap(df_all, "modern57_h40_htf_heatmap_by_length_all.pdf", highlight_historical = TRUE)

# === 🟢 Figure 2: Modern53 subset ===
if (file.exists(modern53_file)) {
  modern53_list <- fread(modern53_file, header = FALSE)[[1]]
  df_mod53 <- df_all %>% filter(Strain %in% modern53_list)
  plot_htf_heatmap(df_mod53, "modern53_htf_heatmap_by_length.pdf", highlight_historical = FALSE)
}

# === 🟢 Figure 3: Historical 34 (excluding 64.GBR*) ===
df_hist34 <- df_all %>%
  filter(!grepl("^p", Strain) & 
           !is.na(`By Kmer`))
plot_htf_heatmap(df_hist34, "hist34_htf_heatmap_by_length.pdf", highlight_historical = TRUE)

# === 🟢 Figure 4: Historical 40 (auto-detected) ===
df_hist40 <- df_all %>%
  filter(!grepl("^p", Strain))
plot_htf_heatmap(df_hist40, "hist40_htf_heatmap_by_length.pdf", highlight_historical = TRUE)
