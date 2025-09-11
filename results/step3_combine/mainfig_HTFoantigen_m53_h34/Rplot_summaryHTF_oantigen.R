# === Libraries ===
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(RColorBrewer)

# === Paths ===
modern_keep_file <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_july_figHTFoantigen/step3_combine/modern53-infig.txt"
file <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step3_combine/combined_HTF_oantigen_m57_h35.txt"

# === Load data ===
df <- fread(file)
modern_keep <- fread(modern_keep_file, header = FALSE)$V1

# Keep only: all historical except 64. + modern_keep list
df <- df[( !grepl("^p", sample) & !grepl("^64\\.", sample) ) | sample %in% modern_keep]


# Modern/Historical status: ONLY modern_keep is Modern
mh_status <- ifelse(df$sample %in% modern_keep, "Modern", "Historical")
mh_status <- factor(mh_status, levels = c("Modern", "Historical"))
names(mh_status) <- df$sample

# === LengthGroup from HTF_haplotype (text inside parentheses) ===
# === LengthGroup from HTF_haplotype (text inside parentheses) ===
df[, LengthGroup := sub(".*\\(([^)]+)\\).*", "\\1", HTF_haplotype)]
df$LengthGroup[is.na(df$LengthGroup) | df$LengthGroup == ""] <- "NA"

# === Fixed LengthGroup order: 1830, 1383, 1803, 1245 (NA last) ===
length_order <- c("1830", "1383", "1803", "1245", "NA")

# === Custom HTF order (already in your script) ===
htf_levels <- c(
  "HTF_p7.G11 (1830)",
  "HTF_p25.A12 (1383)",
  "HTF_p25.C2 (1803)",
  "HTF_p5.D5 (1803)",
  "HTF_p26.D6 (1803)",
  "HTF_p23.B8 (1803)",
  "HTF_p21.F9 (1245)"
)

# Use factors so sorting respects your orders
df[, HTF_haplotype := factor(HTF_haplotype, levels = htf_levels)]

# === Sort samples: by LengthGroup (fixed), then HTF (custom), then Modern→Historical, then name ===
df <- df[
  order(
    factor(LengthGroup, levels = length_order),
    
    mh_status,      # Modern first (because you set levels = c("Modern","Historical"))
    rev(sample)
  )
]

# === Extract 6-gene binary matrix ===
gene_cols <- c("wfgD_PA", "rmlC1_PA", "tagG2_PA", "tagH2_PA", "spsA_PA", "espE4_PA")
gene_mat <- as.matrix(df[, ..gene_cols])
rownames(gene_mat) <- df$sample
colnames(gene_mat) <- gsub("_PA$", "", colnames(gene_mat))
gene_mat <- apply(gene_mat, 2, function(x) ifelse(x == "NA", NA, as.numeric(x)))
gene_mat <- t(gene_mat)

# === HTF haplotype annotation (kept as before) ===
# === HTF haplotype annotation with custom order ===
htf_anno <- df$HTF_haplotype
names(htf_anno) <- df$sample

# Custom order
htf_levels <- c(
  "HTF_p7.G11 (1830)",
  "HTF_p25.A12 (1383)",
  "HTF_p25.C2 (1803)",
  "HTF_p5.D5 (1803)",
  "HTF_p26.D6 (1803)",
  "HTF_p23.B8 (1803)",
  "HTF_p21.F9 (1245)"
)

# Convert to factor so ordering is respected
htf_anno <- factor(htf_anno, levels = htf_levels)

# Colors assigned in the same order
htf_colors <- structure(
  c("#f6d6ff", "#638ccc", "#800233", "#800233", "#800233", "#800233","#f9d42a" ),
  names = htf_levels
)

  
# === Sample name color ===
sample_col <- ifelse(grepl("^p", df$sample), "black", "brown")
names(sample_col) <- df$sample

# === Annotations , keep HTF bar) ===
top_anno <- HeatmapAnnotation(
  HTF = htf_anno[colnames(gene_mat)],
  col = list(
    
    HTF = htf_colors
  ),
  annotation_name_side = "left"
)

# === Color function ===
gene_col_fun <- function(val, row, col) {
  if (is.na(val)) {
    return("#bbbbbb")
  } else if (val == 1) {
    return("#666666")
  } else {
    return("#f9f9f9")
  }
}
# Define your bad espE2 samples
bad_espE2_samples <- c("86.NOR_1911_S7", "27.ESP_1975", "PL0139", "p25.C11")

# === Draw main heatmap ===
ht <- Heatmap(
  gene_mat,
  name = "Gene Presence",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 16, fontface = "italic"),
  column_names_gp = gpar(col = sample_col[colnames(gene_mat)], fontsize = 10),
  top_annotation = top_anno,
  show_heatmap_legend = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    val <- gene_mat[i, j]
    gene_name <- rownames(gene_mat)[i]
    sample_name <- colnames(gene_mat)[j]
    
    if (gene_name == "espE4" && sample_name %in% bad_espE2_samples) {
      # Special color for bad espE2 samples
      grid.rect(
        x = x, y = y, width = width, height = height,
        gp = gpar(
          fill = ifelse(is.na(val), "#bbbbbb", "#969696"),
          col = "grey"
        )
      )
      grid.text("*", x = x, y = y, gp = gpar(fontsize = 8, col = "black", fontface = "bold"))
    } else {
      grid.rect(
        x = x, y = y, width = width, height = height,
        gp = gpar(
          fill = gene_col_fun(val, rownames(gene_mat)[i], colnames(gene_mat)[j]),
          col = "grey"
        )
      )
    }
  },
  column_title = "HTF haplotype - O-antigen OBC genes presence/absence (m53 + h34)"
)
pdf("/Users/cuijiajun/Desktop/others/tmphernan/2025_july_figHTFoantigen/step3_combine/m53_h34_heatmap_sixgene_only_with_HTFanno.pdf",
    width = ncol(gene_mat) * 0.25 + 3, height = 5)
draw(ht)
dev.off()

cat("✅ Done: 6-gene binary heatmap generated.\n")

# === Extra Modern/Historical color bar for PPT (aligned to heatmap columns) ===
# align M/H to heatmap columns
mh_status_vec <- mh_status[colnames(gene_mat)]

mh_matrix <- matrix(1, nrow = 1, ncol = length(mh_status_vec))
rownames(mh_matrix) <- "Status"
colnames(mh_matrix) <- colnames(gene_mat)

mh_colors <- c(
  Modern = "#1b9e77",      # Teal green
  Historical = "#d95f02"   # Earthy orange
)
mh_col_fun <- function(j, i, x, y, width, height, fill) {
  cn <- colnames(mh_matrix)[j]                   # j = column index
  st <- as.character(mh_status_vec[cn])          # "Modern" or "Historical"
  grid.rect(x = x, y = y, width = width, height = height,
            gp = gpar(fill = mh_colors[st], col = NA))
}

mh_ht <- Heatmap(
  mh_matrix,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  show_row_names = FALSE,
  column_names_gp = gpar(fontsize = 1),
  cell_fun = mh_col_fun,
  show_heatmap_legend = FALSE
)

pdf("/Users/cuijiajun/Desktop/others/tmphernan/2025_july_figHTFoantigen/step3_combine/mh_colorbar.pdf",
    width = ncol(mh_matrix) * 0.25 + 3, height = 1)
draw(mh_ht)
dev.off()

cat("✅ Done: Modern/Historical color bar generated.\n")



