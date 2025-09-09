
# === Load libraries ===
library(ggplot2)

# === Set working directory and input paths ===
pswd <- "/Users/cuijiajun/Desktop/2023-2024 PhD ucl/2024_aMeta/wholepipeAt_Ps/2024_233_analysis/mapdamage/mapdamage2024_29/tops/"
setwd(pswd)

# Load first dataset
tablePs <- read.table(paste0(pswd, "allinoneCtoT.txt"), header = FALSE)
colnames(tablePs) <- c("pos", "value", "group", "Pos", "total", "base", "proportion", "pos_from_5end", "sample")

# Load second dataset
tablePs_17 <- read.table(
  "/Users/cuijiajun/Desktop/2023-2024 PhD ucl/2024_aMeta/wholepipeAt_Ps/2024_233_analysis/mapdamage/mapdamage2024_22/for17Hold/usedforCFML_ps/allinoneCtoT.txt",
  header = FALSE
)
colnames(tablePs_17) <- colnames(tablePs)

# Merge datasets
tablePs_all <- rbind(tablePs, tablePs_17)

# Read 46 historical sample names
samples_46 <- read.table(
  "/Users/cuijiajun/Desktop/others/tmphernan/2025_Aug_balgiumtalk_materials_jiajun/material1_mapdamage_h34/v1_h46/step1_h46ps/samples_46.txt",
  header = FALSE
)[[1]]

# Subset for the 46 historical samples
library(stringr)


# Combine into a regex OR pattern
pattern_46 <- paste(samples_46, collapse = "|")

# Filter rows where `sample` column contains any of the 46 substrings
tablePs_46 <- tablePs_all[str_detect(tablePs_all$sample, pattern_46), ]
# Save filtered data
write.table(
  tablePs_46,
  "/Users/cuijiajun/Desktop/others/tmphernan/2025_Aug_balgiumtalk_materials_jiajun/material1_mapdamage_h34/v1_h46/step1_h46ps/h46_CtoTfreq_Ps.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# === Plot: C-to-T damage frequency from 5' end ===
p <- ggplot(tablePs_46, aes(x = pos_from_5end, y = proportion, group = sample)) +
  geom_line(color = "#D04F4F", size = 0.8, alpha = 0.85) +
  labs(
    title = expression("5′ C to T deamination in 46 historical plant " * italic("Pseudomonas viridiflava") * " ATUE5 genomes"),
    x = "Distance from 5′ end (bp)",
    y = "C to T Frequency"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  ylim(0, 0.04)

# Save high-resolution figure
ggsave(
  filename = "/Users/cuijiajun/Desktop/others/tmphernan/2025_Aug_balgiumtalk_materials_jiajun/material1_mapdamage_h34/v1_h46/step1_h46ps/topsallinone_46.png",
  plot = p,
  width = 10,
  height = 6,
  dpi = 600
)
sort(unique(tablePs_46$sample))