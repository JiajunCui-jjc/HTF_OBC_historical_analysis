library(purrr)
library(dplyr)
library(ggplot2)
library(readr)

read_lgdist <- function(file, sample, genome){
  dat <- suppressWarnings(read_tsv(file, comment="#", show_col_types = FALSE))
  
  # Handle different formats
  if(all(c("Std","Length","Occurences") %in% colnames(dat))){
    dat <- dat %>% select(Length, Occurences)
  } else if(all(c("Length","Occurences") %in% colnames(dat))){
    dat <- dat %>% select(Length, Occurences)
  } else {
    dat <- dat[, (ncol(dat)-1):ncol(dat)]
    colnames(dat) <- c("Length","Occurences")
  }
  
  # Force numeric + combine over strands if needed
  dat <- dat %>%
    mutate(
      Length = suppressWarnings(as.numeric(Length)),
      Occurences = suppressWarnings(as.numeric(Occurences))
    ) %>%
    group_by(Length) %>%
    summarise(Count = sum(Occurences, na.rm=TRUE), .groups="drop") %>%
    mutate(
      Sample = sample,
      Genome = genome
    )
  
  return(dat)
}

# === directories ===
at_dirs <- c(
  "/Users/cuijiajun/Desktop/2023-2024 PhD ucl/2024_aMeta/wholepipeAt_Ps/2024_233_analysis/mapdamage/mapdamage2024_29/toAt",
  "/Users/cuijiajun/Desktop/2023-2024 PhD ucl/2024_aMeta/wholepipeAt_Ps/2024_233_analysis/mapdamage/mapdamage2024_22/for17Hold/usedmaptoAt"
)

ps_dirs <- c(
  "/Users/cuijiajun/Desktop/2023-2024 PhD ucl/2024_aMeta/wholepipeAt_Ps/2024_233_analysis/mapdamage/mapdamage2024_29/tops",
  "/Users/cuijiajun/Desktop/2023-2024 PhD ucl/2024_aMeta/wholepipeAt_Ps/2024_233_analysis/mapdamage/mapdamage2024_22/for17Hold/usedforCFML_ps"
)

# === collect files ===
at_files <- unlist(lapply(at_dirs, function(d)
  list.files(d, pattern="(lgdistribution|length_distribution)\\.txt$", recursive=TRUE, full.names=TRUE)
))
ps_files <- unlist(lapply(ps_dirs, function(d)
  list.files(d, pattern="(lgdistribution|length_distribution)\\.txt$", recursive=TRUE, full.names=TRUE)
))

# === read ===
at_all <- map_df(at_files, ~read_lgdist(.x, basename(dirname(.x)), "Arabidopsis thaliana"))
ps_all <- map_df(ps_files, ~read_lgdist(.x, basename(dirname(.x)), "Pseudomonas viridiflava"))


## Keep only HB* and PL* samples in both datasets
at_all <- at_all %>%
  filter(grepl("(HB|PL)", Sample))

ps_all <- ps_all %>%
  filter(grepl("(HB|PL)", Sample))
unique(at_all$Sample)
#29 PL + 7 HB no 10 lopez


# === convert counts to proportions per sample ===
at_all <- at_all %>%
  group_by(Sample) %>%
  mutate(Prop = Count / sum(Count, na.rm=TRUE)) %>%
  ungroup()

ps_all <- ps_all %>%
  group_by(Sample) %>%
  mutate(Prop = Count / sum(Count, na.rm=TRUE)) %>%
  ungroup()

# === plots ===
library(scales)   # for percent_format

# Arabidopsis thaliana
p_at <- ggplot(at_all %>% filter(Length > 0),
               aes(x=Length, y=Prop, group=Sample)) +
  geom_line(alpha=0.5, color="#4575b4") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +   # <-- convert to %
  labs(x="Fragment length (bp)", y="Reads (%)",
       title=expression(paste("Fragment length distribution - ", italic("Arabidopsis thaliana")))) +
  theme_bw(base_size=14)

# Pseudomonas viridiflava
p_ps <- ggplot(ps_all %>% filter(Length > 0),
               aes(x=Length, y=Prop, group=Sample)) +
  geom_line(alpha=0.5, color="#d73027") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +   # <-- convert to %
  labs(x="Fragment length (bp)", y="Reads (%)",
       title=expression(paste("Fragment length distribution - ", italic("Pseudomonas "),'sp.'))) +
  theme_bw(base_size=14)

# Save

out_dir <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step2_mapdamageplot/step5_lengthdistribution_insertsize"
ggsave(file.path(out_dir, "PLHBFragmentLength_Athaliana_proportion.pdf"), p_at, width=7, height=5)
ggsave(file.path(out_dir, "PLHBFragmentLength_Ps_proportion.pdf"), p_ps, width=7, height=5)



#checked, 10 lopez sample
# === Find dominant peak (max Prop) for each sample ===
#ps_top12 <- ps_all %>%
##  group_by(Sample) %>%
#  slice_max(order_by = Prop, n = 1, with_ties = FALSE) %>%
#  ungroup() %>%
#  arrange(desc(Prop)) %>%
#  slice_head(n = 12)

#print(ps_top12)
#out_file <- file.path(out_dir, "Pseudomonas_Top12_Peaks.tsv")
#write.table(ps_top12, out_file, sep = "\t", quote = FALSE, row.names = FALSE)

# ============================
# Mixed plot: Arabidopsis + Pseudomonas
# ============================

library(patchwork)
library(scales)

# ============================
# Clean Nature-style Supp Figure
# ============================

# Arabidopsis thaliana panel
p_at <- ggplot(at_all %>% filter(Length > 0 ),
               aes(x = Length, y = Prop, group = Sample)) +
  geom_line(alpha = 0.6, color = "#4575b4") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 0.15)) +  # set same scale across panels
  labs(x = "Fragment length (bp)", y = "Reads (%)",
       title = expression(italic("Arabidopsis thaliana"))) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11, face = "bold"),
    axis.line = element_line(size = 0.6),
    axis.ticks = element_line(size = 0.6)
  )

# Pseudomonas viridiflava panel
p_ps <- ggplot(ps_all %>% filter(Length > 0),
               aes(x = Length, y = Prop, group = Sample)) +
  geom_line(alpha = 0.6, color = "#d73027") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 0.15)) +  # match scale
  labs(x = "Fragment length (bp)", y = "Reads (%)",
       title = expression(italic(Pseudomonas) ~ "sp.")) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11, face = "bold"),
    axis.line = element_line(size = 0.6),
    axis.ticks = element_line(size = 0.6)
  )

# Combine into two clean panels
p_combined <- p_at + p_ps +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = 'A',
                  tag_prefix = "", 
                  tag_suffix = "",
                  theme = theme(plot.tag = element_text(size = 14, face = "bold")))

# Save combined supplementary figure
ggsave(file.path(out_dir, "SuppFig_FragmentLength_AT_PS.pdf"),
       plot = p_combined, width = 12, height = 5, dpi = 600)

