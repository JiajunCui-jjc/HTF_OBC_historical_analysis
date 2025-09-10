# Load libraries
library(Biostrings)
library(ggplot2)
library(dplyr)

# Set paths
fasta_file <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_April/keykmer_hHTFs/v2wholgenome_filterwithwholgenomekmerdepth_andhamming/step1_cal_hamming_foreachrefpair/runlocationofkmers/7HTF_dna.txt"
kmer_dir <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_April/keykmer_hHTFs/v2wholgenome_filterwithwholgenomekmerdepth_andhamming/step1_cal_hamming_foreachrefpair/kmers_unique/"
filtered_kmer_dir <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_April/keykmer_hHTFs/v2wholgenome_filterwithwholgenomekmerdepth_andhamming/step1_cal_hamming_foreachrefpair/kmer_unique_diverse/"
out_dir <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_April/keykmer_hHTFs/v2wholgenome_filterwithwholgenomekmerdepth_andhamming/step1_cal_hamming_foreachrefpair/runlocationofkmers/"

# Load sequences and compute lengths
htf_seqs <- readDNAStringSet(fasta_file)
names(htf_seqs) <- gsub("\\s.*", "", names(htf_seqs))
htf_lengths <- sapply(htf_seqs, length)

# Function to find kmer positions
find_kmer_positions <- function(htf_name, kmer_file, status_label) {
  if (!file.exists(kmer_file)) return(NULL)
  kmers <- readLines(kmer_file)
  seq <- htf_seqs[[htf_name]]
  positions <- sapply(kmers, function(k) {
    pos <- matchPattern(k, seq)
    if (length(pos) > 0) start(pos)[1] else NA
  })
  df <- data.frame(
    kmer = kmers,
    position = positions,
    HTF = htf_name,
    Status = status_label,
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$position), ]
  return(df)
}

# Load HTF names
htf_names <- names(htf_seqs)

# Load k-mer positions
df_before <- do.call(rbind, lapply(htf_names, function(htf) {
  find_kmer_positions(htf, file.path(kmer_dir, paste0(htf, "_unique.txt")), "Before")
}))

df_after <- do.call(rbind, lapply(htf_names, function(htf) {
  find_kmer_positions(htf, file.path(filtered_kmer_dir, paste0(htf, "_unique_filtered_iterative2.txt")), "After")
}))

# Count filtered kmers for y-axis labeling
kmer_counts <- df_after %>%
  group_by(HTF) %>%
  summarise(filtered_kmers = n(), .groups = "drop")

# Merge and label
df_combined <- bind_rows(df_before, df_after) %>%
  mutate(
    Length = htf_lengths[HTF],
    filtered_kmers = kmer_counts$filtered_kmers[match(HTF, kmer_counts$HTF)],
    HTF_label = ifelse(Status == "After",
                       paste0(HTF, " (", Length, "bp) [", filtered_kmers, " kmers]"),
                       paste0(HTF, " (", Length, "bp)")),
    Status = factor(Status, levels = c("Before", "After"))
  )

# Define HTF order and colors
htf_levels <- c(
  "HTF_p7.G11 (1830bp)",
  "HTF_p25.A12 (1383bp)",
  "HTF_p25.C2 (1803bp)",
  "HTF_p5.D5 (1803bp)",
  "HTF_p26.D6 (1803bp)",
  "HTF_p23.B8 (1803bp)",
  "HTF_p21.F9 (1245bp)"
)
htf_colors <- structure(
  c("#f6d6ff", "#638ccc", "#800233", "#800233", "#800233", "#800233", "#f9d42a"),
  names = htf_levels
)

# Set plot factor levels
df_combined$HTF_label <- factor(df_combined$HTF_label,
                                levels = unique(df_combined$HTF_label[
                                  order(df_combined$HTF, df_combined$Status)
                                ]))

# === PLOT 1: BEFORE ===
df_b <- df_combined %>% filter(Status == "Before")
p_b <- ggplot(df_b, aes(x = position, y = HTF_label)) +
  geom_point(size = 1.2, alpha = 0.5, color = htf_colors[as.character(df_b$HTF_label)]) +
  labs(
    title = "K-mer Positions Before Filtering",
    x = "Position", y = "HTF Haplotype"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.y = element_text(size = 9))

ggsave(file.path(out_dir, "HTF_kmer_positions_before.pdf"), p_b, width = 9, height = 5)

# === PLOT 2: AFTER ===
df_a <- df_combined %>% filter(Status == "After")
p_a <- ggplot(df_a, aes(x = position, y = HTF_label)) +
  geom_point(size = 1.4, alpha = 0.85, color = htf_colors[as.character(gsub(" \\[.*", "", df_a$HTF_label))]) +
  labs(
    title = "K-mer Positions After Filtering",
    x = "Position", y = "HTF Haplotype"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(size = 10, hjust = 0.5),
        axis.text.y = element_text(size = 9))

ggsave(file.path(out_dir, "HTF_kmer_positions_after.pdf"), p_a, width = 9, height = 5)


#combine
# Load libraries
library(Biostrings)
library(ggplot2)
library(dplyr)

# Set paths
fasta_file <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_April/keykmer_hHTFs/v2wholgenome_filterwithwholgenomekmerdepth_andhamming/step1_cal_hamming_foreachrefpair/runlocationofkmers/7HTF_dna.txt"
kmer_dir <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_April/keykmer_hHTFs/v2wholgenome_filterwithwholgenomekmerdepth_andhamming/step1_cal_hamming_foreachrefpair/kmers_unique/"
filtered_kmer_dir <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_April/keykmer_hHTFs/v2wholgenome_filterwithwholgenomekmerdepth_andhamming/step1_cal_hamming_foreachrefpair/kmer_unique_diverse/"
out_dir <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_April/keykmer_hHTFs/v2wholgenome_filterwithwholgenomekmerdepth_andhamming/step1_cal_hamming_foreachrefpair/runlocationofkmers/"

# Load sequences
htf_seqs <- readDNAStringSet(fasta_file)
names(htf_seqs) <- gsub("\\s.*", "", names(htf_seqs))

# Function to load kmers and find positions
find_kmer_positions <- function(htf_name, kmer_file, status_label) {
  kmers <- readLines(kmer_file)
  seq <- htf_seqs[[htf_name]]
  positions <- sapply(kmers, function(k) {
    pos <- matchPattern(k, seq)
    if (length(pos) > 0) start(pos)[1] else NA
  })
  data.frame(
    kmer = kmers,
    position = positions,
    HTF = htf_name,
    Status = status_label
  ) |> filter(!is.na(position))
}

# Get HTF names
htf_names <- names(htf_seqs)

# Original kmers
df_all <- do.call(rbind, lapply(htf_names, function(htf) {
  kmer_file <- file.path(kmer_dir, paste0(htf, "_unique.txt"))
  if (file.exists(kmer_file)) {
    find_kmer_positions(htf, kmer_file, status_label = "beforehamming")
  }
}))

# Filtered kmers
df_filtered <- do.call(rbind, lapply(htf_names, function(htf) {
  kmer_file <- file.path(filtered_kmer_dir, paste0(htf, "_unique_filtered_iterative2.txt"))
  if (file.exists(kmer_file)) {
    find_kmer_positions(htf, kmer_file, status_label = "afterhamming")
  }
}))

# Combine with new labels for y-axis
df_all$HTFstatus <- paste0(df_all$HTF, "_beforehamming")
df_filtered$HTFstatus <- paste0(df_filtered$HTF, "_afterhamming")

df_combined <- bind_rows(df_all, df_filtered)

# Plot combined
p_combined <- ggplot(df_combined, aes(x = position, y = HTFstatus, color = Status)) +
  geom_point(alpha = 0.7, size = 1.1) +
  scale_color_manual(values = c("beforehamming" = "darkblue", "afterhamming" = "darkred")) +
  theme_minimal() +
  labs(
    title = "K-mer Position in HTF References (Before vs After Filtering)",
    x = "Position",
    y = "HTF (Status)",
    color = "K-mer Set"
  )

# Save
ggsave(file.path(out_dir, "HTF_kmer_positions_combined.pdf"), p_combined, width = 10, height = 6)

cat("✅ Combined plot saved to:", out_dir, "\n")