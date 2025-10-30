# ============================================================
# Safe k-mer mapping with reverse-complement check (31 bp span)
# ============================================================

# ---- Libraries ----
library(Biostrings)
library(ggplot2)
library(dplyr)
library(stringr)

# ---- Paths ----
fasta_file <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_April/keykmer_hHTFs/v2wholgenome_filterwithwholgenomekmerdepth_andhamming/step1_cal_hamming_foreachrefpair/runlocationofkmers/7HTF_dna.txt"
kmer_dir   <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_April/keykmer_hHTFs/v2wholgenome_filterwithwholgenomekmerdepth_andhamming/step1_cal_hamming_foreachrefpair/kmers_unique/"
filtered_kmer_dir <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_April/keykmer_hHTFs/v2wholgenome_filterwithwholgenomekmerdepth_andhamming/step1_cal_hamming_foreachrefpair/kmer_unique_diverse/"
out_dir    <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_April/keykmer_hHTFs/v2wholgenome_filterwithwholgenomekmerdepth_andhamming/step1_cal_hamming_foreachrefpair/runlocationofkmers/"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load reference sequences ----
htf_seqs <- readDNAStringSet(fasta_file)
names(htf_seqs) <- gsub("\\s.*", "", names(htf_seqs))
htf_lengths <- sapply(htf_seqs, length)
kmer_length <- 31L

# ---- Safe function to find k-mer positions ----
find_kmer_positions <- function(htf_name, kmer_file, status_label) {
  if (!file.exists(kmer_file)) return(NULL)
  
  # Read and clean k-mers
  kmers <- readLines(kmer_file)
  kmers <- toupper(trimws(kmers))
  kmers <- kmers[kmers != ""]
  
  seq <- htf_seqs[[htf_name]]
  
  positions <- sapply(kmers, function(k) {
    # Forward strand
    pos <- matchPattern(k, seq)
    if (length(pos) > 0) {
      start(pos)[1]
    } else {
      # Try reverse complement
      pos_rc <- matchPattern(as.character(reverseComplement(DNAString(k))), seq)
      if (length(pos_rc) > 0) start(pos_rc)[1] else NA
    }
  })
  
  data.frame(
    kmer = kmers,
    start = positions,
    end = positions + kmer_length - 1,
    HTF = htf_name,
    Status = status_label,
    stringsAsFactors = FALSE
  )
}

# ---- HTF names ----
htf_names <- names(htf_seqs)

# ---- Before filtering ----
df_before <- do.call(rbind, lapply(htf_names, function(htf) {
  kmer_file <- file.path(kmer_dir, paste0(htf, "_unique.txt"))
  find_kmer_positions(htf, kmer_file, "Before")
}))

# ---- After filtering ----
df_after <- do.call(rbind, lapply(htf_names, function(htf) {
  kmer_file <- file.path(filtered_kmer_dir, paste0(htf, "_unique_filtered_iterative2.txt"))
  find_kmer_positions(htf, kmer_file, "After")
}))

# ---- Count kmers (from file) ----
count_kmers_in_file <- function(htf, dir_path, suffix) {
  f <- file.path(dir_path, paste0(htf, suffix))
  if (file.exists(f)) length(readLines(f)) else 0
}
kmer_counts <- data.frame(
  HTF = htf_names,
  filtered_kmers = sapply(htf_names, count_kmers_in_file,
                          dir_path = filtered_kmer_dir,
                          suffix = "_unique_filtered_iterative2.txt"),
  stringsAsFactors = FALSE
)

# ---- Merge ----
df_combined <- bind_rows(df_before, df_after) %>%
  mutate(
    Length = htf_lengths[HTF],
    filtered_kmers = kmer_counts$filtered_kmers[match(HTF, kmer_counts$HTF)],
    HTF_label = ifelse(Status == "After",
                       paste0(HTF, " (", Length, "bp) [", filtered_kmers, " kmers]"),
                       paste0(HTF, " (", Length, "bp)")),
    Status = factor(Status, levels = c("Before", "After"))
  )

# ---- Colors ----
htf_levels <- c(
  "HTF_p7.G11 (1830bp)",
  "HTF_p25.A12 (1383bp)",
  "HTF_p25.C2 (1803bp)",
  "HTF_p5.D5 (1803bp)",
  "HTF_p26.D6 (1803bp)",
  "HTF_p23.B8 (1803bp)",
  "HTF_p21.F9 (1245bp)"
)
htf_colors <- setNames(
  c("#f6d6ff", "#638ccc", "#800233", "#800233", "#800233", "#800233", "#f9d42a"),
  htf_levels
)

# ---- Plot BEFORE ----
df_b <- df_combined %>% filter(Status == "Before" & !is.na(start))
p_b <- ggplot(df_b, aes(x = start, xend = end, y = HTF_label, yend = HTF_label)) +
  geom_segment(aes(color = HTF_label), alpha = 0.5, linewidth = 1.0) +
  scale_color_manual(values = htf_colors, guide = "none") +
  labs(
    title = "K-mer Positions Before Filtering (31 bp spans)",
    x = "Position along HTF reference",
    y = "HTF Haplotype"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 9)
  )

ggsave(file.path(out_dir, "HTF_kmer_positions_before.pdf"), p_b, width = 9, height = 5)

# ---- Plot AFTER ----
df_a <- df_combined %>% filter(Status == "After" & !is.na(start))
# color key stripping “[N kmers]”
color_key <- gsub(" \\[.*", "", df_a$HTF_label)
p_a <- ggplot(df_a, aes(x = start, xend = end, y = HTF_label, yend = HTF_label)) +
  geom_segment(aes(color = color_key), alpha = 0.85, linewidth = 1.2) +
  scale_color_manual(values = htf_colors, guide = "none") +
  labs(
    title = "K-mer Positions After Filtering (31 bp spans)",
    x = "Position along HTF reference",
    y = "HTF Haplotype"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 9)
  )

ggsave(file.path(out_dir, "HTF_kmer_positions_after.pdf"), p_a, width = 9, height = 5)

# ---- Combined (Before + After) ----
df_b$Panel <- "Before"
df_a$Panel <- "After"
df_comb_plot <- bind_rows(df_b, df_a)

p_combined <- ggplot(df_comb_plot,
                     aes(x = start, xend = end, y = HTF_label, yend = HTF_label)) +
  geom_segment(aes(color = gsub(" \\[.*", "", HTF_label)), alpha = 0.75, linewidth = 1.0) +
  scale_color_manual(values = htf_colors, guide = "none") +
  facet_wrap(~ Panel, ncol = 1, scales = "free_y") +
  labs(
    title = "K-mer Positions in HTF References (Before vs After Filtering)",
    x = "Position (31 bp span)",
    y = "HTF Haplotype"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "HTF_kmer_positions_combined.pdf"), p_combined, width = 10, height = 7)

# ---- Summary ----
cat("\n✅ Finished.\nSaved plots to:\n", out_dir, "\n\n")