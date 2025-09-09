#!/usr/bin/env Rscript
  library(data.table)
  library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)



args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript script.R <task_id>")

## --- I/O ---
dump_dir <- "/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step2_isolate_kmers/isolate_jf/dump"
out_dir  <- file.path(dirname(dump_dir), "distribution")
out_dir_filtered <- file.path(out_dir, "filtered")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_filtered, showWarnings = FALSE, recursive = TRUE)

task_id <- as.integer(args[[1]])
if (is.na(task_id) || task_id < 1) stop("task_id must be a positive integer.")
dump_files <- sort(list.files(dump_dir, pattern = "_dump\\.txt$", full.names = TRUE))
if (length(dump_files) == 0) stop("No *_dump.txt in: ", dump_dir)
if (task_id > length(dump_files)) stop("task_id=", task_id, " > #files=", length(dump_files))

file <- dump_files[task_id]
sample <- sub("_dump\\.txt$", "", basename(file))

## --- Read ---
dt <- fread(file, header = FALSE, col.names = c("kmer","count"))
if (nrow(dt) == 0L) {
  message("Empty file: ", file)
  quit(status = 0)
}

## --- Core logic (matches awk) ---
x <- dt$count

# 1) Left-trim: keep only counts > 2
x_gt2 <- x[x > 2]
if (!length(x_gt2)) {
  message(sample, ": no counts > 2; skipping")
  quit(status = 0)
}

# 2) Stats on >2
m  <- mean(x_gt2)
sd <- sd(x_gt2)
if (is.na(sd)) sd <- 0
thr <- m + 3*sd

# 3) Trim right tail: keep <= threshold
x_kept <- x_gt2[x_gt2 <= thr]
if (!length(x_kept)) {
  message(sprintf("%s: filtered out (mean=%.3f sd=%.3f thr=%.3f)", sample, m, sd, thr))
  quit(status = 0)
}

# 4) Frequency table for plotting
freq <- as.data.table(table(x_kept))
setnames(freq, c("count","frequency"))
freq[, count := as.numeric(count)]

## --- Plot ---
p <- ggplot(freq, aes(x = count, y = frequency)) +
  geom_col(fill = "#0072B2", alpha = 0.95, width = 0.9) +
  labs(
    title    = sprintf("Filtered k-mer depth (>%s & ≤ mean+3·sd): %s", 2, sample),
    subtitle = sprintf("mean=%.2f, sd=%.2f, thr=%.2f; final mean≈%d",
                       m, sd, thr, round(mean(x_kept))),
    x = "k-mer depth", y = "Number of k-mers"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle    = element_text(hjust = 0.5, colour = "grey30")
  )

pdf_path <- file.path(out_dir_filtered, paste0(sample, "_filtered_kmer_depth_dist.pdf"))
ggsave(pdf_path, plot = p, width = 7, height = 5,
       device = grDevices::pdf, useDingbats = FALSE)

message("Saved: ", pdf_path)
