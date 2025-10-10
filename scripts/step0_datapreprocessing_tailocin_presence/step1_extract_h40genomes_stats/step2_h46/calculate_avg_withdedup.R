# ============================
# Merge dedup stats with metadata & summarize including read proportions
# ============================

library(dplyr)

# ---- Paths ----
setwd("/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step1_h40metadata/")

meta_file <- "allh46_metadata.txt"   # before dedup
dedup_file <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/scripts/step0_datapreprocessing_tailocin_presence/step1_extract_h40genomes_stats/step2_h46/dedupstats_At_Ps/allh46_markdup_covdepth_summary.txt"

outfile_merged  <- "allh46_metadata_withdedup.txt"
outfile_summary <- "Summary_AtPs_BeforeAfterDedup.txt"

# ---- Load ----
meta  <- read.table(meta_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
dedup <- read.table(dedup_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# ---- Merge ----
merged <- meta %>%
  left_join(dedup, by = "samplename")

# ---- Save merged metadata ----
write.table(merged, outfile_merged, sep = "\t", quote = FALSE, row.names = FALSE)
cat("✅ Merged table written to:", outfile_merged, "\n")

# ---- Helper for summary ----
summ_stat <- function(x) {
  if (all(is.na(x))) return("NA")
  sprintf("%.3f (%.3f–%.3f)", mean(x, na.rm = TRUE), min(x, na.rm = TRUE), max(x, na.rm = TRUE))
}

# ---- Exclude 6 non-ATUE5 ----
excluded_samples <- c("HB0828", "HB0863", "PL0066", "PL0108", "PL0203", "PL0258")

# ---- Summary function ----
calc_summary <- function(df) {
  list(
    At_percent = summ_stat(df$At_percent),
    Ps_percent = summ_stat(df$Ps_percent),
    At_cov_pre = summ_stat(df$At_covered),
    At_cov_post = summ_stat(df$At_cov_markdup),
    Ps_cov_pre = summ_stat(df$Ps_covered),
    Ps_cov_post = summ_stat(df$Ps_cov_markdup),
    At_depth_pre = summ_stat(df$At_depth),
    At_depth_post = summ_stat(df$At_depth_markdup),
    Ps_depth_pre = summ_stat(df$Ps_depth),
    Ps_depth_post = summ_stat(df$Ps_depth_markdup)
  )
}

# ---- All 46 ----
stats_all <- calc_summary(merged)

# ---- Trimmed 40 ----
merged_trim <- merged %>% filter(!samplename %in% excluded_samples)
stats_trim  <- calc_summary(merged_trim)

# ---- Write summary ----
sink(outfile_summary)
cat("========================================================\n")
cat(" Summary of Arabidopsis & Pseudomonas Before vs After Deduplication\n")
cat("========================================================\n\n")

cat(">> All 46 samples:\n")
cat(sprintf("At_percent (read proportion): %s\n", stats_all$At_percent))
cat(sprintf("Ps_percent (read proportion): %s\n\n", stats_all$Ps_percent))

cat(sprintf("At_covered:  %s → %s\n", stats_all$At_cov_pre, stats_all$At_cov_post))
cat(sprintf("At_depth:    %s → %s\n", stats_all$At_depth_pre, stats_all$At_depth_post))
cat(sprintf("Ps_covered:  %s → %s\n", stats_all$Ps_cov_pre, stats_all$Ps_cov_post))
cat(sprintf("Ps_depth:    %s → %s\n\n", stats_all$Ps_depth_pre, stats_all$Ps_depth_post))

cat(">> hATUE5 (40 samples, excluding HB0828, HB0863, PL0066, PL0108, PL0203, PL0258):\n")
cat(sprintf("At_percent (read proportion): %s\n", stats_trim$At_percent))
cat(sprintf("Ps_percent (read proportion): %s\n\n", stats_trim$Ps_percent))

cat(sprintf("At_covered:  %s → %s\n", stats_trim$At_cov_pre, stats_trim$At_cov_post))
cat(sprintf("At_depth:    %s → %s\n", stats_trim$At_depth_pre, stats_trim$At_depth_post))
cat(sprintf("Ps_covered:  %s → %s\n", stats_trim$Ps_cov_pre, stats_trim$Ps_cov_post))
cat(sprintf("Ps_depth:    %s → %s\n", stats_trim$Ps_depth_pre, stats_trim$Ps_depth_post))
sink()

cat("✅ Summary saved to:", outfile_summary, "\n")