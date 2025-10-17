# ============================
# Replace At/Ps columns in metadata using new dedup stats (safe reorder + summary to file)
# ============================

library(dplyr)

# ---- Full paths ----
meta_file <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/scripts/step0_datapreprocessing_tailocin_presence/step1_extract_h40genomes_stats/step2_h46/tmph49save/old/SuppTable_h46_dedup_ready.txt"
newstats_file <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/scripts/step0_datapreprocessing_tailocin_presence/step1_extract_h40genomes_stats/step2_h46/tmph49save/h46afterdedup_newlopez10.txt"

outfile  <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/scripts/step0_datapreprocessing_tailocin_presence/step1_extract_h40genomes_stats/step2_h46/tmph49save/SuppTable_h46_dedup_updated.txt"
outfile2 <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step1_h40metadata/SuppTable_h46_dedup_updated.txt"
outfile3<- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/figures_tables/supptables/SuppTable_h46_dedup_updated.txt"

summary_file <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/scripts/step0_datapreprocessing_tailocin_presence/step1_extract_h40genomes_stats/step2_h46/tmph49save/Summary_AtPs_stats.txt"
summary_file2 <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step1_h40metadata/Summary_AtPs_stats.txt"

# ---- Load ----
meta <- read.table(meta_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
newstats <- read.table(newstats_file, header = TRUE, sep = " ", stringsAsFactors = FALSE, check.names = FALSE)

# ---- Harmonize column names ----
colnames(newstats)[1] <- "SAMPLE"

# ---- Drop old At/Ps columns ----
meta_clean <- meta %>%
  select(-At_percent, -At_covered, -At_depth,
         -Ps_percent, -Ps_covered, -Ps_depth)

# ---- Merge ----
meta_updated <- meta_clean %>%
  left_join(newstats, by = "SAMPLE")

# ---- Desired order ----
desired_order <- c(
  "SAMPLE", "TYPE", "YEAR", "COUNTRY",
  "At_percent", "At_covered", "At_depth",
  "Ps_percent", "Ps_covered", "Ps_depth",
  "group", "source"
)

# ---- Reorder safely ----
available_cols <- intersect(desired_order, colnames(meta_updated))
meta_updated <- meta_updated[, c(available_cols, setdiff(colnames(meta_updated), available_cols))]

# ---- Check missing samples ----
missing <- setdiff(meta$SAMPLE, newstats$SAMPLE)
if (length(missing) > 0) {
  cat("⚠️ Warning: samples in metadata without updated stats:\n")
  print(missing)
}

# ---- Save updated metadata ----
write.table(meta_updated, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(meta_updated, outfile2, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(meta_updated, outfile3, sep = "\t", quote = FALSE, row.names = FALSE)

cat("✅ Updated metadata written to:\n", outfile, "\n", outfile2, "\n")

# ============================
# Summary block: avg (min–max)
# ============================
summ_stat <- function(x) {
  if (all(is.na(x))) return("NA")
  sprintf("%.3f (%.3f–%.3f)", mean(x, na.rm = TRUE), min(x, na.rm = TRUE), max(x, na.rm = TRUE))
}

stat_cols <- c("At_percent", "At_covered", "At_depth",
               "Ps_percent", "Ps_covered", "Ps_depth")

# ---- Build summary text ----
summary_lines <- c(
  "========================================================",
  " Summary of At/Ps statistics (mean ± range)",
  "========================================================",
  ""
)

for (col in stat_cols) {
  if (col %in% names(meta_updated)) {
    summary_lines <- c(summary_lines,
                       sprintf("%-12s: %s", col, summ_stat(meta_updated[[col]])))
  }
}
summary_lines <- c(summary_lines, "========================================================")

# ---- Write summary to file ----
writeLines(summary_lines, summary_file)
writeLines(summary_lines, summary_file2)

cat("✅ Summary written to:", summary_file, "\n\n")

# ---- Also print on screen ----
cat(paste(summary_lines, collapse = "\n"), "\n")