# ============================
# Replace At/Ps columns in metadata using new dedup stats
# Compute summaries for all49 and hATUE5_43 (exclude six)
# ============================

library(dplyr)

# ---- Set working directory ----
setwd("/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step1_h49metadata/")

# ---- Full paths ----
meta_file <-  "tmp489_withdatesandlocs_uniq.txt"
newstats_file <- "h49dedup_maptoOTU5withouthaplotype_stats.txt"


outfile2 <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step1_h49metadata/SuppTable_h49.txt"
outfile3 <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/figureandtable/supptables/SuppTable_h49.txt"

summary_file2 <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step1_h49metadata/Summary_h49AtPs_stats.txt"

# ---- Load ----
meta <- read.table(meta_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
newstats <- read.table(newstats_file, header = TRUE, sep = " ", stringsAsFactors = FALSE, check.names = FALSE)

# ---- Harmonize column names ----

# ---- Drop old At/Ps columns ----
meta_clean <- meta %>%
  select(-At_percent, -At_covered, -At_depth,
         -Ps_percent, -Ps_covered, -Ps_depth,-group)
meta_clean$source[meta_clean$source == 'PL'] <- 'This study'
meta_clean$source[meta_clean$source == '(Latorre et al., 2022)'] <- '(Lang et al., 2024)'
meta_clean$source[meta_clean$source == '(Lopez et al., 2022)'] <- '(Lopez et al., 2025)'


# ---- Merge ----meta_updated <- meta_clean %>%
meta_updated <- meta_clean %>%
  right_join(newstats, by = "samplename") %>%
  rename(
    SAMPLE  = samplename,
    TYPE    = type,
    YEAR    = year,
    COUNTRY = country
  )


# ---- Desired order ----
desired_order <- c(
  "SAMPLE", "TYPE", "YEAR", "COUNTRY",
  "At_percent", "At_covered", "At_depth",
  "Ps_percent", "Ps_covered", "Ps_depth", "source"
)
available_cols <- intersect(desired_order, colnames(meta_updated))
meta_updated <- meta_updated[, c(available_cols, setdiff(colnames(meta_updated), available_cols))]

# ---- Check missing samples ----
missing <- setdiff(meta$SAMPLE, newstats$SAMPLE)
if (length(missing) > 0) {
  cat("⚠️ Warning: samples in metadata without updated stats:\n")
  print(missing)
}

# ---- Save updated metadata ----
write.table(meta_updated, outfile2, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(meta_updated, outfile3, sep = "\t", quote = FALSE, row.names = FALSE)
cat("✅ Updated metadata written to:\n", normalizePath(outfile), "\n", outfile2, "\n", outfile3, "\n")

# ============================
# Summary block: avg (min–max)
# ============================

summ_stat <- function(x) {
  if (all(is.na(x))) return("NA")
  sprintf("%.3f (%.3f–%.3f)", mean(x, na.rm = TRUE), min(x, na.rm = TRUE), max(x, na.rm = TRUE))
}

stat_cols <- c("At_percent", "At_covered", "At_depth",
               "Ps_percent", "Ps_covered", "Ps_depth")

# ---- Exclude six non-ATUE5 ----
excluded_samples <- c("HB0828", "HB0863", "PL0066", "PL0108", "PL0203", "PL0258")
meta_trim <- meta_updated %>% filter(!SAMPLE %in% excluded_samples)

# ---- Build summary text ----
summary_lines <- c(
  "========================================================",
  " Summary of At/Ps statistics (mean ± range)",
  "========================================================",
  ""
)

# --- All 49 samples ---
summary_lines <- c(summary_lines, ">> All 49 samples:")
for (col in stat_cols) {
  if (col %in% names(meta_updated)) {
    summary_lines <- c(summary_lines,
                       sprintf("%-12s: %s", col, summ_stat(meta_updated[[col]])))
  }
}
summary_lines <- c(summary_lines, "")

# --- hATUE5 43 samples (excluded six) ---
summary_lines <- c(summary_lines, ">> hATUE5 (43 samples, excluding HB0828, HB0863, PL0066, PL0108, PL0203, PL0258):")
for (col in stat_cols) {
  if (col %in% names(meta_trim)) {
    summary_lines <- c(summary_lines,
                       sprintf("%-12s: %s", col, summ_stat(meta_trim[[col]])))
  }
}
summary_lines <- c(summary_lines, "========================================================")

# ---- Write summary to file ----
writeLines(summary_lines, summary_file2)

cat("✅ Summary written to:\n", normalizePath(summary_file), "\n", summary_file2, "\n\n")

# ---- Also print on screen ----
cat(paste(summary_lines, collapse = "\n"), "\n")