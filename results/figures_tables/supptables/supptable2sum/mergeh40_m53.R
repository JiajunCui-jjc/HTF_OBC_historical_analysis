# === Input paths ===
setwd('/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/figures_tables/supptables/supptable2sum')
# === Input paths ===
cov_file <- "tailocin_coverage_summary.tsv"
bykmer_file <- "combined_HTF_oantigen_m57_h35.txt"
bylocal_file <- "final_h40_m57_HTFhaplotype_bylocalassembly.txt"

# === Read tables ===
cov <- read.delim(cov_file, header = TRUE, stringsAsFactors = FALSE)
bykmer <- read.delim(bykmer_file, header = TRUE, stringsAsFactors = FALSE)
bylocal <- read.delim(bylocal_file, header = TRUE, stringsAsFactors = FALSE)

# === Drop 'note' column from bykmer if it exists ===
bykmer$note <- NULL

# === Rename for consistency ===
names(bykmer)[names(bykmer) == "HTF_haplotype"] <- "HTF_haplotype_bykmer"
names(bylocal)[names(bylocal) == "top_HTFhaplotype"] <- "HTF_haplotype_bylocalassembly"

# === Merge stepwise ===
merged <- merge(cov, bykmer, by = "sample", all.x = TRUE)
merged <- merge(merged, bylocal, by = "sample", all.x = TRUE)

# === Remove specific samples ===
exclude <- c("p12.F2", "p13.C7", "p6.A10", "p9.C4")
merged <- merged[!(merged$sample %in% exclude), ]

# === Reorder columns ===
final_cols <- c(
  "sample", "tailocin_cov_prop", "tailocin_depth",
  "HTF_haplotype_bykmer", "HTF_haplotype_bylocalassembly",
  "HTFgroup_Oantigen_PA", "wfgD_PA", "rmlC1_PA",
  "tagG1_PA", "tagH1_PA", "spsA_PA", "espE2_PA", "espE2length"
)

merged_final <- merged[, final_cols]


# === Manual corrections to specific samples ===

# 1️⃣ Mark espE2-related columns (espE2_PA, espE2length) as NA for four isolates
na_espE2 <- c("p25.C11", "PL0139", "86.NOR_1911_S7", "27.ESP_1975")
merged_final[merged_final$sample %in% na_espE2, c("espE2_PA", "espE2length")] <- NA

# 2️⃣ Mark HTF-related columns as NA for 64.GBR_1933b_S36
merged_final[merged_final$sample == "64.GBR_1933b_S36",
             c("HTF_haplotype_bykmer", "HTF_haplotype_bylocalassembly", "HTFgroup_Oantigen_PA")] <- NA

# === Write output ===
write.table(merged_final, file = "merged_final_HTF_tailocin.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("[✔] Output saved to: merged_final_HTF_tailocin.txt\n")
cat("Rows:", nrow(merged_final), "\n")


