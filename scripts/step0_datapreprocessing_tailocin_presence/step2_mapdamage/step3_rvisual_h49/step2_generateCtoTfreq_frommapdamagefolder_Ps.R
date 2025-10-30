#!/usr/bin/env Rscript
# -------------------------------------------------------------------------
# Extract 5′ C→T and 3′ G→A misincorporation frequencies from all mapDamage folders
# Output two separate files:
#   - allinoneCtoT.txt  (5′ end)
#   - allinoneGtoA.txt  (3′ end)
# -------------------------------------------------------------------------

setwd("/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step2_mapdamageplot/tops")

suppressPackageStartupMessages({
  library(reshape2)
  library(dplyr)
})

# === Output files ===
outfile_CtoT <- "allinoneCtoT.txt"
outfile_GtoA <- "allinoneGtoA.txt"
file.remove(c(outfile_CtoT, outfile_GtoA))

# Create header for both files
header <- data.frame(
  pos = integer(),
  ref_base = character(),
  value = integer(),
  group = character(),
  
  total = integer(),
  proportion = numeric(),
  pos_from_end = integer(),
  sample = character()
)
#pos	ref_base	value	group	total	proportion	pos_from_end	sample
write.table(header, outfile_CtoT, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(header, outfile_GtoA, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# === Find all mapDamage folders ===
folders <- list.dirs(path = ".", recursive = FALSE)
#27.ESP_1975.mapped_to_Pseudomonas.dd.q20.markeddup.mapDamage
folders <- folders[grep("mapped_to_Pseudomonas.dd.q20.*mapDamage$", folders)]
cat("Found", length(folders), "mapDamage folders\n")

# === Main loop ===
for (indir in folders) {
  samplename <- basename(indir)
  cat("\n🔹 Processing:", samplename, "\n")
  
  infile <- file.path(indir, "misincorporation.txt")
  if (!file.exists(infile)) {
    cat("⚠️  Skipping (no misincorporation.txt):", indir, "\n")
    next
  }
  
  # ---- Load table ----
  tableAt <- tryCatch(read.table(infile, header = TRUE), error = function(e) NULL)
  if (is.null(tableAt)) {
    cat("⚠️  Error reading", infile, "\n")
    next
  }
  
  for (end_type in c("5p", "3p")) {
    
    subtable <- subset(tableAt, End == end_type)[, 3:22]
    bases <- c("A","C","G","T")
    total_sum <- aggregate(. ~ Pos, data = subtable[, c("Pos", bases)], sum)
    total_long <- reshape2::melt(total_sum, id.vars = "Pos",
                                 variable.name = "base", value.name = "total")
    
    mut_cols <- grep("\\.", colnames(subtable), value = TRUE)
    dfall <- do.call(rbind, lapply(mut_cols, function(col) {
      tmp <- aggregate(subtable[[col]] ~ subtable$Pos, data = subtable, sum)
      data.frame(pos = tmp[,1], value = tmp[,2], group = col)
    }))
    
    dfall$ref_base <- sub("^(.)\\..*$", "\\1", dfall$group)
    dfall <- merge(dfall, total_long, by.x = c("pos","ref_base"),
                   by.y = c("Pos","base"), all.x = TRUE)
    dfall$proportion <- dfall$value / dfall$total
    dfall$pos_from_end <- dfall$pos
    
    # ---- Select appropriate substitution ----
    if (end_type == "5p") {
      dfone <- subset(dfall, group == "C.T")
      dfone$sample <- samplename
      write.table(dfone, outfile_CtoT, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
      cat("✅ Done:", samplename, "→ 5′ C→T\n")
    } else {
      dfone <- subset(dfall, group == "G.A")
      dfone$sample <- samplename
      write.table(dfone, outfile_GtoA, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
      cat("✅ Done:", samplename, "→ 3′ G→A\n")
    }
  }
}

cat("\n🎯 Finished! Outputs saved to:\n  -", outfile_CtoT, "\n  -", outfile_GtoA, "\n")
