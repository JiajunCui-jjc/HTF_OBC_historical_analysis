library(dplyr)

# === Set working directory ===
setwd("/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step1_cal_hamming_pairwise_matrix/")
#mkdir -p /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/kmers_unique_hamming2/
# === Load data ===
df <- read.table("pairwise_hamming_distance_longformatbypy.tsv", header=TRUE, sep="\t")


#iterative filtering from smallest largest to  to rm any <2 distance kmers and keep as much as possible in small HTF kmer pool like p26.d6 and p23.b8
# ✅ manually define the order:
ref_order <- c("HTF_p26.D6","HTF_p23.B8", "HTF_p21.F9", "HTF_p25.A12", "HTF_p5.D5", "HTF_p7.G11", "HTF_p25.C2")
cat("Processing order (manual):", ref_order, "\n")

# === 2️⃣ Initialize storage for kept kmers
kept_kmers <- list()

# === 3️⃣ Process each ref in order
for (ref in ref_order) {
  cat("\n️ Processing", ref, "\n")
  
  # all kmers in this ref
  kmers_this_ref <- unique(df$QueryKmer[df$QueryRef == ref])
  
  if (length(kept_kmers) == 0) {
    # first ref → no filtering needed
    kept_kmers[[ref]] <- kmers_this_ref
    cat("kept all", length(kmers_this_ref), "kmers from", ref, "\n")
  } else {
    # cumulative kept kmers so far
    cumulative_kept <- unlist(kept_kmers)
    
    # find conflicts: kmers with dist <2 to any cumulative kept kmer
    conflicts <- df %>%
      filter(QueryRef == ref, TargetKmer %in% cumulative_kept, HammingDist < 2) %>%
      pull(QueryKmer) %>%
      unique()
    
    # exclude conflicts
    filtered_kmers <- setdiff(kmers_this_ref, conflicts)
    
    kept_kmers[[ref]] <- filtered_kmers
    
    cat(" kept", length(filtered_kmers), "kmers from", ref,
        "(removed", length(conflicts), "conflicting kmers)\n")
  }
}

# === 4️⃣ Write filtered kmers to files
for (ref in names(kept_kmers)) {
  outfile <- paste0("../kmers_unique_hamming2/", ref, "_unique_filtered_iterative2.txt")
  writeLines(kept_kmers[[ref]], outfile)
  cat(" saved", length(kept_kmers[[ref]]), "kmers to", outfile, "\n")
}

cat("\n Iterative filtering (manual order) complete!\n")




