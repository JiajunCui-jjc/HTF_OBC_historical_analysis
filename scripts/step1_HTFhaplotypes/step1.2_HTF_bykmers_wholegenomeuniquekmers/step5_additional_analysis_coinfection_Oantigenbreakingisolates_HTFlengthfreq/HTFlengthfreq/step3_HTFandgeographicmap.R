# ================================
# Historical HTF/O-antigen map
# ================================

library(sf)
library(dplyr)
library(ggplot2)
library(readr)
library(rnaturalearth)
library(rnaturalearthdata)

# === Paths ===
samples_path <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step5_additional_analysis_coinfection_Oantigenbreakingisolates_HTFlengthfreq/HTFlengthfreq/tmp489_withdatesandlocs_uniq.txt"
dominant_path <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step4_additionalanalysis/step3_HTFfreq/step2_3_HTF_oantigen_dominant_table.tsv"
out_plot <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step4_additionalanalysis/step3_HTFfreq/step3_h34_m57_OTU5candidates_map_HTFshape_yearfill.pdf"

# === Load and clean sample metadata ===
samples <- read.table(samples_path, header = TRUE, sep = '\t')
samples[samples$country == 'Congo', "country"] <- 'Democratic Republic of the Congo'
samples[samples$country == 'UK', "country"] <- 'United Kingdom'
# Add p1.G2 and p23.B2 manually
new_samples <- data.frame(
  samplename = c("p1.G2", "p23.B2"),
  country = "Germany",
  stringsAsFactors = FALSE
)

# Fill in missing columns with NA (matching the structure of 'samples')
for (col in setdiff(colnames(samples), colnames(new_samples))) {
  new_samples[[col]] <- NA
}

# Ensure column order matches
new_samples <- new_samples[, colnames(samples)]

# Append to samples
samples <- rbind(samples, new_samples)
#rm  p12.F2 p13.C7 p6.A10 p9.C4
# List of samples to remove
samples_to_remove <- c("p12.F2", "p13.C7", "p6.A10", "p9.C4")

# Filter out those samples
samples <- samples[!samples$samplename %in% samples_to_remove, ]

#match the sample names
samples$samplename <- ifelse(
  grepl("^p", samples$samplename),
  gsub("_", ".", samples$samplename),
  samples$samplename
)
# === Add HTF/O-antigen info and year ===
dominant_df <- read.table(dominant_path, header = TRUE, sep = "\t")
samples <- left_join(samples, dominant_df[, c("Isolate", "OAntigenStatus", "Year")],
                     by = c("samplename" = "Isolate")) %>%
  filter(!is.na(OAntigenStatus))

# === Load high-resolution world map ===
world <- ne_countries(scale = "medium", returnclass = "sf")

# === Assign static coordinates per country ===
# Use country centroids (so it's reproducible and accurate)
country_coords <- world %>%
  st_centroid(of_largest_polygon = TRUE) %>%
  st_coordinates() %>%
  as.data.frame()
country_coords$country <- world$name

samples <- left_join(samples, country_coords[, c("country", "X", "Y")], by = "country") %>%
  rename(lon = X, lat = Y)

# === Plot map ===
p <- ggplot() +
  geom_sf(data = world, fill = "#D7E8F4", color = "#CCE1EF", size = 0.1) +  # map base with thin border
  geom_point(data = samples,
             aes(x = lon, y = lat, shape = OAntigenStatus, fill = Year),
             position = position_jitter(width = 1.5, height = 1.5),  # gentle jitter
             size = 0.5, color = "grey30", stroke = 0.1) +
  scale_shape_manual(values = c(
    "HTF group with O-antigen" = 21,
    "HTF group without O-antigen" = 24
  )) +
  scale_fill_gradient(low = "darkgreen", high = "yellow", name = "Year", na.value = "grey90") +
  theme_minimal() +
  labs(
    title = "HTF Group Shape & Year Fill by Historical Sample Location",
    x = NULL, y = NULL, shape = "O-antigen Group"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )

# === Save ===
ggsave(out_plot, plot = p, width = 10, height = 6)
cat("✅ Saved map to:", out_plot, "\n")
