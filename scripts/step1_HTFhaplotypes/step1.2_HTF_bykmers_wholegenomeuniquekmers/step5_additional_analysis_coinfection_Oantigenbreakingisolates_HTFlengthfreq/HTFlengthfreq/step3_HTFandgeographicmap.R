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
# Create a copy so you don't modify the original
samples_plot <- samples

# Override lon/lat for the Russian sample
samples_plot$lon[samples_plot$country == "Russia"] <- 40   # move westward
samples_plot$lat[samples_plot$country == "Russia"] <- 60   # optionally tweak north/south
# === Plot: O-antigen group colored on map (manuscript ready) ===
p <- ggplot() +
  geom_sf(data = world, fill = "#F5F5F5", color = "#CCCCCC", size = 0.1) +
  geom_point(
    data = samples_plot,  # <-- use adjusted coordinates
    aes(x = lon, y = lat, color = OAntigenStatus),
    position = position_jitter(width = 0.7, height = 1),
    size = 1.5, alpha = 0.9
  ) +
  scale_color_manual(
    name = "O-antigen Group",
    values = c(
      "HTF group with O-antigen" = "#E69F00",
      "HTF group without O-antigen" = "#56B4E9"
    )
  ) +
  coord_sf(
    xlim = c(-15, 55), 
    ylim = c(30, 70), 
    expand = FALSE
  ) +
  theme_minimal(base_size = 14) +
  labs(x = NULL, y = NULL) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "top",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    plot.margin = margin(10, 10, 10, 10)
  )
# === Save ===
ggsave(out_plot, plot = p, width = 10, height = 6)
cat("✅ Saved map to:", out_plot, "\n")
