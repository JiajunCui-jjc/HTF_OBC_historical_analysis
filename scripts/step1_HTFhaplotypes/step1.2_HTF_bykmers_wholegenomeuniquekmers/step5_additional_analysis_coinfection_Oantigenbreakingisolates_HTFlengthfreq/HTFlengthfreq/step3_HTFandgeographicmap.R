# ================================
# Historical + Modern HTF haplotypes map
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
out_plot <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step4_additionalanalysis/step3_HTFfreq/step3_h34_m57_OTU5candidates_map_HTFhaplotypes_modernTransparent.pdf"

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

for (col in setdiff(colnames(samples), colnames(new_samples))) {
  new_samples[[col]] <- NA
}
new_samples <- new_samples[, colnames(samples)]
samples <- rbind(samples, new_samples)

# Remove unwanted samples
samples_to_remove <- c("p12.F2", "p13.C7", "p6.A10", "p9.C4")
samples <- samples[!samples$samplename %in% samples_to_remove, ]

# Match naming
samples$samplename <- ifelse(
  grepl("^p", samples$samplename),
  gsub("_", ".", samples$samplename),
  samples$samplename
)

# === Add haplotype info ===
dominant_df <- read.table(dominant_path, header = TRUE, sep = "\t")
dominant_df <- dominant_df %>%
  mutate(LengthGroup = case_when(
    Reference == "HTF_p7.G11" ~ "1830",
    Reference == "HTF_p25.A12" ~ "1383",
    Reference %in% c("HTF_p5.D5", "HTF_p23.B8", "HTF_p26.D6", "HTF_p25.C2") ~ "1803",
    Reference == "HTF_p21.F9" ~ "1245",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(LengthGroup))

samples <- left_join(samples,
                     dominant_df[, c("Isolate", "LengthGroup", "Year", "IsolateType")],
                     by = c("samplename" = "Isolate")) %>%
  filter(!is.na(LengthGroup))

# === Map data ===
world <- ne_countries(scale = "medium", returnclass = "sf")

# Country centroids
country_coords <- world %>%
  st_centroid(of_largest_polygon = TRUE) %>%
  st_coordinates() %>%
  as.data.frame()
country_coords$country <- world$name

samples <- left_join(samples, country_coords[, c("country", "X", "Y")], by = "country") %>%
  rename(lon = X, lat = Y)
set.seed(19)  # reproducible jitter

# === Adjust lon/lat centers ===
samples_plot <- samples

# Russia center
samples_plot$lon[samples_plot$country == "Russia"] <- 32
samples_plot$lat[samples_plot$country == "Russia"] <- 60

# Norway center
samples_plot$lon[samples_plot$country == "Norway"] <- 10
samples_plot$lat[samples_plot$country == "Norway"] <- 62

# Sweden center
samples_plot$lon[samples_plot$country == "Sweden"] <- 15
samples_plot$lat[samples_plot$country == "Sweden"] <- 60

# Lithuania center (tight jitter)
samples_plot$lon[samples_plot$country == "Lithuania"] <- 25
samples_plot$lat[samples_plot$country == "Lithuania"] <- 55

# UK, Spain, Germany centers (will use wide jitter)



# === Plot with country-specific jitter ===
p <- ggplot() +
  geom_sf(data = world, fill = "#F5F5F5", color = "#CCCCCC", size = 0.1) +
  
  # Historical + Modern for UK, Spain, Germany (wider jitter)
  # ---- Modern FIRST (all Germany, faded, goes behind) ----
  geom_point(
  data = filter(samples_plot, country == "Germany" & IsolateType == "Modern"),
  aes(x = lon, y = lat, color = LengthGroup),
  position = position_jitter(width = 2, height = 1.5),  # wider spread for Germany
  size = 2, alpha = 0.1
  ) +
  
  # ---- Historical SECOND (bright, on top) ----
  # Spain + Germany (wide jitter)
  geom_point(
  data = filter(samples_plot, country %in% c("Spain", "Germany") & IsolateType == "Historical"),
  aes(x = lon, y = lat, color = LengthGroup),
  position = position_jitter(width = 2, height = 1.5),
  size = 2, alpha = 1
  ) +
  
  
  # Lithuania (tight jitter)
  geom_point(
    data = filter(samples_plot, country == "Lithuania"),
    aes(x = lon, y = lat, color = LengthGroup, alpha = IsolateType),
    position = position_jitter(width = 0.6, height = 0.6),
    size = 2
  ) +
  
  # Other countries (default jitter)
  geom_point(
    data = filter(samples_plot, !country %in% c( "Spain", "Germany", "Lithuania")),
    aes(x = lon, y = lat, color = LengthGroup, alpha = IsolateType),
    position = position_jitter(width = 0.8, height = 1.5),
    size = 2
  ) +
  
  scale_color_manual(values = length_colors) +
  scale_alpha_manual(values = c("Historical" = 1, "Modern" = 0.1)) +
  coord_sf(xlim = c(-15, 55), ylim = c(30, 70), expand = FALSE) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "top"
  )

# Save
ggsave(out_plot, plot = p, width = 10, height = 6)
cat("✅ Saved inside-polygon haplotype map to:", out_plot, "\n")




#the legend
library(ggplot2)
library(dplyr)

# === Define colors ===
length_colors <- c(
  "1830" = "#f6d6ff",
  "1383" = "#638ccc",
  "1803" = "#800233",
  "1245" = "#f9d42a"
)

# === Desired order ===
length_order <- c("1830", "1383", "1803", "1245")

# === Build dataset ===
legend_df <- expand.grid(
  LengthGroup = length_order,
  IsolateType = c("Historical", "Modern")
) %>%
  mutate(
    LengthGroup = factor(LengthGroup, levels = rev(length_order)), # top-to-bottom
    IsolateType = factor(IsolateType, levels = c("Historical", "Modern")), # left=H, right=M
    x = as.numeric(IsolateType),
    y = as.numeric(LengthGroup)
  )

# === Plot ===
p_legend <- ggplot(legend_df, aes(x = x, y = y, color = LengthGroup, alpha = IsolateType)) +
  geom_point(size = 6) +
  scale_color_manual(values = length_colors, name = "HTF haplotypes (bp)") +
  scale_alpha_manual(values = c("Historical" = 1, "Modern" = 0.4), name = "Sample type") +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("Historical", "Modern"),
    expand = expansion(mult = c(0.2, 0.2))
  ) +
  scale_y_continuous(
    breaks = 1:4,
    labels = rev(length_order),
    expand = expansion(mult = c(0.2, 0.2))
  ) +
  ggtitle("HTF haplotypes (bp)") +   # <-- main title at the top
  
  theme_minimal(base_size = 16) +
  theme(
    axis.text = element_text(size = 18, face = "bold"),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5), # centered bold title
    
    plot.margin = margin(30, 60, 30, 60)  # more padding
  )

# Save
ggsave(
  "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step4_additionalanalysis/step3_HTFfreq/legend_only_clean.pdf",
  p_legend, width = 4, height = 2.5
)