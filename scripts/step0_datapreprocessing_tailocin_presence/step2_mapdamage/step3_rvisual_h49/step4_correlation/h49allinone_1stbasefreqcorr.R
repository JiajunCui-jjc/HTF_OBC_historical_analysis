#!/usr/bin/env Rscript
# -------------------------------------------------------------------------
# Compare first-base 5′ C to T frequencies between Pseudomonas and Arabidopsis
# across 49 historical samples, and test correlation + time regression.
# Updated for 2025_summerpaper_all file structure.
# -------------------------------------------------------------------------

# ====== Libraries ======
library(ggplot2)
library(ggpmisc)
library(reshape2)
library(scales)
library(dplyr)
library(stringr)

if (!requireNamespace("ggrepel", quietly = TRUE)) {
  message("Tip: install.packages('ggrepel') for better non-overlapping labels.")
}

# ====== Directories ======
setwd("/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step2_mapdamageplot/suppfig3mapdamage")

# Input C to T frequency files
filePs <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step2_mapdamageplot/tops/allinoneCtoT.txt"
fileAt <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/step0_datapreprocessing_tailocin_presence/step2_mapdamageplot/toAt/allinoneCtoT.txt"

# Historical sample list (short names)
file_samples <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/data/historical49.txt"

# Year metadata file
file_dates <- "/Users/cuijiajun/Desktop/others/tmphernan/2025_summerpaper_all/2025_summer_paperfig_m57/results/figureandtable/supptables/Supptable_h49.txt" # temp PL0087 is NA

# ====== Column names for consistency ======
colnames_v <- c("pos", "ref_base", "value", "group", "total",
                "proportion", "pos_from_end", "sample")

# ====== Load data ======
tablePs <- read.table(filePs, header = TRUE, sep = "\t")
tableAt <- read.table(fileAt, header = FALSE, sep = "\t", col.names = colnames_v)

# ====== Clean sample names to match historical49.txt ======
trim_sample_name <- function(x) {
  x <- sub("\\.mapped_to.*", "", x)
  x <- sub("\\.markeddup.*", "", x)
  x <- sub("_mapped.*", "", x)
  x <- sub("_markdup.*", "", x)
  x
}
tablePs$sample <- trim_sample_name(tablePs$sample)
tableAt$sample <- trim_sample_name(tableAt$sample)

# ====== Extract first base (position 1) ======
tablePs1st <- tablePs %>% filter(pos == 1) %>% arrange(sample)
tableAt1st <- tableAt %>% filter(pos == 1) %>% arrange(sample)

tmp1 <- data.frame(
  sample = tablePs1st$sample,
  Ps_1stbaseCtoTfreq = as.numeric(tablePs1st$proportion),
  At_1stbaseCtoTfreq = as.numeric(tableAt1st$proportion)
)

# ====== Add sample years ======
dates <- read.table(file_dates, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
dates$YEAR <- suppressWarnings(as.numeric(dates$YEAR))

samples49 <- read.table(file_samples, header = FALSE, stringsAsFactors = FALSE)[[1]]

dates49 <- dates %>%
  filter(SAMPLE %in% tmp1$sample & SAMPLE %in% samples49) %>%
  select(SAMPLE, YEAR) %>%
  arrange(SAMPLE)

filtered_tmp1 <- tmp1 %>%
  filter(sample %in% dates49$SAMPLE) %>%
  arrange(sample)

all_combined <- cbind(filtered_tmp1, dates49)
colnames(all_combined) <- c("sample", "Ps_1stbaseCtoTfreq", "At_1stbaseCtoTfreq", "samplename", "dates")

write.table(all_combined, "49samples_firstbase_CtoT_with_dates.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("✅ Combined data (after removing invalid YEARs):", nrow(all_combined), "samples\n")

# ====== Safe helper to fit linear model ======
safe_lm <- function(formula, data) {
  if (nrow(data) < 3) return(NULL)
  fit <- tryCatch(lm(formula, data = data), error = function(e) NULL)
  return(fit)
}

# ===========================================================
# 1️⃣ Correlation: Pseudomonas vs Arabidopsis
# ===========================================================
lm_pa <- safe_lm(At_1stbaseCtoTfreq ~ Ps_1stbaseCtoTfreq, all_combined)
if (!is.null(lm_pa)) {
  summ_pa <- summary(lm_pa)
  R2_pa <- summ_pa$r.squared
  p_pa <- coef(summ_pa)[2, "Pr(>|t|)"]
} else {
  R2_pa <- NA; p_pa <- NA
}

p1 <- ggplot(all_combined, aes(Ps_1stbaseCtoTfreq, At_1stbaseCtoTfreq)) +
  geom_point(size = 4, alpha = 1) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
  scale_x_continuous(labels = label_number(accuracy = 0.01), limits = c(0, 0.04)) +
  scale_y_continuous(labels = label_number(accuracy = 0.01), limits = c(0, 0.04)) +
  labs(
    title = expression("First-base C to T frequency: " * italic("Pseudomonas") * " sp. vs " * italic("A. thaliana")),
    subtitle = sprintf("Linear fit: R² = %.2f, p = %.3g", R2_pa, p_pa),
    x = expression(italic("Pseudomonas") * " sp. (C to T at 5 end)"),
    y = expression(italic("Arabidopsis thaliana") * " (C to T at 5 end)")
  ) +
  theme_bw(base_size = 18, base_family = "Helvetica") +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 18),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 20, face = "bold"),
    panel.grid = element_blank()
  )

ggsave("49_CtoTcor_Pseudomonas_vs_Arabidopsis.pdf", p1, width = 9.6, height = 7.6)

# Remove NA or invalid years
all_combined <- all_combined %>%
  filter(!is.na(dates) & is.finite(dates))

# ===========================================================
# 2️⃣ Arabidopsis vs year
# ===========================================================
lm_at <- lm(At_1stbaseCtoTfreq ~ dates, data = all_combined)
summ_at <- summary(lm_at)
R2_at <- summ_at$r.squared
p_at_lm <- coef(summ_at)[2, "Pr(>|t|)"]

p_at <- ggplot(all_combined, aes(dates, At_1stbaseCtoTfreq)) +
  geom_point(color = "#D04F4F", size = 4, alpha = 1) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    title = expression(italic("Arabidopsis thaliana") * " : first-base C to T vs year"),
    subtitle = sprintf("Linear fit: R² = %.2f, p = %.3g", R2_at, p_at_lm),
    x = "Year", y = "C to T frequency"
  ) +
  theme_bw(base_size = 18, base_family = "Helvetica") +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 20),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 20, face = "bold"),
    panel.grid = element_blank()
  )

ggsave("CtoT_vs_year_Arabidopsis.pdf", p_at, width = 8.6, height = 7.6)

# ===========================================================
# 3️⃣ Pseudomonas vs year
# ===========================================================
lm_ps <- lm(Ps_1stbaseCtoTfreq ~ dates, data = all_combined)
summ_ps <- summary(lm_ps)
R2_ps <- summ_ps$r.squared
p_ps_lm <- coef(summ_ps)[2, "Pr(>|t|)"]

p_ps <- ggplot(all_combined, aes(dates, Ps_1stbaseCtoTfreq)) +
  geom_point(color = "#D04F4F", size = 4, alpha = 1) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    title = expression(italic("Pseudomonas") * " sp.: first-base C to T vs year"),
    subtitle = sprintf("Linear fit: R² = %.2f, p = %.3g", R2_ps, p_ps_lm),
    x = "Year", y = "C to T frequency"
  ) +
  theme_bw(base_size = 18, base_family = "Helvetica") +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 20),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 20, face = "bold"),
    panel.grid = element_blank()
  )

ggsave("CtoT_vs_year_Pseudomonas.pdf", p_ps, width = 8.6, height = 7.6)

# ===========================================================
# 4️⃣ Save regression summary
# ===========================================================
reg_summary <- data.frame(
  model = c("Arabidopsis thaliana ~ year", "Pseudomonas sp. ~ year"),
  slope = c(coef(lm_at)[["dates"]], coef(lm_ps)[["dates"]]),
  intercept = c(coef(lm_at)[["(Intercept)"]], coef(lm_ps)[["(Intercept)"]]),
  r_squared = c(R2_at, R2_ps),
  p_value = c(p_at_lm, p_ps_lm),
  n = c(nobs(lm_at), nobs(lm_ps))
)

write.table(reg_summary, "CtoT_regression_summaries.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n✅ Done — Saved all plots and tables in:", getwd(), "\n")
