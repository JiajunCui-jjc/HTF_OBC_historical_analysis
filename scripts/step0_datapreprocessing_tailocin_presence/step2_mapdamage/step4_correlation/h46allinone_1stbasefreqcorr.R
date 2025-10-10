# ====== Libraries ======
library(ggplot2)
library(ggpmisc)    # stat_poly_eq
library(reshape2)   # melt
library(scales)     # percent_format
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  message("Tip: install.packages('ggrepel') for better non-overlapping labels.")
}
setwd('/Users/cuijiajun/Desktop/others/tmphernan/2025_Aug_balgiumtalk_materials_jiajun/material1_mapdamage_h34/v1_h46/step2_correlation/')
tablePs<-read.table('../step1_h46ps/h46_CtoTfreq_Ps.txt',header = FALSE)
tableAt<-read.table('../step1_h46ps/h46_CtoTfreq_At.txt',header = FALSE)
#you need make sure samplename is the same in two tables
colnames(tablePs)<-c("pos"         ,  "value"        , "group"    ,     "Pos"    ,      "total"     ,    "base"        ,  "proportion",    "pos_from_5end"
                     ,"sample"  )
colnames(tableAt)<-c("pos"         ,  "value"        , "group"    ,     "Pos"    ,      "total"     ,    "base"        ,  "proportion",    "pos_from_5end"
                     ,"sample"   )
# For Pseudomonas table
tablePs$sample <- dplyr::case_when(
  # Historical samples: extract anything between "results_" and "_mapped"
  grepl("^results_[0-9]", tablePs$sample) ~ 
    sub("^results_([^_]+.*?)(?:_removalAt.*)?$", "\\1", tablePs$sample),
  
  # HB samples with prefix "results_HB"
  grepl("^results_HB[0-9]+", tablePs$sample) ~ 
    sub("^results_(HB[0-9]+).*", "\\1", tablePs$sample),
  
  
  # Modern PL samples
  grepl("^PL[0-9]{4}", tablePs$sample) ~ 
    sub("^(PL[0-9]{4}).*", "\\1", tablePs$sample),
  
  TRUE ~ tablePs$sample
)
# For Arabidopsis table
tableAt$sample <- dplyr::case_when(
  # Historical samples: extract anything between "results_" and "_mapped"
  grepl("^results_[0-9]", tableAt$sample) ~ 
    sub("^results_([^_]+.*?)(?:_mapped.*)?$", "\\1", tableAt$sample),
  
  # HB samples with prefix "results_HB"
  grepl("^results_HB[0-9]+", tableAt$sample) ~ 
    sub("^results_(HB[0-9]+).*", "\\1", tableAt$sample),
  
  
  # Modern PL samples
  grepl("^PL[0-9]{4}", tableAt$sample) ~ 
    sub("^(PL[0-9]{4}).*", "\\1", tableAt$sample),
  
  TRUE ~ tableAt$sample
)

tablePs1st<-tablePs[tablePs$pos==1,]
tableAt1st<-tableAt[tableAt$pos==1,]
tablePs1st <- tablePs1st[order(tablePs1st$sample), ]
tableAt1st <- tableAt1st[order(tableAt1st$sample), ]
tmp1<-data.frame(sample=tablePs1st$sample,Ps_1stbaseCtoTfreq=tablePs1st$proportion,At_1stbaseCtoTfreq=tableAt1st$proportion)
tmp1$Ps_1stbaseCtoTfreq <- as.numeric(tmp1$Ps_1stbaseCtoTfreq)
tmp1$At_1stbaseCtoTfreq <- as.numeric(tmp1$At_1stbaseCtoTfreq)

#then add dates:
dates<-read.table('/Users/cuijiajun/Desktop/2023-2024\ PhD\ ucl/2024_aMeta/wholepipeAt_Ps/2024_233_analysis/stats/tmp489_withdatesandlocs_uniq.txt',sep='\t',header = TRUE)
tmp1$sample
dates46<-data.frame(sample=dates[dates$samplename%in% tmp1$sample,]$samplename,dates=dates[dates$samplename%in% tmp1$sample,]$year)
#46 samples

tmp1[tmp1$sample%in% dates46$sample,]
# Sort dates46 by sample to ensure the order matches
dates46 <- dates46[order(dates46$sample),]

# Filter and sort tmp1 to match the samples and order in dates46
filtered_tmp1 <- tmp1[tmp1$sample %in% dates46$sample,]
filtered_tmp1 <- filtered_tmp1[order(filtered_tmp1$sample),]

# Combine the two data frames with the order of sample column
all_combined <- cbind(filtered_tmp1, dates46)
all<-all_combined[,c(1,2,3,5)]
write.table(all, "ctot_firstbaseseq_withdates.txt", sep = "\t", row.names = FALSE, quote = FALSE)
# ====== Colors & theme ======
col_pseudo <- "#D55E00"  # blue
  col_arabi  <- "#D55E00"  # vermillion
    col_fit    <- "grey20"
      
    theme_conf <- function(base_size = 12){
      theme_bw(base_size = base_size) +
        theme(
          panel.grid.major = element_line(size = 0.2, colour = "grey90"),
          panel.grid.minor = element_blank(),
          panel.border     = element_blank(),
          axis.line        = element_line(colour = "black"),
          plot.title       = element_text(face = "bold", hjust = 0.5, margin = margin(b = 6), size = 14),
          plot.subtitle    = element_text(hjust = 0.5, colour = "grey30", margin = margin(b = 6), size = 14),
          legend.position  = "top",
          legend.title     = element_blank()
        )
    }
    
    # ===========================================================
    # 1) Correlation: Pseudomonas vs Arabidopsis (first base)
    # ===========================================================
    # --- Linear regression for the Ps-vs-Arabidopsis scatter ---
    # Model: Arabidopsis (Y) ~ Pseudomonas (X)
    lm_pa   <- lm(At_1stbaseCtoTfreq ~ Ps_1stbaseCtoTfreq, data = tmp1)
    summ_pa <- summary(lm_pa)
    R2_pa   <- summ_pa$r.squared
    p_pa    <- coef(summ_pa)[2, "Pr(>|t|)"]   # p-value for slope
    
    # --- Plot: same style as your year regressions, with R² & p in subtitle ---
    p1 <- ggplot(tmp1, aes(Ps_1stbaseCtoTfreq, At_1stbaseCtoTfreq)) +
      geom_point(size = 3, alpha = 1) +
      geom_smooth(method = "lm", se = F, color = col_fit, linewidth = 0.6, alpha = 0.12) +
      scale_x_continuous(labels = scales::label_number(accuracy = 0.01), limits = c(0, 0.04)) +
      scale_y_continuous(labels = scales::label_number(accuracy = 0.01), limits = c(0, 0.04)) +
      labs(
        title    = expression("First-base C to T frequency: " * italic("Pseudomonas") * " sp. vs " * italic("Arabidopsis thaliana")),
        subtitle = sprintf("Linear fit: R² = %.2f, p = %.3g", R2_pa, p_pa),
        x        = expression(italic("Pseudomonas") * " sp. (C to T at 5' base)"),
        y        = expression(italic("Arabidopsis thaliana")  * "  (C to T at 5' base)")
      ) +
      theme_conf()+
      theme(
      plot.title    = element_text(size = 18, face = "bold", hjust = 1.5),
      plot.subtitle = element_text(size = 16, hjust = 0.5),
      axis.title    = element_text(size = 16, face = "bold"),
      axis.text     = element_text(size = 14)
    )
    
    ggsave("46_CtoTcor_Pseudomonas_vs_Arabidopsis.pdf", p1, width = 8.6, height = 7.6)
    # ===========================================================
    # 3) Arabidopsis vs year — regression & stats
    # ===========================================================
    lm_at   <- lm(At_1stbaseCtoTfreq ~ dates, data = all)
    summ_at <- summary(lm_at)
    R2_at   <- summ_at$r.squared
    p_at_lm <- coef(summ_at)[2, "Pr(>|t|)"]
    
    p_at <- ggplot(all, aes(dates, At_1stbaseCtoTfreq)) +
      geom_point(color = col_arabi, size = 3, alpha = 1) +
      geom_smooth(method = "lm", se = F, color = col_fit, linewidth = 0.7, alpha = 0.12) +
      scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
      labs(
        title    = expression(italic("Arabidopsis thaliana") * " : first-base C to T vs year"),
        subtitle = sprintf("Linear fit: R² = %.2f, p = %.3g", R2_at, p_at_lm),
        x = "Year", y = "C to T frequency"
      ) +
      theme_conf()+
      theme(
      plot.title    = element_text(size = 20, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 16, hjust = 0.5),
      axis.title    = element_text(size = 16, face = "bold"),
      axis.text     = element_text(size = 14)
    )
    ggsave("CtoT_vs_year_Arabidopsis.pdf", p_at, width = 8.6, height = 7.6)

    # ===========================================================
    # 4) Pseudomonas vs year — regression & stats
    # ===========================================================
    lm_ps   <- lm(Ps_1stbaseCtoTfreq ~ dates, data = all)
    summ_ps <- summary(lm_ps)
    R2_ps   <- summ_ps$r.squared
    p_ps_lm <- coef(summ_ps)[2, "Pr(>|t|)"]
    
    p_ps <- ggplot(all, aes(dates, Ps_1stbaseCtoTfreq)) +
      geom_point(color = col_pseudo, size = 3, alpha = 1) +
      geom_smooth(method = "lm", se = F, color = col_fit, linewidth = 0.7, alpha = 0.12) +
      scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
      labs(
        title    = expression(italic("Pseudomonas") * " sp.: first-base C to T vs year"),
        subtitle = sprintf("Linear fit: R² = %.2f, p = %.3g", R2_ps, p_ps_lm),
        x = "Year", y = "C to T frequency"
      ) +
      theme_conf()+
     theme(
      plot.title    = element_text(size = 20, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 16, hjust = 0.5),
      axis.title    = element_text(size = 16, face = "bold"),
      axis.text     = element_text(size = 14)
    )
    ggsave("CtoT_vs_year_Pseudomonas.pdf", p_ps, width = 8.6, height = 7.6)

    # ===========================================================
    # 5) Save regression stats table
    # ===========================================================
    reg_summary <- data.frame(
      model     = c("Arabidopsis thaliana ~ year", "Pseudomonas sp. ~ year"),
      slope     = c(coef(lm_at)[["dates"]], coef(lm_ps)[["dates"]]),
      intercept = c(coef(lm_at)[["(Intercept)"]], coef(lm_ps)[["(Intercept)"]]),
      r_squared = c(R2_at, R2_ps),
      p_value   = c(p_at_lm, p_ps_lm),
      n         = c(nobs(lm_at), nobs(lm_ps))
    )
    write.table(reg_summary, "CtoT_regression_summaries.tsv", sep = "\t", row.names = FALSE, quote = FALSE)