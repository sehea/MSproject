#Code used in the manuscript: Negative impacts of global change stressors permeate into deep soils 
#March 2025
#This file contains code for the analytical methods and visualization of the main figure results

#### Load packages ####
library(janitor)
library(tidyverse)
library(ggplot2)
library(MuMIn)
library(dplyr)
library(ggpmisc)
library(reshape2)
library(relaimpo)
library(MASS)
library(simpleboot)
library(corrplot)
library(lmerTest)#
library(readxl)
library(Matrix)#
library(lme4)
library(emmeans)
library(gvlma)
library(MASS)
library(caret)
library(car)
library(vegan)
library(randomForest)
library(ggalluvial)
library(viridis)
library(RColorBrewer)
library(psych)
library(ComplexHeatmap)# This is a bioconducter package and you need to install it using devtools
#library(devtools)
#install_github("jokergoo/ComplexHeatmap")
library(circlize)
library(tidyr)
library(purrr)

#### Load data ####
raw_data<-read.csv("data/DATA_WEN.csv")

# Scale Multistressor(MS) values, calculate MS intensity and number of MS (nMS) over three thresholds
process_raw_data <- function(raw_data) {
  raw_data %>%
    janitor::clean_names() %>%
    group_by(depth) %>%
    group_split(.keep = TRUE) %>%
    set_names(c("A_Topsoil", "B_Subsoil", "C_Deepsoil")) %>%
    map(~{
      df <- .x %>%
        filter(!is.na(weighted_emf))
      ms_columns <- df[, 7:14]
      
      # Scale MS
      ms_scaled <- as.data.frame(lapply(ms_columns, function(x) {
        100 * (x - min(x)) / (max(x) - min(x))
      }))
      
      # Calculate MS_intensity (mean of scaled MS)
      ms_scaled$ms_intensity <- rowMeans(ms_scaled, na.rm = TRUE)/100
      
      # Calculate nMS exceeding thresholds
      ms_scaled <- ms_scaled %>%
        mutate(
          nMS25 = rowSums(ms_scaled[, 1:8] >= 25, na.rm = TRUE),
          nMS50 = rowSums(ms_scaled[, 1:8] >= 50, na.rm = TRUE),
          nMS75 = rowSums(ms_scaled[, 1:8] >= 75, na.rm = TRUE)
        )
      
      # Update dataframe with metadata and MS variables
      df_combined <- bind_cols(df[, c(1:6, 15:21)], ms_scaled)
      
      list(
        metadata = list(
          depth = unique(df$depth),
          n_samples = nrow(df),
          ms_columns = names(ms_columns),
          df_columns = names(df_combined)),
        DF = df_combined,
        MS = ms_scaled)
    })}

# Processed dataset for further analyses
DATAFRAMES <- process_raw_data(raw_data)


################################# Figure 1 #########################################
# Combine data from all depth layers (df_f1:dataframe for fig.1)
df_f1 <- bind_rows(lapply(c("A_Topsoil", "B_Subsoil", "C_Deepsoil"), function(layer) DATAFRAMES[[layer]]$DF))

# Fig.1b and 1c boxplot
plot_boxplot <- function(data, x_var, y_var, colors,y_axis_label, y_limits = c(0.1, 0.8)) {
  ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], colour = .data[[x_var]])) +
    geom_boxplot(alpha = 1, linewidth = 1.5, width = 0.4) +
    stat_summary(fun = mean, geom = "point", shape = "diamond", size = 5, color = "#5E6054") +
    geom_jitter(alpha = 0.4, size = 5, width = 0.4) +
    scale_color_manual(values = colors) +
    theme_classic(base_line_size = 1) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 1),
      axis.title = element_text(size = 20),
      axis.text.x = element_text(color = "black", size = 18),
      axis.text.y = element_text(color = "black", size = 18),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20)
    ) +
    labs(x = "", y = "", color = "Depth") +
    coord_cartesian(ylim = y_limits)
}
save_plot_pdf <- function(plot_obj, filename, width = 7, height = 5) {
  dir.create("output", showWarnings = FALSE)
  ggsave(filename = file.path("output", filename), plot = plot_obj, device = "pdf", width = width, height = height)
}

colors <- c("#808080", "#FF6600", "#61A0D7")

p1b <- plot_boxplot(df_f1, "depth", "ms_intensity", colors, y_axis_label = "Intensity of mutiple stress")
save_plot_pdf(p1b, "Fig1b.pdf")

p1c <- plot_boxplot(df_f1, "depth", "weighted_emf", colors, y_axis_label = "Soil multifunctionality")
save_plot_pdf(p1c, "Fig1c.pdf")

# Statistic test
stat_lme <- function(data, response) {
  formula <- reformulate(termlabels = c("depth", "(1|site)", "(1|veg_type)"), response = response)
  model <- lmer(formula, data)
  print(summary(model))
  print(summary(emmeans(model, pairwise ~ depth)))
  invisible(model)
}

stat_lme(df_f1, "ms_intensity")
stat_lme(df_f1, "weighted_emf")

# Fig 1d linear regression
depth_lm_analysis <- function(data) {
  map_dfr(split(data, data$depth), ~{
    s <- summary(lm(weighted_emf ~ ms_intensity, .x))
    tibble(depth = .x$depth[1], R2_adj = signif(s$adj.r.squared, 2),
           p_value = pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE))
  })
}

lm_summary <- depth_lm_analysis(df_f1)
lm_summary <- lm_summary %>%
  mutate(significance = case_when(p_value < 0.001 ~ "***", p_value < 0.01 ~ "**", p_value < 0.05~ "*", TRUE ~ "ns"),
         label = paste0("R² = ", R2_adj, significance), x_pos = 0.5,
         y_pos = case_when(depth == "A_Topsoil" ~ 0.55, depth == "B_Subsoil" ~ 0.53, depth == "C_Deepsoil" ~ 0.51))

p1d <- ggplot(data = df_f1, aes(x = ms_intensity, y = weighted_emf, group = depth, color = as.factor(depth))) +
  theme_bw() +
  geom_smooth(method = "lm", se = FALSE) +
  ylab("Soil Multifunctionality") + xlab("Intensity of multiple stresses") +
  scale_color_manual(values = c("#808080", "#FF6600", "#61A0D7"), labels = c("A_Topsoil", "B_Subsoil", "C_Deepsoil")) +
  theme(panel.grid = element_blank(),legend.title = element_blank(),legend.position = "top",legend.text = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),axis.title = element_text(size = 16)) +
  geom_text(data = lm_summary,aes(x = x_pos, y = y_pos, label = label),inherit.aes = FALSE,size = 5,color = "black")

ggsave("output/Fig1d.pdf", p1d, width = 6, height = 6, device = "pdf")

# Fig 1e bar plot
# Calculate the frequency of the site with MS (nMS > 1 site (%))
df_f1e <- df_f1 %>%
  pivot_longer(cols = starts_with("nMS"), names_to = "threshold", values_to = "nMS_value") %>%
  mutate(threshold = gsub("nMS", "", threshold)) %>%
  group_by(depth, threshold) %>%
  summarise(frequency = mean(nMS_value > 1, na.rm = TRUE) * 100, .groups = "drop")

p1e <- ggplot(df_f1e, aes(x = threshold, y = frequency, fill = depth)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.7) +
  scale_fill_manual(values = c("A_Topsoil" = "#C6C6C6", "B_Subsoil" = "#FF944D", "C_Deepsoil" = "#77B6ED")) +
  labs(y = "Frequency of the sites (%)", x = "Threshold of multiple stressors (%)", fill = "Depth") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

ggsave("output/Fig1e.pdf", p1e, width = 7, height = 5)


################################# Figure 2 ###############################################
# Relative importance of mutiple stressors
super.MOD <- imap(DATAFRAMES, function(iter, depth) {
  df <- iter$DF
  ms <- iter$MS
  
  # Modeling data
  dfi <- bind_cols(df["weighted_emf"],ms %>% dplyr::select(-ms_intensity))
  
  # Linear model
  mdl <- lm(weighted_emf ~ ., data = dfi)
  assumptions <- gvlma(mdl)
  vif_res <- vif(mdl)
  
  # Relative importance
  RIMP <- calc.relimp(mdl)
  R2 <- RIMP@R2
  RIMP_df <- enframe(RIMP@lmg, name = "variable", value = "lmg.imp") %>%
    mutate(perc.imp = lmg.imp / R2, ori.DF = depth)
  
  # MuMIn model selection
  options(na.action = "na.fail")
  best_models <- dredge(mdl, rank = "BIC", extra = "R^2") %>%
    subset(delta < 4) %>%
    mutate(ori.DF = depth)
  
  # Return results
  iter$DF <- df
  iter$MS <- ms
  
  list(
    assumptions = list(assumptions = assumptions, vif = vif_res),
    relative.imp = RIMP_df,
    R2 = R2,
    MuMIN.mod = as.data.frame(best_models),
    updated = iter
  )
})

# Combine and export results
rimp_result <- bind_rows(map(super.MOD, "relative.imp")) #result for fig.2
mumin_result <- bind_rows(map(super.MOD, "MuMIN.mod")) #result for fig.S2


# Sankey diagram_fig.2
sankey_data <- data.frame(
  Source = rimp_result$variable,
  Target = rimp_result$ori.DF,
  Value = rimp_result$perc.imp
)

sankey_data$Source<-factor(sankey_data$Source, levels = c("nMS25", "nMS50", "nMS75", "aridity","warming","maxt","extreme_events","seasonality","human_influence","salinity","p_h"))
sankey_data$Target <- factor(sankey_data$Target, levels = unique(rimp_result$ori.DF))

custom_colors <- c("#FDE724","#B4DD2B","#6BCD58","#35B778","#1E9C88","#25828D","#30678D","#3D4989","#6A5ACD","#472777","#440154")

p2<-ggplot(sankey_data, aes(axis1 = Source, axis2 = Target, y = Value)) +
  geom_alluvium(aes(fill = Source), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey90", color = "grey50") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_fill_manual(values = custom_colors) +
  theme_void() +
  theme(legend.position = "none")

ggsave("output/Fig2.pdf", p2, width = 5, height = 6, device = "pdf")
# The importance of each stressor group is calculated by summing its corresponding stress factors in rimp_result


############################# Figure 3 ###############################
# Combine data from all soil layers
all_df <- bind_rows(lapply(c("A_Topsoil", "B_Subsoil", "C_Deepsoil"), function(layer) DATAFRAMES[[layer]]$DF))

# Convert threshold indicators to long format
df_f3 <- melt(
  all_df,
  id.vars = c("sample_id", "weighted_emf", "depth"),
  measure.vars = c("nMS25", "nMS50", "nMS75"),
  variable.name = "threshold",
  value.name = "nMS")

# Statistics
stat_labels <- df_f3 %>%
  group_by(depth, threshold) %>%
  summarise(adj_r2 = summary(lm(weighted_emf ~ nMS))$adj.r.squared,
            pval = summary(lm(weighted_emf ~ nMS))$coefficients[2, 4],
            .groups = "drop") %>%
  mutate(stars = case_when(pval < 0.001 ~ "***",pval < 0.01 ~ "**",pval < 0.05  ~ "*",TRUE ~ ""),
    label = paste0("R² = ", signif(adj_r2, 2), stars))

# Merge statistical labels back into main data
df_f3 <- left_join(df_f3, stat_labels, by = c("depth", "threshold"))

# Plot
p3<-ggplot(df_f3, aes(x = nMS, y = weighted_emf, color = threshold)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_grid(threshold ~ depth) +
  geom_text(aes(label = label),
            x = Inf, y = -Inf, hjust = 1.05, vjust = -1.1,
            inherit.aes = FALSE, size = 4) +
  scale_color_manual(values = c("#808080", "#FF6600", "#61A0D7")) +
  labs(x = "Number of stressors exceeding threshold",y = "Soil Multifunctionality",olor = NULL) +
  theme_bw() +
  theme(legend.position = "top",legend.text = element_text(size = 14), axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16),strip.text = element_text(size = 14),panel.grid = element_blank())
 
ggsave("output/Fig3.pdf", p3, width = 6, height = 6, device = "pdf")


############################# Figure 4 ###############################
df_f4 <- read.csv("data/EMF&MS_heatmap_v21.csv")
layers <- c("A_Topsoil", "B_Subsoil", "C_Deepsoil")
ms_vars <- 6:21
emf_vars <- c(3:5, 28:35)

# Compute Spearman correlation coefficients and adjusted p-values
results <- map(layers, ~{
  data <- df_f4 %>% filter(Depth == .x)
  corr.test(data[, ms_vars], data[, emf_vars], method = "spearman", adjust = "holm")
  })

# Combine correlation and p-value matrices
res_mat <- do.call(cbind, map(results, "r"))
pval_mat <- do.call(cbind, map(results, "p"))

# Convert p-values to significance stars
pval_star <- matrix("", nrow = nrow(pval_mat), ncol = ncol(pval_mat))
pval_star[pval_mat < 0.001] <- "***"
pval_star[pval_mat >= 0.001 & pval_mat < 0.01] <- "**"
pval_star[pval_mat >= 0.01 & pval_mat < 0.05] <- "*"

col_fun <- colorRamp2(c(-1, 0, 0.5), c("#00BFFF", "white", "red"))
col_split <- rep(layers, each = 11)

# Plot correlation heatmap
p4<-Heatmap(res_mat, col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
        column_split = col_split,
        column_title_gp = gpar(fontsize = 12, col = c("#666666", "#06959C", "#2D7BEE")),
        column_names_gp = gpar(fontsize = 10),
        cell_fun = function(j, i, x, y, w, h, fill) {
          grid.text(pval_star[i, j], x, y, gp = gpar(fontsize = 10))
        })
pdf("output/Fig4.pdf", width = 9, height = 6)
draw(p4)
dev.off()


####################################Figure 5#########################################
# reorganize data (df_f5:dataframe for fig.5)
df_f5 <- all_df %>%
  dplyr::select(sample_id, depth, weighted_emf, aridity, warming, nMS25, nMS50, nMS75)

# Convert to long format for easier iteration
df_f5_long <- df_f5 %>%
  pivot_longer(cols = starts_with("nMS"), names_to = "MS_level", values_to = "MS_value") %>%
  mutate(threshold = case_when(MS_level == "nMS25" ~ "MS>25%", MS_level == "nMS50" ~ "MS>50%", MS_level == "nMS75" ~ "MS>75%"))

# Classify into high nMS, low nMS, and general environments
df_f5_classified <- df_f5_long %>%
  group_by(MS_level) %>%
  mutate(type = case_when(
    MS_value > median(MS_value, na.rm = TRUE) ~ "High nMS environment",
    MS_value <= median(MS_value, na.rm = TRUE) ~ "Low nMS environment"
  )) %>%
  ungroup() %>%
  bind_rows(df_f5_long %>% 
            mutate(type = "General environment (all sites)"))

# Plot
plot_stress_effect <- function(df_f5, stress_var, stress_label) {
  ggplot(df_f5, aes(x = .data[[stress_var]], y = weighted_emf, color = type)) +
    geom_smooth(method = "lm", se = FALSE) +
    facet_grid(threshold ~ depth) +
    labs(x = stress_label, y = "Soil multifunctionality") +
    scale_color_manual(values = c("General environment (all sites)" = "#686868", "High nMS environment" = "#E12F47", "Low nMS environment" = "#2D7BEE")) +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.title = element_blank(), legend.position = "top")
}

p5a <- plot_stress_effect(df_f5_classified, "aridity", "Stress level of aridity (%)")
ggsave("output/Fig5a.pdf", p5a, width = 6, height = 6, device = "pdf")

p5b <- plot_stress_effect(df_f5_classified, "warming", "Stress level of warming (%)")
ggsave("output/Fig5b.pdf", p5b, width = 6, height = 6, device = "pdf")


####################################Table 1#########################################
# Function to fit a linear model and extract coefficient for the single stressor (i.e., effect size of aridity/warming) under high nMS, low nMS, and general environments
get_lm_result <- function(data, stressor, label) {
  data <- data %>%
    mutate(across(all_of(stressor), ~ . / 100))
  
  formula <- reformulate(stressor, response = "weighted_emf")
  model <- lm(formula = formula, data = data)
  
  result <- as.data.frame(summary(model)$coefficients)
  result <- result[rownames(result) == stressor, , drop = FALSE]
  
  result$lm <- label
  return(result)
}

# Analyze single climatic stressor (i.e., aridity and warming) effect across general, high nMS, and low nMS environment for a specific depth
analyze_stressor <- function(df, depth, stressor, ms_type_label) {
  df_depth <- df[df$depth == depth, ]
  results <- list()
  
  # General environment
  label_all <- paste(stressor, "under General environment")
  results[["ALL"]] <- get_lm_result(df_depth[df_depth$type == "General environment (all sites)", ], stressor, label_all)
  
  # High nMS environment
  label_high <- paste(stressor, "under High", ms_type_label)
  results[["HIGH"]] <- get_lm_result(df_depth[df_depth$type == "High nMS environment", ], stressor, label_high)
  
  # Low nMS environment
  label_low <- paste(stressor, "under Low", ms_type_label)
  results[["LOW"]] <- get_lm_result(df_depth[df_depth$type == "Low nMS environment", ], stressor, label_low)
  
  combined <- do.call(rbind, results)
  combined$depth <- depth
  combined$stressor <- stressor
  combined$nms_level <- ms_type_label
  return(combined)
}

# Main function to iterate over all combinations of stressors × depths × MS thresholds
run_all_models <- function(df, stressors = c("aridity", "warming"), depths = unique(df$depth), ms_levels = c("MS>25%", "MS>50%", "MS>75%")) {
  all_results <- list()
  
  for (stressor in stressors) {
    for (ms_label in ms_levels) {
      df_sub <- df[df$threshold == ms_label, ]
      for (d in depths) {
        res <- analyze_stressor(df_sub, d, stressor, ms_label)
        all_results[[paste(stressor, ms_label, d, sep = "_")]] <- res
      }
    }
  }
  
  return(do.call(rbind, all_results))
}

final_results_f5 <- run_all_models(df_f5_classified)

write.csv(final_results_f5, "output/Table1_effect_size_result.csv", row.names = T)
