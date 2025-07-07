
library(dplyr)
library(ggplot2)
library(viridis)


# READ IN 
Master_df_GY<- read.csv("Definitive_GY_36_GBM.fit_MASTER_METRICS_DF_27.03.2025.csv")
Master_df_TGW <- read.csv("Definitive_TGW_36_GBM.fit_MASTER_METRICS_DF_27.03.2025.csv")
Master_df_GPC<- read.csv("Definitive_GPC_36_GBM.fit_MASTER_METRICS_DF_27.03.2025.csv")
Master_df_HET <-  read.csv("Definitive_HET_36_GBM.fit_MASTER_METRICS_DF_27.03.2025.csv")

####  MASTER DATAFRAME FOR PLOTTING ####


Master_df <-  rbind(Master_df_GY,Master_df_GPC, Master_df_TGW, Master_df_HET)


#### ICC ####
# Filter the data for ICC metric
icc_data  <- Master_df %>% 
  filter(Metric =="ICC_Prediction_Stability")

# Summarize the data by computing the mean, min, and max for each combination of BoostRounds, LearningRate, and Phenotype
# Summarize the data: calculating mean, standard error, and min/max
summary_df <- icc_data %>%
  group_by(BoostRounds, LearningRate, Phenotype) %>%
  summarize(
    mean_ICC = mean(Value), # Mean value of ICC - 5 for each combination of BoostRounds, LearningRate
    sd_ICC = sd(Value),
    se_ICC = sd(Value) / sqrt(n()),
    min_ICC = min(Value),
    max_ICC = max(Value),
    .groups = "drop"
  )




# Plot the data with error bars and a trend line
ggplot(summary_df, aes(x = BoostRounds, y = mean_ICC, color = as.factor(LearningRate))) +
  geom_line(linewidth = 0.05, alpha= 0.5, linetype = "dashed") +  # Line connecting the mean values
  geom_errorbar(aes(
    ymin = mean_ICC -  se_ICC,
    ymax = mean_ICC +  se_ICC
  ), width = 100) +  
  geom_smooth(data = icc_data, aes(x = as.numeric(BoostRounds), y = Value, color = as.factor(LearningRate), fill =as.factor(LearningRate)) ,
              method = "loess", level = 0.95, span = 1.0,    linetype = "dotted", alpha = 0.1, linewidth = 0.5, show.legend = FALSE ) +
  geom_point(data = icc_data, aes(x = BoostRounds, y = Value, color = as.factor(LearningRate)), 
             position = position_jitter(width = 0.1), alpha = 0.8, size =0.75) + 
  facet_wrap(summary_df$Phenotype <- factor(summary_df$Phenotype, levels = c("GY_1", "TGW_1", "GPC_1", "HET_1"))) + 
  labs(title = "ICC vs BoostRounds by LearningRate",
       x = "BoostRounds",
       y = "ICC",
       color = "Learning Rate") +
  theme_minimal()



# Filter the data for ICC metric
icc_data<-  Master_df %>% 
  filter(Metric== "ICC_Prediction_Stability")

# Summarize the data by computing the mean ICC for each combination of BoostRounds, LearningRate, and Phenotype
summary_df <- icc_data %>%
  group_by(BoostRounds, LearningRate, Phenotype) %>%
  summarize(
    mean_ICC = mean(Value),
    .groups = "drop"
  )

# Create a heatmap: x-axis is BoostRounds, y-axis is LearningRate, fill is mean ICC.
ggplot(summary_df, aes(x = BoostRounds, y = as.factor(LearningRate), fill = mean_ICC)) +
  geom_tile() +
  facet_wrap(summary_df$Phenotype <- factor(summary_df$Phenotype, levels = c("GY_1", "TGW_1", "GPC_1"))) + # "HET_1"
  scale_fill_viridis_c(option = "plasma", name = "Mean ICC") +
  labs(title = "Heatmap: ICC vs BoostRounds by LearningRate",
       x = "BoostRounds",
       y = "Learning Rate") +
  theme_minimal()


#### FLEISS KAAPPA ####

# Filter and summarize for Fleiss Kappa metric
FLEISS_data<- Master_df %>% 
  filter(Metric =="Fleiss_Kappa_Prediction_Stability")

summary_df <- FLEISS_data %>%
  group_by(BoostRounds, LearningRate, Phenotype) %>%
  summarize(
    mean_FLEISS = mean(Value),
    sd_FLEISS = sd(Value),
    se_FLEISS = sd(Value) / sqrt(n()),
    min_FLEISS = min(Value),
    max_FLEISS = max(Value),
    .groups = "drop"
  )

# Adjusted Plot 1: Fleiss Kappa vs BoostRounds by LearningRate
ggplot(summary_df, aes(x = BoostRounds, y = mean_FLEISS, color = as.factor(LearningRate))) +
  geom_line(linewidth = 0.05, alpha = 0.5, linetype = "dashed") +  
  geom_errorbar(aes(
    ymin = mean_FLEISS - se_FLEISS,
    ymax = mean_FLEISS + se_FLEISS
  ), width = 100) +
  geom_smooth(data = FLEISS_data, 
              aes(x = as.numeric(BoostRounds), y = Value, color = as.factor(LearningRate), 
                  fill = as.factor(LearningRate)),
              method = "loess", level = 0.95, span = 1.0,
              linetype = "dotted", alpha = 0.1, linewidth = 0.5, show.legend = FALSE) +
  geom_point(data = FLEISS_data, 
             aes(x = BoostRounds, y = Value, color = as.factor(LearningRate)),
             position = position_jitter(width = 0.1), alpha = 0.8, size = 0.75) +
  facet_wrap(summary_df$Phenotype <- factor(summary_df$Phenotype, levels = c("GY_1", "TGW_1", "GPC_1", "HET_1"))) + # 
  labs(title = "FLEISS vs BoostRounds by LearningRate",
       x = "BoostRounds",
       y = "FLEISS",
       color = "Learning Rate") +
  theme_minimal()





# Filter the data for FLEISS metric
FLEISS_data <-  Master_df %>% 
  filter(Metric == "Fleiss_Kappa_Prediction_Stability")

# Summarize the data by computing the mean FLEISS for each combination of BoostRounds, LearningRate, and Phenotype
summary_df <- FLEISS_data %>%
  group_by(BoostRounds, LearningRate, Phenotype) %>%
  summarize(
    mean_FLEISS = mean(Value),
    .groups = "drop"
  )

# Create a heatmap: x-axis is BoostRounds, y-axis is LearningRate, fill is mean FLEISS.
ggplot(summary_df, aes(x = BoostRounds, y = as.factor(LearningRate), fill = mean_FLEISS)) +
  geom_tile() +
  facet_wrap(~Phenotype) +
  scale_fill_viridis_c(option = "plasma", name = "Mean FLEISS") +
  labs(title = "Heatmap: FLEISS vs BoostRounds by LearningRate",
       x = "BoostRounds",
       y = "Learning Rate") +
  theme_minimal()




#### Pearson ####
# Plot 3 
# Filter the data for ICC metric
pearson_data <- Master_df %>% 
  filter(Metric == "Pearson")

# Summarize the data by computing the mean, min, and max for each combination of BoostRounds, LearningRate, and Phenotype
# Summarize the data: calculating mean, standard error, and min/max
summary_df <- pearson_data %>%
  group_by(BoostRounds, LearningRate, Phenotype) %>%
  summarize(
    mean_pearson = mean(Value),
    sd_pearson = sd(Value),
    se_pearson = sd(Value) / sqrt(n()),
    min_pearson = min(Value),
    max_pearson = max(Value),
    .groups = "drop"
  ) 

summary_df_subset <- pearson_data %>%
  group_by(BoostRounds, LearningRate, Phenotype, Subset) %>%
  summarize(
    mean_pearson_subset = mean(Value),
    sd_pearson_subset = sd(Value),
    se_pearson_subset = sd(Value) / sqrt(n()),
    min_pearson_subset = min(Value),
    max_pearson_subset = max(Value),
    .groups = "drop"
  )


summary_df$Phenotype <- factor(summary_df$Phenotype, levels = c("GY_1", "TGW_1", "GPC_1", "HET_1"))

#### PRESENT ME
ggplot(summary_df, aes(x = BoostRounds, y = mean_pearson, color = as.factor(LearningRate))) +
  geom_line(position = position_dodge(width = 50)) +  # Line connecting the mean values
  geom_errorbar(aes(ymin = mean_pearson - se_pearson, ymax = mean_pearson + se_pearson), width = 500, position = position_dodge(width = 50)) +  # Error bars SE
  geom_smooth(method = "loess", se = FALSE, linetype = "dotted", alpha = 0.1, linewidth = 0.5) +  # Overall trend line (line of best fit)
    facet_wrap(~Phenotype, scale = "free_y") + 
  labs(title = "Pearson Correlation vs BoostRounds by LearningRate",
       x = "BoostRounds",
       y = "Pearson Correlation",
       color = "Learning Rate",
       subtitle = "Mean Pearson Correlation per Subset (5)") +
  theme_minimal()









#### R Squared ####

# Filter the data for R_squared metric
r2_data<- Master_df %>% 
  filter(Metric == "R_squared") 

# Summarize the data by computing the mean, standard error, and min/max for each combination of BoostRounds, LearningRate, and Phenotype
summary_df <- r2_data %>%
  group_by(BoostRounds, LearningRate, Phenotype) %>%
  summarize(
    mean_r2 = mean(Value),
    sd_r2 = sd(Value),
    se_r2 = sd(Value) / sqrt(n()),
    min_r2 = min(Value),
    max_r2 = max(Value),
    .groups = "drop"
  )

summary_df$Phenotype <- factor(summary_df$Phenotype, levels = c("GY_1", "TGW_1", "GPC_1", "HET_1"))

# Plot the data with error bars and a trend line
ggplot(summary_df, aes(x = BoostRounds, y = mean_r2, color = as.factor(LearningRate))) +
  geom_line() +  # Line connecting the mean R_squared values
  geom_errorbar(aes(ymin = mean_r2 - se_r2, ymax = mean_r2 + se_r2), width = 100) +  # Error bars SE
  geom_smooth(method = "loess", se = FALSE, linetype = "dotted", alpha = 0.1, linewidth = 0.5) +  # Overall trend line (line of best fit)
  geom_point(data = summary_df, aes(x = BoostRounds, y = mean_r2, color = as.factor(LearningRate))) + 
  facet_wrap(~Phenotype, scale = "free_y") +
  labs(title = "R-squared vs BoostRounds by LearningRate",
       x = "BoostRounds",
       y = "R-squared",
       color = "Learning Rate") +
  theme_minimal()



#### AUC ####

# Filter the data for AUC metric
auc_data <- Master_df %>% 
  filter(Metric == "AUC")

# Summarize the data by computing the mean, standard error, and min/max for each combination of BoostRounds, LearningRate, and Phenotype
summary_df <- auc_data %>%
  group_by(BoostRounds, LearningRate, Phenotype) %>%
  summarize(
    mean_auc = mean(Value),
    sd_auc = sd(Value),
    se_auc = sd(Value) / sqrt(n()),
    min_auc = min(Value),
    max_auc = max(Value),
    .groups = "drop"
  )

summary_df$Phenotype <- factor(summary_df$Phenotype, levels = c("GY_1", "TGW_1", "GPC_1", "HET_1"))


# Plot the data with error bars and a trend line (Plot 1)
ggplot(summary_df, aes(x = BoostRounds, y = mean_auc, color = as.factor(LearningRate))) +
  geom_line() +  # Line connecting the mean AUC values
  geom_errorbar(aes(ymin = mean_auc - se_auc, ymax = mean_auc + se_auc), width = 100) +  # Error bars SE
  geom_smooth(method = "loess", se = FALSE, linetype = "dotted", alpha = 0.1, linewidth = 0.5) +  # Overall trend line
  geom_point(data = summary_df, aes(x = BoostRounds, y = mean_auc, color = as.factor(LearningRate))) +
  facet_wrap(~Phenotype, scale = "free_y" ) + 
  labs(title = "AUC vs BoostRounds by LearningRate",
       x = "BoostRounds",
       y = "AUC",
       color = "Learning Rate") +
  theme_minimal()







#### NDCG 20% ####


# Filter the data for NDCG_at_20_percent metric
ndcg_data <-  Master_df %>% 
  filter(Metric== "NDCG_at_20_percent")

# Summarize the data by computing the mean, standard error, and min/max for each combination of BoostRounds, LearningRate, and Phenotype
summary_df <- ndcg_data %>%
  group_by(BoostRounds, LearningRate, Phenotype) %>%
  summarize(
    mean_ndcg = mean(Value),
    sd_ndcg = sd(Value),
    se_ndcg = sd(Value) / sqrt(n()),
    min_ndcg = min(Value),
    max_ndcg = max(Value),
    .groups = "drop"
  )

summary_df$Phenotype <- factor(summary_df$Phenotype, levels = c("GY_1", "TGW_1", "GPC_1", "HET_1"))

# Plot the data with error bars and a trend line (Plot 1)
ggplot(summary_df, aes(x = BoostRounds, y = mean_ndcg, color = as.factor(LearningRate))) +
  geom_line() +  # Line connecting the mean NDCG values
  geom_errorbar(aes(ymin = mean_ndcg - se_ndcg, ymax = mean_ndcg + se_ndcg ), width = 100) +  # Error bars for the range 
  geom_smooth(method = "loess", se = FALSE, linetype = "dotted", alpha = 0.1, linewidth = 0.5) +  # Overall trend line
  geom_point(data = summary_df, aes(x = BoostRounds, y = mean_ndcg, color = as.factor(LearningRate))) +
  facet_wrap(~Phenotype, scale = "free_y") + 
  labs(title = "NDCG@20% vs BoostRounds by LearningRate",
       x = "BoostRounds",
       y = "NDCG@20%",
       color = "Learning Rate") +
  theme_minimal()



 



#### ALL TRAITS ####
# ICC 
# Filter the data for ICC metric
icc_data <- Master_df %>% 
  filter(Metric == "ICC_Prediction_Stability")

# Summarize the data by computing the mean, min, and max for each combination of BoostRounds, LearningRate, and Phenotype
# Summarize the data: calculating mean, standard error, and min/max
summary_df <- icc_data %>%
  group_by(BoostRounds, LearningRate) %>%
  summarize(
    mean_ICC = mean(Value), # Mean value of ICC - 5 for each combination of BoostRounds, LearningRate
    sd_ICC = sd(Value),
    se_ICC = sd(Value) / sqrt(n()),
    min_ICC = min(Value),
    max_ICC = max(Value),
    .groups = "drop"
  )




# Plot the data with error bars and a trend line
pALL_ICC <- ggplot(summary_df, aes(x = BoostRounds, y = mean_ICC, color = as.factor(LearningRate))) +
  geom_line(linewidth = 0.05, alpha= 0.5, linetype = "dashed") +  # Line connecting the mean values
  
  geom_errorbar(aes(
    ymin = mean_ICC -  se_ICC,
    ymax = mean_ICC +  se_ICC
  ), width = 100) +  
  geom_smooth(data = icc_data, aes(x = as.numeric(BoostRounds), y = Value, color = as.factor(LearningRate), fill =as.factor(LearningRate)) ,
              method = "loess", level = 0.95, span = 1.0,    linetype = "dotted", alpha = 0.1, linewidth = 0.5, show.legend = FALSE ) +
  
  geom_point(data = icc_data, aes(x = BoostRounds, y = Value, color = as.factor(LearningRate)), 
             position = position_jitter(width = 0.1), alpha = 0.8, size =0.75) + 
  # facet_wrap(summary_df$Phenotype <- factor(summary_df$Phenotype, levels = c("GY_1", "TGW_1", "GPC_1"))) + # "HET_1"
  
  labs(title = "ICC vs BoostRounds by LearningRate",
       x = "BoostRounds",
       y = "ICC",
       color = "Learning Rate") +
  theme_minimal()


# FLEISS

# Filter and summarize for Fleiss Kappa metric
FLEISS_data <-Master_df %>% 
  filter(Metric == "Fleiss_Kappa_Prediction_Stability")

summary_df <- FLEISS_data %>%
  group_by(BoostRounds, LearningRate) %>%
  summarize(
    mean_FLEISS = mean(Value),
    sd_FLEISS = sd(Value),
    se_FLEISS = sd(Value) / sqrt(n()),
    min_FLEISS = min(Value),
    max_FLEISS = max(Value),
    .groups = "drop"
  )

# Adjusted Plot 1: Fleiss Kappa vs BoostRounds by LearningRate
pALL_FLEISS <-ggplot(summary_df, aes(x = BoostRounds, y = mean_FLEISS, color = as.factor(LearningRate))) +
  geom_line(linewidth = 0.05, alpha = 0.5, linetype = "dashed") +  
  geom_errorbar(aes(
    ymin = mean_FLEISS - se_FLEISS,
    ymax = mean_FLEISS + se_FLEISS
  ), width = 100) +
  geom_smooth(data = FLEISS_data, 
              aes(x = as.numeric(BoostRounds), y = Value, color = as.factor(LearningRate), 
                  fill = as.factor(LearningRate)),
              method = "loess", level = 0.95, span = 1.0,
              linetype = "dotted", alpha = 0.1, linewidth = 0.5, show.legend = FALSE) +
  geom_point(data = FLEISS_data, 
             aes(x = BoostRounds, y = Value, color = as.factor(LearningRate)),
             position = position_jitter(width = 0.1), alpha = 0.8, size = 0.75) +
  # facet_wrap(summary_df$Phenotype <- factor(summary_df$Phenotype, levels = c("GY_1", "TGW_1", "GPC_1"))) + # "HET_1"
  labs(title = "FLEISS vs BoostRounds by LearningRate",
       x = "BoostRounds",
       y = "FLEISS",
       color = "Learning Rate") +
  theme_minimal()



# Pearson

# Filter the data for ICC metric
pearson_data <- Master_df %>% 
  filter(Metric == "Pearson")

# Summarize the data by computing the mean, min, and max for each combination of BoostRounds, LearningRate, and Phenotype
# Summarize the data: calculating mean, standard error, and min/max
summary_df <- pearson_data %>%
  group_by(BoostRounds, LearningRate) %>%
  summarize(
    mean_pearson = mean(Value),
    sd_pearson = sd(Value),
    se_pearson = sd(Value) / sqrt(n()),
    min_pearson = min(Value),
    max_pearson = max(Value),
    .groups = "drop"
  ) 

summary_df_subset <- pearson_data %>%
  group_by(BoostRounds, LearningRate, Phenotype, Subset) %>%
  summarize(
    mean_pearson_subset = mean(Value),
    sd_pearson_subset = sd(Value),
    se_pearson_subset = sd(Value) / sqrt(n()),
    min_pearson_subset = min(Value),
    max_pearson_subset = max(Value),
    .groups = "drop"
  )




# Plot the data with error bars and a trend line
#### PRESENT ME
pALL_PEARSON <-ggplot(summary_df, aes(x = BoostRounds, y = mean_pearson, color = as.factor(LearningRate))) +
  geom_line(position = position_dodge(width = 50)) +  # Line connecting the mean values
  geom_errorbar(aes(ymin = mean_pearson - se_pearson, ymax = mean_pearson + se_pearson), width = 500, position = position_dodge(width = 50)) +  # Error bars SE
  geom_smooth(method = "loess", se = FALSE, linetype = "dotted", alpha = 0.1, linewidth = 0.5) +  # Overall trend line (line of best fit)
  
  labs(title = "Pearson Correlation vs BoostRounds by LearningRate",
       x = "BoostRounds",
       y = "Pearson Correlation",
       color = "Learning Rate") + #,
  theme_minimal()


# R SQUARED 

# Filter the data for R_squared metric
r2_data <- Master_df %>% 
  filter(Metric == "R_squared")

# Summarize the data by computing the mean, standard error, and min/max for each combination of BoostRounds, LearningRate, and Phenotype
summary_df <- r2_data %>%
  group_by(BoostRounds, LearningRate) %>%
  summarize(
    mean_r2 = mean(Value),
    sd_r2 = sd(Value),
    se_r2 = sd(Value) / sqrt(n()),
    min_r2 = min(Value),
    max_r2 = max(Value),
    .groups = "drop"
  )

# Plot the data with error bars and a trend line
pALL_Rsquared <-ggplot(summary_df, aes(x = BoostRounds, y = mean_r2, color = as.factor(LearningRate))) +
  geom_line() +  # Line connecting the mean R_squared values
  geom_errorbar(aes(ymin = mean_r2 - se_r2, ymax = mean_r2 + se_r2), width = 100) +  # Error bars for SE
  geom_smooth(method = "loess", se = FALSE, linetype = "dotted", alpha = 0.1, linewidth = 0.5) +  # Overall trend line (line of best fit)
  geom_point(data = summary_df, aes(x = BoostRounds, y = mean_r2, color = as.factor(LearningRate))) + 
  # facet_wrap(~Phenotype) +
  labs(title = "R-squared vs BoostRounds by LearningRate",
       x = "BoostRounds",
       y = "R-squared",
       color = "Learning Rate") +
  theme_minimal()


# AUC
# Filter the data for AUC metric
auc_data <- Master_df %>% 
  filter(Metric == "AUC")

# Summarize the data by computing the mean, standard error, and min/max for each combination of BoostRounds, LearningRate, and Phenotype
summary_df <- auc_data %>%
  group_by(BoostRounds, LearningRate) %>%
  summarize(
    mean_auc = mean(Value),
    sd_auc = sd(Value),
    se_auc = sd(Value) / sqrt(n()),
    min_auc = min(Value),
    max_auc = max(Value),
    .groups = "drop"
  )

# Plot the data with error bars and a trend line (Plot 1)
pALL_AUC <-ggplot(summary_df, aes(x = BoostRounds, y = mean_auc, color = as.factor(LearningRate))) +
  geom_line() +  # Line connecting the mean AUC values
  geom_errorbar(aes(ymin = mean_auc - se_auc, ymax = mean_auc + se_auc), width = 100) +  # Error bars SE
  geom_smooth(method = "loess", se = FALSE, linetype = "dotted", alpha = 0.1, linewidth = 0.5) +  # Overall trend line
  geom_point(data = summary_df, aes(x = BoostRounds, y = mean_auc, color = as.factor(LearningRate))) +
  # facet_wrap(summary_df$Phenotype <- factor(summary_df$Phenotype, levels = c("GY_1", "TGW_1", "GPC_1"))) + # "HET_1"
  labs(title = "AUC vs BoostRounds by LearningRate",
       x = "BoostRounds",
       y = "AUC",
       color = "Learning Rate") +
  theme_minimal()


# NDCG at 20 % 

# Filter the data for NDCG_at_20_percent metric
ndcg_data <- Master_df %>% 
  filter(Metric == "NDCG_at_20_percent")

# Summarize the data by computing the mean, standard error, and min/max for each combination of BoostRounds, LearningRate, and Phenotype
summary_df <- ndcg_data %>%
  group_by(BoostRounds, LearningRate) %>%
  summarize(
    mean_ndcg = mean(Value),
    sd_ndcg = sd(Value),
    se_ndcg = sd(Value) / sqrt(n()),
    min_ndcg = min(Value),
    max_ndcg = max(Value),
    .groups = "drop"
  )

# Plot the data with error bars and a trend line (Plot 1)
pALL_NDCG20 <-ggplot(summary_df, aes(x = BoostRounds, y= mean_ndcg, color = as.factor(LearningRate))) +
  geom_line() +  # Line connecting the mean NDCG values
  geom_errorbar(aes(ymin = mean_ndcg - se_ndcg, ymax = mean_ndcg + se_ndcg ), width = 100) +  # Error bars for the range 
  geom_smooth(method = "loess", se = FALSE, linetype = "dotted", alpha = 0.1, linewidth = 0.5) +  # Overall trend line
  geom_point(data = summary_df, aes(x = BoostRounds, y = mean_ndcg, color =as.factor(LearningRate))) +
  #facet_wrap(summary_df$Phenotype <- factor(summary_df$Phenotype, levels = c("GY_1", "TGW_1", "GPC_1"))) + # "HET_1"
  labs(title = "NDCG@20% vs BoostRounds by LearningRate",
       x = "BoostRounds",
       y = "NDCG@20%",
       color = "Learning Rate") +
  theme_minimal()


# Plot all together 
# pALL_ICC + pALL_FLEISS + pALL_PEARSON + pALL_Rsquared + pALL_AUC + pALL_NDCG20
library(gridExtra) 
grid.arrange(pALL_ICC, pALL_FLEISS, pALL_PEARSON, pALL_Rsquared, pALL_AUC, pALL_NDCG20, ncol = 2)
# Save the plots to a PDF file

# ggsave("TWR_Definitive_GY_Super_Dataset_Metrics_Plots.pdf", 
#        arrangeGrob(pALL_ICC, pALL_FLEISS, pALL_PEARSON, pALL_Rsquared, pALL_AUC, pALL_NDCG20, ncol = 2), 
#        width = 12, height = 10)
