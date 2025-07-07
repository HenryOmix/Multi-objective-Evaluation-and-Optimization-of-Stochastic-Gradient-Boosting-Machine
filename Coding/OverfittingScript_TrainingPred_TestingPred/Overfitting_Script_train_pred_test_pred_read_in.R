
# Plot Overfitting: 
# former name: TWR_DEFINITIVE_OVERFITTING_PLOT_all_traits_1.0.0.0.6.hetfix

library(dplyr)
library(tidyverse)  # data wrangling
library(viridis)    # color scales in heatmaps
library(ggplot2)  



#### ROW BIND ALL CORRELATION RESULTS ###
# Combine all correlation results into one data frame
correlation_results_TEST_GY <- read.csv("correlation_results_TEST_GY_36_model.csv")
correlation_results_TRAIN_GY <- read.csv("correlation_results_TRAIN_GY_36_model.csv")
correlation_results_TEST_TGW <- read.csv("correlation_results_TEST_TGW_36_model.csv")
correlation_results_TRAIN_TGW <- read.csv("correlation_results_TRAIN_TGW_36_model.csv")
correlation_results_TEST_GPC <- read.csv("correlation_results_TEST_GPC_36_model.csv")
correlation_results_TRAIN_GPC <- read.csv("correlation_results_TRAIN_GPC_36_model.csv")
correlation_results_TEST_HET <- read.csv("correlation_results_TEST_HET_36_model.csv")
correlation_results_TRAIN_HET <- read.csv("correlation_results_TRAIN_HET_36_model.csv")

# Combine all correlation results into one data frame
all_correlation_results_TEST <- bind_rows(
  correlation_results_TEST_GY %>% mutate(trait = "GY"),
  correlation_results_TEST_TGW %>% mutate(trait  = "TGW"),
correlation_results_TEST_GPC %>% mutate(trait = "GPC"),
  correlation_results_TEST_HET %>% mutate(trait = "HET")
)
all_correlation_results_TEST <- all_correlation_results_TEST %>% mutate(dataset = "Test")

all_correlation_results_TRAIN <- bind_rows(
  correlation_results_TRAIN_GY %>% mutate(trait =  "GY"),
  correlation_results_TRAIN_TGW %>% mutate(trait= "TGW"),
 correlation_results_TRAIN_GPC %>% mutate(trait = "GPC"),
  correlation_results_TRAIN_HET %>% mutate(trait = "HET")
)
all_correlation_results_TRAIN <- all_correlation_results_TRAIN %>% mutate(dataset = "Train")


# Combine the data frames
all_df <- bind_rows(all_correlation_results_TEST, all_correlation_results_TRAIN)

# Aggregate by eta (learning rate), nround (boost rounds), and dataset
agg_df <- all_df %>%
  group_by(eta, nround, dataset) %>%
  summarise(mean_correlation = mean(correlation), .groups = 'drop')


# 2. Faceted Line Plots (Train vs. Test by eta)


p1 <- ggplot(agg_df, aes(x = nround, y = mean_correlation, color = dataset)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ eta, scales = "free", ncol  = 2) +
  labs(title = "Accuracy vs Boost Rounds by Learning Rate",
       x = "Boost Rounds (nround)",
       y = "Mean Correlation",
       color =  "Dataset") +
ylim(0, 1) +
  theme_minimal()

print(p1)