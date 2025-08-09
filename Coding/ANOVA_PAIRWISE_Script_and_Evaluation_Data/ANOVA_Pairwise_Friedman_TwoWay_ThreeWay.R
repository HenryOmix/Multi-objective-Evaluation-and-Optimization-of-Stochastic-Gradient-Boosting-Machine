# this script runs ANOVA, emmeans estimates and post hoc pairwise tests
## NOTE!! "Phenotype" is technically: "Trait"
# original name: "TWR_Definitive_STATISTICS_all_traits_1.0.0.0.9.6.1"

install.packages("pbkrtest")
install.packages("lmerTest") 
install.packages("emmeans")
install.packages("effectsize")   

library(dplyr)
library(tidyr)
library(lme4)
library(emmeans)
library(pbkrtest)
library(lmerTest) 
library(effectsize)             

# 31.07.25
options(contrasts = c("contr.sum", "contr.poly"))
emm_options(pbkrtest.limit = 7200) #  

# more info on type III , partial eta squared: https://easystats.github.io/effectsize/articles/anovaES.html


# read & Combine All Trait-Specific CSVs
GY_Master_metrics<- read.csv("Definitive_GY_36_GBM.fit_MASTER_METRICS_DF_27.03.2025.csv")
TGW_Master_metrics<-  read.csv("Definitive_TGW_36_GBM.fit_MASTER_METRICS_DF_27.03.2025.csv")
GPC_Master_metrics<- read.csv("Definitive_GPC_36_GBM.fit_MASTER_METRICS_DF_27.03.2025.csv")
HET_Master_metrics <- read.csv("Definitive_HET_36_GBM.fit_MASTER_METRICS_DF_27.03.2025.csv")

#remove index colums - unwanted
GY_Master_metrics$X <- NULL
TGW_Master_metrics$X <- NULL
GPC_Master_metrics$X<- NULL
HET_Master_metrics$X<- NULL

#combine all trait metrics into one master dataframe
Master_df <- rbind(
  GY_Master_metrics,
  TGW_Master_metrics,
  GPC_Master_metrics,
  HET_Master_metrics
)

# 31.07.25
Master_df <- Master_df %>% 
  mutate(
    LearningRate = factor(LearningRate),
    BoostRounds  = factor(BoostRounds),
    Phenotype    = factor(Phenotype)
  )

lmm_metrics <- c("Pearson", "R_squared", "AUC", "NDCG_at_20_percent")
stability_metrics <- c("ICC_Prediction_Stability", "Fleiss_Kappa_Prediction_Stability")

thesis_metrics <- c("Pearson", "R_squared", "AUC", "NDCG_at_20_percent", "ICC_Prediction_Stability", "Fleiss_Kappa_Prediction_Stability")

# Per Phenotype LMM(Metrics with 50 obs/param)
#####
Phenotypes <- unique(Master_df$Phenotype) 

for (Phenotype_i in Phenotypes) {
  
  # Filter data for current phenotype
  Phenotype_data <- Master_df %>% filter(Phenotype ==Phenotype_i)
  
  for (metric in lmm_metrics) {
    
    cat("Trait:", Phenotype_i, "_____", "Metric:", metric) 
    
    df_sub <- Phenotype_data %>% filter(Metric==metric)
    
    #### MODEL: Using LearningRate * BoostRounds 
    lmer_lr_br <-lmer(Value ~ LearningRate * BoostRounds + (1|Subset) + (1|Subset:Iteration), data= df_sub)
    print(summary(lmer_lr_br))
    
    # qq Plot for Residuals -normality
    qqnorm(resid(lmer_lr_br), main =paste("Q-Q Plot: Model 2 for", Phenotype_i, "-", metric))
    qqline(resid(lmer_lr_br))
    
    # Residual vs. Fitted Plot -homoscedasticity
    plot(fitted(lmer_lr_br), resid(lmer_lr_br),
         main= paste("Residuals vs Fitted: Model 2 for", Phenotype_i, "-", metric),
         xlab = "Fitted Values",
         ylab= "Residuals")
    abline(h= 0, col = "red")
    
    anova_res2 <- anova(lmer_lr_br, type = 3)  # # 31.07.25
    anova_df2 <- as.data.frame(anova_res2)
    
    ### EFFECT SIZE start
    eta2_df <- effectsize::eta_squared(anova_res2, partial = TRUE)
    anova_df2$Partial_Eta_Sq <- eta2_df$Eta2_partial[match(rownames(anova_df2), eta2_df$Parameter)]
    ### EFFECT SIZE end
    
    write.csv(anova_df2, file= paste0("ANOVA_Model2_",Phenotype_i,"_",metric,".csv"), row.names = TRUE)
    
    lr_levels <- sort(unique(df_sub$LearningRate))
    br_levels <- sort(unique(df_sub$BoostRounds))
    emm_lr_br <-emmeans(lmer_lr_br, ~ LearningRate * BoostRounds,
                        at= list(LearningRate = lr_levels, BoostRounds = br_levels))
    emm_lr_br_df <-as.data.frame(summary(emm_lr_br))
    write.csv(emm_lr_br_df, file = paste0("EMMEANS_Model2_",Phenotype_i,"_",metric,".csv"), row.names= FALSE)
    
    pairwise_lr_br<- contrast(emm_lr_br, method = "pairwise", adjust = "tukey")
    pairwise_lr_br_df <- as.data.frame(summary(pairwise_lr_br))
    write.csv(pairwise_lr_br_df, file= paste0("Pairwise_Model2_",Phenotype_i,"_",metric,".csv"), row.names= FALSE)
  }
}


# Friedman Test for Stability Metrics (5 obs/param)
for (Phenotype_i in Phenotypes) {
  
  Phenotype_data <-Master_df %>% filter(Phenotype== Phenotype_i)
  
  for (metric in stability_metrics) {
    
    cat("Trait:", Phenotype_i, "_____", "Metric:", metric) 
    
    df_sub <- Phenotype_data %>% filter(Metric == metric)
    
    ftest <- friedman.test(Value ~ ParameterLabel | Subset, data = df_sub)
    print(ftest)
    
    # save Friedman test results to CSV 
    friedman_out <- data.frame(
      Phenotype = Phenotype_i,
      Metric= metric,
      ChiSq= ftest$statistic,
      df = ftest$parameter,
      p.value = ftest$p.value
    )
    
    ### Kendal W ~ effect size start
    n_blocks  <- length(unique(df_sub$Subset))
    k_groups  <- length(unique(df_sub$ParameterLabel))
    friedman_out$Kendall_W <- as.numeric(ftest$statistic) / (n_blocks * (k_groups - 1))
    ### Kendal W ~ effect size end
    
    write.csv(friedman_out, file = paste0("Friedman_",Phenotype_i,"_",metric,".csv"), row.names= FALSE)
  }
}





## 6.01) Overall LMM Including All Phenotypes with Three-Way Interaction
######
for (metric in c(lmm_metrics)) {
  
  df_sub<- Master_df %>% filter(Metric== metric)
  
  cat("Across traits three-way:", metric) 
  
  overall_lmer_601 <- lmer(
    Value ~ Phenotype + LearningRate * BoostRounds +
      Phenotype:LearningRate + Phenotype:BoostRounds +
      Phenotype:LearningRate:BoostRounds +
      (1|Phenotype:Subset) + (1|Phenotype:Subset:Iteration),
    data= df_sub
  )
  print(summary(overall_lmer_601))
  
  #qq plot for Residuals
  qqnorm(resid(overall_lmer_601), main = paste("Q-Q Plot: Overall Model 6.01 for", metric))
  qqline(resid(overall_lmer_601))
  
  # Residual vs Fitted Plot
  plot(fitted(overall_lmer_601), resid(overall_lmer_601),
       main = paste("Residuals vs Fitted: Overall Model 6.01 for", metric),
       xlab= "Fitted Values",
       ylab = "Residuals")
  abline(h= 0, col = "red")
  
  overall_anova <- anova(overall_lmer_601, type = 3) # 31.07.25
  overall_anova_df <- as.data.frame(overall_anova)
  
  ### EFFECT SIZE start
  eta2_df <- effectsize::eta_squared(overall_anova, partial = TRUE)
  overall_anova_df$Partial_Eta_Sq <- eta2_df$Eta2_partial[match(rownames(overall_anova_df), eta2_df$Parameter)]
  ### EFFECT SIZE end
  
  write.csv(overall_anova_df, file = paste0("Overall_ANOVA_601_", metric, ".csv"), row.names = TRUE)
  
  #EMMEANS for the LearningRate * BoostRounds interaction 
  lr_levels_overall <- sort(unique(df_sub$LearningRate))
  br_levels_overall <- sort(unique(df_sub$BoostRounds))
  
  emm_overall <- emmeans(overall_lmer_601, ~ LearningRate * BoostRounds | Phenotype,
                         at= list(LearningRate= lr_levels_overall, BoostRounds= br_levels_overall))
  emm_overall_df <- as.data.frame(summary(emm_overall))
  write.csv(emm_overall_df, file = paste0("EMMEANS_by_trait_overallmodel_601_", metric, ".csv"), row.names = FALSE)
  
  pairwise_overall <- contrast(emm_overall, method= "pairwise", adjust= "tukey")
  pairwise_overall_df <- as.data.frame(summary(pairwise_overall))
  
  write.csv(pairwise_overall_df, file = paste0("Pairwise_by_trait_overallmodel_601_", metric, ".csv"), row.names= FALSE)
}


###########
## 6.21) Overall LMM for Stability Metrics with Three-Way Interaction
####
for (metric in stability_metrics) {
  
  cat("Across traits three-way stability:  ", metric) 
  
  df_sub <- Master_df %>% filter(Metric== metric)
  
  overall_lmer_stab_621 <- lmer(
    Value ~ Phenotype + LearningRate * BoostRounds +
      Phenotype:LearningRate + Phenotype:BoostRounds +
      Phenotype:LearningRate:BoostRounds +
      (1|Phenotype:Subset),
    data= df_sub
  )
  print(summary(overall_lmer_stab_621))
  
  # qq Plot for Residuals
  qqnorm(resid(overall_lmer_stab_621), main = paste("Q-Q Plot: Overall Model 6.21 for", metric))
  qqline(resid(overall_lmer_stab_621))
  
  # Residual vs Fitted Plot
  plot(fitted(overall_lmer_stab_621), resid(overall_lmer_stab_621),
       main= paste("Residuals vs Fitted: Overall Model 6.21 for", metric),
       xlab = "Fitted Values",
       ylab = "Residuals")
  abline(h= 0, col= "red")
  
  overall_anova_stab <- anova(overall_lmer_stab_621, type = 3) # 31.07.25
  overall_anova_stab_df <- as.data.frame(overall_anova_stab)
  
  ### EFFECT SIZE start
  eta2_df <- effectsize::eta_squared(overall_anova_stab, partial = TRUE)
  overall_anova_stab_df$Partial_Eta_Sq <- eta2_df$Eta2_partial[match(rownames(overall_anova_stab_df), eta2_df$Parameter)]
  ### EFFECT SIZE end
  
  write.csv(overall_anova_stab_df, file = paste0("Overall_ANOVA_621_", metric, ".csv"), row.names = TRUE)
  
  # EMMEANS for the LearningRate * BoostRounds interaction
  lr_levels_overall <- sort(unique(df_sub$LearningRate))
  br_levels_overall<- sort(unique(df_sub$BoostRounds))
  
  emm_overall <- emmeans(overall_lmer_stab_621, ~ LearningRate * BoostRounds | Phenotype,
                         at= list(LearningRate= lr_levels_overall, BoostRounds = br_levels_overall))
  emm_overall_df <- as.data.frame(summary(emm_overall))
  write.csv(emm_overall_df, file = paste0("EMMEANS_by_trait_overallmodel_621_", metric, ".csv"), row.names = FALSE)
  
  pairwise_overall <-  contrast(emm_overall, method = "pairwise", adjust= "tukey")
  pairwise_overall_df<- as.data.frame(summary(pairwise_overall))
  
  write.csv(pairwise_overall_df, file= paste0("Pairwise_by_trait_overallmodel_621_", metric, ".csv"), row.names = FALSE)
}




## 6.3) Adjusted Overall LMM for lmm_metrics (Pooled Across Phenotypes)
for (metric in lmm_metrics) {
  
  df_sub <- Master_df %>% filter(Metric == metric)
  
  cat("Aggregated traits two-way performance:  ", metric) 
  
  overall_lmer_adj <- lmer(
    Value ~ LearningRate * BoostRounds + Phenotype + (1|Phenotype:Subset) + (1|Phenotype:Subset:Iteration),
    data = df_sub
  )
  print(summary(overall_lmer_adj))
  
  # qq Plot for Residuals
  qqnorm(resid(overall_lmer_adj), main = paste("Q-Q Plot: Adjusted Overall Model for", metric))
  qqline(resid(overall_lmer_adj))
  
  # Residual vs. Fitted Plot
  plot(fitted(overall_lmer_adj), resid(overall_lmer_adj),
       main = paste("Residuals vs Fitted: Adjusted Overall Model for", metric),
       xlab= "Fitted Values",
       ylab= "Residuals")
  abline(h= 0, col = "red")
  
  overall_anova_adj <- anova(overall_lmer_adj, type = 3) # 31.07.25
  overall_anova_adj_df <- as.data.frame(overall_anova_adj)
  
  ### EFFECT SIZE start
  eta2_df <- effectsize::eta_squared(overall_anova_adj, partial = TRUE)
  overall_anova_adj_df$Partial_Eta_Sq <- eta2_df$Eta2_partial[match(rownames(overall_anova_adj_df), eta2_df$Parameter)]
  ### EFFECT SIZE end
  
  write.csv(overall_anova_adj_df, file = paste0("Overall_ANOVA_Adj_", metric, ".csv"), row.names = TRUE)
  
  #perform emmeans post-hoc analysis
  lr_levels_overall <- sort(unique(df_sub$LearningRate))
  br_levels_overall <- sort(unique(df_sub$BoostRounds))
  
  emm_overall_adj <- emmeans(overall_lmer_adj, ~ LearningRate * BoostRounds,
                             at= list(LearningRate= lr_levels_overall, BoostRounds = br_levels_overall))
  emm_overall_adj_df <- as.data.frame(summary(emm_overall_adj))
  write.csv(emm_overall_adj_df, file = paste0("EMMEANS_overallmodel_Adj_", metric, ".csv"), row.names= FALSE)
  
  pairwise_overall_adj <-contrast(emm_overall_adj, method= "pairwise", adjust = "tukey")
  pairwise_overall_adj_df <- as.data.frame(summary(pairwise_overall_adj))
  write.csv(pairwise_overall_adj_df, file = paste0("Pairwise_overallmodel_Adj_", metric, ".csv"), row.names= FALSE)
}


## 6.4) Adjusted Overall LMM for Stability Metrics (aggregated/pooled Across Phenotypes)
for (metric in stability_metrics) { 
  
  df_sub <- Master_df %>% filter(Metric== metric)
  
  cat("Aggregated traits two-way performance stability:  ", metric) 
  
  overall_lmer_stab_adj <- lmer(
    Value ~ LearningRate * BoostRounds + Phenotype + (1|Subset),
    data = df_sub
  )
  print(summary(overall_lmer_stab_adj))
  
  # qq Plot for Residuals
  qqnorm(resid(overall_lmer_stab_adj), main = paste("Q-Q Plot: Adjusted Overall Model for", metric))
  qqline(resid(overall_lmer_stab_adj))
  
  # Residual vs. Fitted Plot
  plot(fitted(overall_lmer_stab_adj), resid(overall_lmer_stab_adj),
       main= paste("Residuals vs Fitted: Adjusted Overall Model for", metric),
       xlab= "Fitted Values",
       ylab = "Residuals")
  abline(h = 0, col = "red")
  
  overall_anova_stab_adj <- anova(overall_lmer_stab_adj, type = 3) # 31.07.25
  overall_anova_stab_adj_df <- as.data.frame(overall_anova_stab_adj)
  
  ### EFFECT SIZE start
  eta2_df <- effectsize::eta_squared(overall_anova_stab_adj, partial = TRUE)
  overall_anova_stab_adj_df$Partial_Eta_Sq <- eta2_df$Eta2_partial[match(rownames(overall_anova_stab_adj_df), eta2_df$Parameter)]
  ### EFFECT SIZE end
  
  write.csv(overall_anova_stab_adj_df, file= paste0("Overall_ANOVA_Stability_Adj_", metric, ".csv"), row.names= TRUE)
  
  #perform emmeans post-hoc analysis
  lr_levels_overall <- sort(unique(df_sub$LearningRate))
  br_levels_overall <- sort(unique(df_sub$BoostRounds))
  
  emm_overall_adj <- emmeans(overall_lmer_stab_adj, ~ LearningRate * BoostRounds,
                             at = list(LearningRate = lr_levels_overall, BoostRounds = br_levels_overall))
  emm_overall_adj_df <- as.data.frame(summary(emm_overall_adj))
  write.csv(emm_overall_adj_df, file= paste0("EMMEANS_overallmodel_Adj_", metric, ".csv"), row.names = FALSE)
  
  pairwise_overall_adj <- contrast(emm_overall_adj, method= "pairwise", adjust = "tukey")
  pairwise_overall_adj_df<- as.data.frame(summary(pairwise_overall_adj))
  write.csv(pairwise_overall_adj_df, file = paste0("Pairwise_overallmodel_Adj_", metric, ".csv"), row.names= FALSE)
}
