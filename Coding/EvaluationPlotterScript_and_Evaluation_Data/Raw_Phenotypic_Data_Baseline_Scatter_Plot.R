library(lme4)     
library(dplyr)    
library(tidyr)     
install.packages("sommer")
library(sommer)


folder_name <-   "AllTraitsMasterLoopTrainingTestingMissRangerImputed"


# Construct file paths for the fixed training and testing sets
fixed_train_file_path <- paste0(folder_name, "/training_set_", 1, ".csv")
fixed_test_file_path   <- paste0(folder_name, "/testing_set_",1, ".csv")

# Load the subset data
fixed_training_set <- read.csv(fixed_train_file_path)
fixed_testing_set  <-read.csv(fixed_test_file_path)

nrow(fixed_training_set)
nrow(fixed_testing_set)

# Row bind training and testing sets
Total_Dataset <- rbind(fixed_training_set, fixed_testing_set)

nrow(Total_Dataset)

# Define the traits for which paired year values exist
traits <- c("GY", "TGW",  "GPC", "HET") #"YR", 

# Create a dataframe to store the results
accuracy_results <- data.frame(Trait = traits, Spearman = NA, R_squared = NA)

# Loop through each trait pair (e.g., GY_1 vs GY_2) to calculate accuracy measures
for (trait in traits) {
  trait_year1 <- paste0(trait, "_1")
  trait_year2 <- paste0(trait, "_2")
  
  # Calculate the Spearman correlation between year 1 and year 2
  #spearman_corr <- cor(Total_Dataset[[trait_year1]], Total_Dataset[[trait_year2]],
  #                     method = "spearman", use = "complete.obs")
  
  spearman_corr <-  cor.test(Total_Dataset[[trait_year1]], Total_Dataset[[trait_year2]],
                             method = "spearman", exact = FALSE)$estimate
  
  # Calculate pseudo R^2 
  # (treating year 2 values as predictions for year 1 observations)
  ss_res <- sum((Total_Dataset[[trait_year1]] - Total_Dataset[[trait_year2]])^2, na.rm = TRUE)
  ss_tot <- sum((Total_Dataset[[trait_year1]] - mean(Total_Dataset[[trait_year1]], na.rm = TRUE))^2, na.rm = TRUE)
  r_squared <- 1 - (ss_res / ss_tot)
  
  # Save the results in the dataframe
  accuracy_results[accuracy_results$Trait == trait, "Spearman"] <- spearman_corr
  accuracy_results[accuracy_results$Trait == trait, "R_squared"] <- r_squared
}

# Display the results
print(accuracy_results)




#### PLOTTING ALL 5 TRAITS ####

# Set up the plotting area: 1 row and 5 columns
par(mfrow = c(1, 4), mar = c(4, 4, 2, 1))  

# Loop through each trait pair (e.g., GY_1 vs GY_2) and create scatter plots
for (trait in traits) {
  trait_year1<- paste0(trait, "_1")
  trait_year2 <- paste0(trait, "_2")
  
  # Identify the maximum and minimum values from both Year 1 and Year 2 columns for this trait
  max_val <- max(Total_Dataset[[trait_year1]], Total_Dataset[[trait_year2]], na.rm = TRUE)
  min_val <- min(Total_Dataset[[trait_year1]], Total_Dataset[[trait_year2]], na.rm =TRUE)
  
  # Create a scatter plot with x and y limits set from min_val to max_val
  plot(Total_Dataset[[trait_year1]], Total_Dataset[[trait_year2]],
       main = trait,
       xlab = paste(trait, "Year 1"),
       ylab = paste(trait, "Year 2"),
     xlim = c(min_val, max_val),
       ylim = c(min_val, max_val),
       pch = 19, col = "blue", cex = 0.25)
  
  # Add the identity (45-degree) line for reference
  abline(0, 1, lty = 2, col = rgb(1,0,0,0.5)) 
  
  # TASK 1:  secondary line – the regression (best-fit) line
  fit <- lm(Total_Dataset[[trait_year2]] ~ Total_Dataset[[trait_year1]])
  abline(fit, col = "darkgreen", lty = 2)
  
  # Calculate statistics for annotation
  
  ## Year 1 statistics
  mean1 <- mean(Total_Dataset[[trait_year1]], na.rm = TRUE)
  sd1 <-sd(Total_Dataset[[trait_year1]], na.rm = TRUE)
  n1 <- sum(!is.na(Total_Dataset[[trait_year1]]))
  se1 <- sd1 / sqrt(n1)
  
  ## Year 2 statistics
  mean2 <- mean(Total_Dataset[[trait_year2]],na.rm= TRUE)
  sd2 <-sd(Total_Dataset[[trait_year2]], na.rm =TRUE)
  n2 <- sum(!is.na(Total_Dataset[[trait_year2]]))
  se2<- sd2 / sqrt(n2)
  
  ## Calculate Spearman correlation
  spearman_corr <-cor.test(Total_Dataset[[trait_year1]], Total_Dataset[[trait_year2]], 
                            method = "pearson", exact = FALSE)$estimate
  
  ## Calculate R^2 (treating Year 2 as prediction for Year 1)
  ss_res <- sum((Total_Dataset[[trait_year1]] - Total_Dataset[[trait_year2]])^2, na.rm = TRUE)
  ss_tot <- sum((Total_Dataset[[trait_year1]] - mean(Total_Dataset[[trait_year1]], na.rm = TRUE))^2, na.rm = TRUE)
  r_squared <- 1 -(ss_res / ss_tot)
  
  # Create annotation text
  stats_text <- paste("Year 1: Mean=", round(mean1,2), ", SD=", round(sd1,2), ", SE=", round(se1,2),
                      "\nYear 2: Mean=", round(mean2,2), ", SD=", round(sd2,2), ", SE=", round(se2,2),
                      "\nPearson=", round(spearman_corr,2),", R²=", round(r_squared,2))
  
  # Annotate the plot at the bottom left
  # Positioning: a little offset from the minimum value
  text(x = min_val + 0.05* (max_val - min_val), 
      y = max_val- 0.05*(max_val -min_val),
       labels = stats_text, pos = 4, cex = 1.25)   
}
