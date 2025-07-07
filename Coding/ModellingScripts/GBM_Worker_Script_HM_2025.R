# Former r file name, "DEFINITIVE_MT_gbm.fit_gbm()_TASK1_GradientBoost_DIRECT_rmFeat_P_matrices_parallel_MODEL_GRID_3.0_MasterRun_1.4" 
# new name: GBM_Worker_Script_HM_2025.R



run_experiment <- function(
    tune_grid_arg,         # pass in the tuneGrid data.frame
 
base_path = "C:/Users/hmunr/Downloads",
results_folder_core_name = "GBM_package_final_model_training_Environment_04.03.2025",
folder_name = "AllTraitsMasterLoopTrainingTestingMissRangerImputed",
iterations_fixed = 10  ,  # Number of iterations per subset for stability
num_subsets = 10,          # Number of subsets to evaluate - assisting with generalizability and stability
FEATURE_SELECTION_THRESHOLD = 0.05,
PHENOTYPE_SELECTION_THRESHOLD = 0.25,

# # Number of cores to use
# num_cores = round(detectCores() * 0.8),  # Using 80% cores # not using it since we parallelize cores for each worker through the tune grid (see intitiator script)

param_label, 

# Define the phenotype of interest
phenotype_of_interest,
phenotype_of_interest_year2, 

timestamp = format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
learning_rate, 
boost_rounds
) {

    # Append the timestamp to the folder name for uniqueness
    results_folder_name <- paste0( param_label, "_", results_folder_core_name, "_", timestamp)
    
    # Update the results folder path with the new unique folder name
    results_folder <- file.path(base_path, results_folder_name)
  
  # Create the folder if it doesn't exist
  dir.create(results_folder, recursive = TRUE)

# Vectors to store RMSE results
results_df <- data.frame(
  SNP = character(),
  Baseline_RMSE = numeric(),
  Permuted_RMSE = numeric(),
  Delta_RMSE = numeric(),
  Subset_no = integer(),
  Iteration_no = integer(),
  stringsAsFactors = FALSE
)

# Initialise vector to store models for each subset
models_gbm_fixed <- vector("list", num_subsets)

# Vectors to store ICC results for predictions and feature importance
icc_pred_values_gbm_fixed   <- numeric(num_subsets)
fleiss_pred_values_fixed <- numeric(num_subsets)

# Initialize results containers by subset and iteration
shapiro_results_gbm_fixed <- vector("list", num_subsets)
spearman_cor_results_gbm_fixed <- vector("list", num_subsets)
pearson_cor_results_gbm_fixed  <- vector("list", num_subsets)
spearman_accuracies_gbm_fixed  <- vector("list", num_subsets)
pearson_accuracies_gbm_fixed   <- vector("list", num_subsets)
r_squared_values_gbm_fixed     <- vector("list", num_subsets)
mse_values_gbm_fixed           <- vector("list", num_subsets)
predicted_vs_observed_gbm_fixed <- vector("list", num_subsets)
predicted_vs_observed_TRAINING_gbm_fixed <- vector("list", num_subsets)
ndcg_results_gbm_fixed         <- vector("list", num_subsets)

# for addressing hyperparemter interactions with expand.grid
best_tune_results <- vector("list", num_subsets)

# Classification metrics
specificity_gbm_fixed <- vector("list", num_subsets)
sensitivity_gbm_fixed <- vector("list", num_subsets)
precision_gbm_fixed   <- vector("list", num_subsets)
f1_score_gbm_fixed     <- vector("list", num_subsets)
auc_gbm_fixed          <- vector("list", num_subsets)

# Inbuilt feature importance
permutation_importance_gbm_fixed <- vector("list", num_subsets)

### Start measuring time ###
start_time <- Sys.time()

# Function for fleiss kappa phenotype stability
assign_phenotype_selection <- function(mat, selected_proportion = PHENOTYPE_SELECTION_THRESHOLD) {
  # For each column (i.e. each iteration), compute the cutoff such that the top selected_proportion values are "selected"
  selection_mat <- apply(mat, 2, function(col) {  # Column-wise operation
    cutoff <- quantile(col, probs = 1 - selected_proportion, na.rm = TRUE)
    ifelse(col >= cutoff, "selected", "not_selected")
  })
  # Ensure that row names are preserved
  rownames(selection_mat) <- rownames(mat)
  return(selection_mat)
}



# Now, initialize a list for each subset's results (iterations contained within)
for (i in 1:num_subsets) {
  shapiro_results_gbm_fixed[[i]] <- vector("list", iterations_fixed)
  spearman_cor_results_gbm_fixed[[i]] <- vector("list", iterations_fixed)
  pearson_cor_results_gbm_fixed[[i]] <- vector("list", iterations_fixed)
  predicted_vs_observed_gbm_fixed[[i]] <- vector("list", iterations_fixed)
  predicted_vs_observed_TRAINING_gbm_fixed[[i]] <- vector("list", iterations_fixed)
  ndcg_results_gbm_fixed[[i]] <- vector("list", iterations_fixed)
  permutation_importance_gbm_fixed[[i]] <- vector("list", iterations_fixed)
  
  # Initialize other metrics as well
  spearman_accuracies_gbm_fixed[[i]]  <- numeric(iterations_fixed)
  pearson_accuracies_gbm_fixed[[i]]   <- numeric(iterations_fixed)
  r_squared_values_gbm_fixed[[i]]     <- numeric(iterations_fixed)
  mse_values_gbm_fixed[[i]]           <- numeric(iterations_fixed)
  specificity_gbm_fixed[[i]]          <- numeric(iterations_fixed)
  sensitivity_gbm_fixed[[i]]          <- numeric(iterations_fixed)
  precision_gbm_fixed[[i]]            <- numeric(iterations_fixed)
  f1_score_gbm_fixed[[i]]             <- numeric(iterations_fixed)
  auc_gbm_fixed[[i]]                  <- numeric(iterations_fixed)
  
  # Initialize the best tune results
  best_tune_results[[i]] <- vector("list", iterations_fixed)
  
  # Save model for each iteration
  models_gbm_fixed[[i]] <- vector("list", iterations_fixed)
}


# Model Stability --> Intraclass correlation Coefficient (ICC)
ICC_matrix_creation_based_on_testing <- paste0(folder_name, "/testing_set_", 1, ".csv")
fixed_testing_set <- read.csv(ICC_matrix_creation_based_on_testing)
num_test_subjects <- nrow(fixed_testing_set)

# Create a matrix to store predictions from the master loop (random subsets)
pred_matrix_gbm_fixed <- matrix(NA, nrow = num_test_subjects, ncol = iterations_fixed)
colnames(pred_matrix_gbm_fixed) <- paste0("Iter_", 1:iterations_fixed)


# Traits to analyze
all_traits <- c("GY_1", "GY_2", "TGW_1", "TGW_2", "YR_1", "YR_2", "GPC_1", "GPC_2", "HET_1", "HET_2")

line_name_variable <- "line_name"

#remove number or underscore from phenotype of interest - only for filename outputs
phenotype_of_interest_config_filename <- gsub("[0-9_]", "", phenotype_of_interest)

# Identify traits to remove (all except the phenotype of interest)
traits_to_remove_training <- setdiff(all_traits, phenotype_of_interest)
traits_to_remove_testing  <- setdiff(all_traits, phenotype_of_interest_year2)

# Master Loop: XGBoost GBM over multiple random subsets
for (fixed_subset_index in 1:num_subsets) {
  
  cat("\n*** Analysis for subset:", fixed_subset_index, "***\n")
  
  # Construct file paths for the fixed training and testing sets
  fixed_train_file_path <- paste0(folder_name, "/training_set_", fixed_subset_index, ".csv")
  fixed_test_file_path  <- paste0(folder_name, "/testing_set_",  fixed_subset_index, ".csv")
  
  # Load the subset data
  fixed_training_set <- read.csv(fixed_train_file_path)
  fixed_testing_set  <- read.csv(fixed_test_file_path)
  
  # Subset for fast runs (FORGET ME NOT DURING ACTUAL RUN!)
   # fixed_training_set <- fixed_training_set[,1:100 ]
   # fixed_testing_set  <- fixed_testing_set[,1:100 ]
  
  # Prepare a matrix to store predictions for this subset across iterations
  pred_matrix_gbm_fixed <- matrix(NA, nrow = nrow(fixed_testing_set), ncol = iterations_fixed)
  colnames(pred_matrix_gbm_fixed) <- paste0("Iter_", 1:iterations_fixed)
  
  # ADDITION: For training data (only for overfitting insights to consider.. not used yet)
  pred_matrix_TRAINING_gbm_fixed <- matrix(NA, nrow = nrow(fixed_training_set), ncol = iterations_fixed)
  colnames(pred_matrix_TRAINING_gbm_fixed) <- paste0("Iter_", 1:iterations_fixed)
  
  # Progress bar for the inner loop (stability iterations)
  pb <- txtProgressBar(min = 0, max = iterations_fixed, style = 3)
  
  # Inner loop: train and evaluate on the same subset multiple times
  for (stab_i in 1:iterations_fixed) {
    setTxtProgressBar(pb, stab_i)
    
    # Setting seed, need to keep both subset and iteration in mind, multiplying by large number to avoid overlap
    set.seed(123 + fixed_subset_index * 1000 + stab_i)
    
    # Remove unwanted trait columns
    fixed_training_sub <- fixed_training_set[, !(names(fixed_training_set) %in% traits_to_remove_training)]
    fixed_testing_sub  <- fixed_testing_set[, !(names(fixed_testing_set) %in% traits_to_remove_testing)]
    
    # Remove line_name column
    fixed_training_sub <- fixed_training_sub[, !(names(fixed_training_sub) %in% line_name_variable)]
    fixed_testing_sub  <- fixed_testing_sub[, !(names(fixed_testing_sub) %in% line_name_variable)]
    
    # Rename phenotype to "PHENOTYPE"
    names(fixed_training_sub)[names(fixed_training_sub) == phenotype_of_interest] <- "PHENOTYPE"
    names(fixed_testing_sub)[names(fixed_testing_sub) == phenotype_of_interest_year2] <- "PHENOTYPE"

    
    # Define predictor columns (SNPs)
    predictor_cols <- setdiff(names(fixed_testing_sub), c("PHENOTYPE", line_name_variable))
    
    # Convert data frames to matrices (forcing a matrix output with drop = FALSE)
    x_train_fixed <- as.matrix(fixed_training_sub[, !(names(fixed_training_sub) %in% "PHENOTYPE"), drop = FALSE])
    y_train_fixed <- fixed_training_sub$PHENOTYPE   
    x_test_fixed  <- as.matrix(fixed_testing_sub[, !(names(fixed_testing_sub) %in% "PHENOTYPE"), drop = FALSE])
    
    # Convert the one-row tuning grid to a parameter list for gbm.fit (using the same parameter names as before)
    nrounds_val <- as.numeric(tune_grid_arg$ntrees)
    
    # --- Train gbm using gbm.fit directly ---
    gbm_model_fixed <- gbm.fit(
      x = x_train_fixed,
      y = y_train_fixed,
      distribution = tune_grid_arg$distribution,      
      n.trees = as.numeric(tune_grid_arg$ntrees),         
      interaction.depth = tune_grid_arg$interaction.depth,  
      shrinkage = as.numeric(tune_grid_arg$shrinkage),    
      n.minobsinnode = tune_grid_arg$n.minobsinnode,      
      bag.fraction = tune_grid_arg$bag.fraction,        
      verbose = tune_grid_arg$verbose             
    )
    
    # Store the model using the same naming convention
    models_gbm_fixed[[fixed_subset_index]][[stab_i]] <- gbm_model_fixed
    
    model_call_text <- paste(deparse(gbm_model_fixed$call), collapse = "\n")
    
    # For logging purposes, store the tuning parameters
    best_params <- tune_grid_arg
    best_tune_results[[fixed_subset_index]][[stab_i]] <- best_params
    
    # Predict on the testing set (convert newdata to a matrix)
    pred_vals_fixed <- predict(
      gbm_model_fixed,
      newdata = x_test_fixed,
      n.trees = nrounds_val
    )
    
    # Store predictions as before
    pred_matrix_gbm_fixed[, stab_i] <- pred_vals_fixed
    predicted_vs_observed_gbm_fixed[[fixed_subset_index]][[stab_i]] <- data.frame(
      Observed = fixed_testing_sub$PHENOTYPE,
      Predicted = pred_vals_fixed
    )
    
    # Predict on the training set (for overfitting checks)
    pred_vals_train <- predict(
      gbm_model_fixed,
      newdata = x_train_fixed,
      n.trees = nrounds_val
    )
    
    pred_matrix_TRAINING_gbm_fixed[, stab_i] <- pred_vals_train
    predicted_vs_observed_TRAINING_gbm_fixed[[fixed_subset_index]][[stab_i]] <- data.frame(
      Observed = fixed_training_sub$PHENOTYPE,
      Predicted = pred_vals_train
    )
  
    
  
    # Shapiro-Wilk test for normality 
    shapiro_results_gbm_fixed[[fixed_subset_index]][[stab_i]] <- if (length(unique(pred_vals_fixed)) > 1) {
      shapiro_result <- shapiro.test(pred_vals_fixed)
    } else {
      warning("Predicted values are constant. Skipping Shapiro-Wilk test.")
      shapiro_result <- NA
    }
    
    # Spearman and Pearson correlations
    spearman_cor_results_gbm_fixed[[fixed_subset_index]][[stab_i]] <- cor.test(fixed_testing_sub$PHENOTYPE, pred_vals_fixed, method = "spearman", exact = FALSE)
    pearson_cor_results_gbm_fixed[[fixed_subset_index]][[stab_i]] <- cor.test(fixed_testing_sub$PHENOTYPE, pred_vals_fixed, method = "pearson")
    
    # Store accuracies
    spearman_accuracies_gbm_fixed[[fixed_subset_index]][[stab_i]] <- spearman_cor_results_gbm_fixed[[fixed_subset_index]][[stab_i]]$estimate
    pearson_accuracies_gbm_fixed[[fixed_subset_index]][[stab_i]] <- pearson_cor_results_gbm_fixed[[fixed_subset_index]][[stab_i]]$estimate
    
    # Calculate R^2
    ss_res <- sum((fixed_testing_sub$PHENOTYPE - pred_vals_fixed)^2)
    ss_tot <- sum((fixed_testing_sub$PHENOTYPE - mean(fixed_testing_sub$PHENOTYPE))^2)
    r_squared_values_gbm_fixed[[fixed_subset_index]][[stab_i]] <- 1 - (ss_res / ss_tot)
    
    # MSE
    mse_values_gbm_fixed[[fixed_subset_index]][[stab_i]] <- mean((fixed_testing_sub$PHENOTYPE - pred_vals_fixed)^2)
    
    # NDCG Calculation
    sorted_predicted_df_gbm_fixed <- predicted_vs_observed_gbm_fixed[[fixed_subset_index]][[stab_i]] %>%
      arrange(desc(Predicted))
    observed_ranking_gbm_fixed <- sorted_predicted_df_gbm_fixed$Observed
    
    log2_factor_dcg_gbm_fixed <- log2(seq_along(observed_ranking_gbm_fixed) + 1)
    ideal_relevance_gbm_fixed <- sort(observed_ranking_gbm_fixed, decreasing = TRUE)
    log2_factor_idcg_gbm_fixed <- log2(seq_along(ideal_relevance_gbm_fixed) + 1)
    
    max_length_gbm_fixed <- min(150, length(observed_ranking_gbm_fixed))
    ndcg_vals_gbm_fixed  <- numeric(max_length_gbm_fixed)
    
    for (j in 1:max_length_gbm_fixed) {
      dcg_j <- sum(observed_ranking_gbm_fixed[1:j] / log2_factor_dcg_gbm_fixed[1:j])
      idcg_j <- sum(ideal_relevance_gbm_fixed[1:j] / log2_factor_idcg_gbm_fixed[1:j])
      ndcg_vals_gbm_fixed[j] <- dcg_j / idcg_j
    }
    
    ndcg_results_gbm_fixed[[fixed_subset_index]][[stab_i]] <- data.frame(Cutoff = 1:max_length_gbm_fixed, NDCG = ndcg_vals_gbm_fixed)
    
    # Candidate Classification & Metric Calculations (threshold-based, top 25%)
    observed_threshold <- quantile(predicted_vs_observed_gbm_fixed[[fixed_subset_index]][[stab_i]]$Observed, 0.75, na.rm = TRUE)
    
    # True Class
    predicted_vs_observed_gbm_fixed[[fixed_subset_index]][[stab_i]]$True_Class <- ifelse(
      predicted_vs_observed_gbm_fixed[[fixed_subset_index]][[stab_i]]$Observed >= observed_threshold,
      "Candidate",
      "Non-Candidate"
    )
    
    # Predicted Class
    predicted_vs_observed_gbm_fixed[[fixed_subset_index]][[stab_i]]$Pred_Class <- ifelse(
      predicted_vs_observed_gbm_fixed[[fixed_subset_index]][[stab_i]]$Predicted >= observed_threshold,
      "Candidate",
      "Non-Candidate"
    )
    
    # Confusion matrix elements
    TP <- sum(
      predicted_vs_observed_gbm_fixed[[fixed_subset_index]][[stab_i]]$True_Class == "Candidate" &
        predicted_vs_observed_gbm_fixed[[fixed_subset_index]][[stab_i]]$Pred_Class == "Candidate"
    )
    TN <- sum(
      predicted_vs_observed_gbm_fixed[[fixed_subset_index]][[stab_i]]$True_Class == "Non-Candidate" &
        predicted_vs_observed_gbm_fixed[[fixed_subset_index]][[stab_i]]$Pred_Class == "Non-Candidate"
    )
    FP <- sum(
      predicted_vs_observed_gbm_fixed[[fixed_subset_index]][[stab_i]]$True_Class == "Non-Candidate" &
        predicted_vs_observed_gbm_fixed[[fixed_subset_index]][[stab_i]]$Pred_Class == "Candidate"
    )
    FN <- sum(
      predicted_vs_observed_gbm_fixed[[fixed_subset_index]][[stab_i]]$True_Class == "Candidate" &
        predicted_vs_observed_gbm_fixed[[fixed_subset_index]][[stab_i]]$Pred_Class == "Non-Candidate"
    )
    
    # Sensitivity (Recall), Specificity, Precision
    # UNDERSTANDING: of all the true “Candidate” lines, what fraction did we successfully identify? TP / (TP + FN)
    sensitivity_gbm_fixed[[fixed_subset_index]][[stab_i]] <- if ((TP + FN) > 0) TP / (TP + FN) else NA 
    # UNDERSTANDING: of all the true “Non-Candidate” lines, what fraction did we successfully identify? TN / (TN + FP)
    specificity_gbm_fixed[[fixed_subset_index]][[stab_i]] <- if ((TN + FP) > 0) TN / (TN + FP) else NA
    # UNDERSTANDING: of all the lines we identified as “Candidate”, what fraction were actually “Candidate”? TP / (TP + FP)
    precision_gbm_fixed[[fixed_subset_index]][[stab_i]] <- if ((TP + FP) > 0) TP / (TP + FP) else NA
    
    # F1 Score: formula as: F1= 2× (Precision+Recall) /  (Precision×Recall) 
    
    if (!is.na(precision_gbm_fixed[[fixed_subset_index]][[stab_i]]) && !is.na(sensitivity_gbm_fixed[[fixed_subset_index]][[stab_i]]) &&
        (precision_gbm_fixed[[fixed_subset_index]][[stab_i]] + sensitivity_gbm_fixed[[fixed_subset_index]][[stab_i]]) > 0) {
      f1_score_gbm_fixed[[fixed_subset_index]][[stab_i]] <- 2 * (precision_gbm_fixed[[fixed_subset_index]][[stab_i]] * sensitivity_gbm_fixed[[fixed_subset_index]][[stab_i]]) /
        (precision_gbm_fixed[[fixed_subset_index]][[stab_i]] + sensitivity_gbm_fixed[[fixed_subset_index]][[stab_i]])
    } else {
      f1_score_gbm_fixed[[fixed_subset_index]][[stab_i]] <- NA
    }
    
    # AUC using pROC
    true_class_factor <- factor(
      predicted_vs_observed_gbm_fixed[[fixed_subset_index]][[stab_i]]$True_Class,
      levels = c("Non-Candidate", "Candidate")
    )
    
    roc_obj <- tryCatch({  # (trycatch) I had some errors, this makes sure it doesnt error and stop.. 
      roc(true_class_factor,
          predicted_vs_observed_gbm_fixed[[fixed_subset_index]][[stab_i]]$Predicted,
          levels = c("Non-Candidate", "Candidate"),
          direction = ">")
    }, error = function(e) { NULL })
    
    auc_gbm_fixed[[fixed_subset_index]][[stab_i]] <- if (!is.null(roc_obj)) auc(roc_obj) else NA
                                                                          

        } # for (stab_i in 1:iterations_fixed)

  
  # Prediction Stability ICC for this subset
  pred_df_gbm_fixed <- as.data.frame(pred_matrix_gbm_fixed)
  icc_result_pred <- icc(pred_df_gbm_fixed, model = "oneway", type = "consistency", unit = "single")
  icc_pred_values_gbm_fixed[fixed_subset_index] <- icc_result_pred$value
  cat("Subset:", fixed_subset_index, "- Target Variable Prediction ICC =", icc_result_pred$value, "\n")
  
  # Fleiss Kappa for ENTRY Prediction SELECTION Stability
  phenotype_selection_matrix <- assign_phenotype_selection(pred_df_gbm_fixed, 
                                                           selected_proportion = PHENOTYPE_SELECTION_THRESHOLD)
  phenotype_selection_numeric <- ifelse(phenotype_selection_matrix == "selected", 1, 0)
  
  fleiss_result_pred <- kappam.fleiss(phenotype_selection_numeric)
  fleiss_pred_values_fixed[fixed_subset_index] <- fleiss_result_pred$value
  fleiss_results_df_prediction_stability <- data.frame(
    Subset = 1:num_subsets,
    Fleiss_Kappa = fleiss_pred_values_fixed,
    stringsAsFactors = FALSE
  )
  
                                                                                                
  } # for (fixed_subset_index in 1:num_subsets) LOOP 
  

# Summary of Prediction ICC across subsets
cat("\nSummary of Target Variable Prediction ICC across subsets:\n")
print(icc_pred_values_gbm_fixed)

cat("\nSummary of Target Variable Prediction Fleiss' Kappa across subsets:\n")
print(fleiss_results_df_prediction_stability)

### End measuring time ###
end_time <- Sys.time()

# Calculate and print the elapsed time
elapsed_time <- end_time - start_time
cat("Total compute time: ", elapsed_time, "\n")

# Unregister parallel backend once done
stopImplicitCluster()


# Exporting Regression Metrics (Correlation Coefficients, R squared, MSE)

# Exporting Correlation Metrics (Spearman and Pearson)
spearman_accuracies_df <- data.frame(
  Subset = rep(1:num_subsets, each = iterations_fixed),
  Iteration = rep(1:iterations_fixed, num_subsets),
  Spearman = unlist(spearman_accuracies_gbm_fixed)
)
write.csv(spearman_accuracies_df, file.path(results_folder, "spearman_accuracies.csv"), row.names = FALSE)

pearson_accuracies_df <- data.frame(
  Subset = rep(1:num_subsets, each = iterations_fixed),
  Iteration = rep(1:iterations_fixed, num_subsets),
  Pearson = unlist(pearson_accuracies_gbm_fixed)
)
write.csv(pearson_accuracies_df, file.path(results_folder, "pearson_accuracies.csv"), row.names = FALSE)

r_squared_df <- data.frame(
  Subset = rep(1:num_subsets, each = iterations_fixed),
  Iteration = rep(1:iterations_fixed, num_subsets),
  R_squared = unlist(r_squared_values_gbm_fixed)
)
write.csv(r_squared_df, file.path(results_folder, "r_squared_values.csv"), row.names = FALSE)

mse_df <- data.frame(
  Subset = rep(1:num_subsets, each = iterations_fixed),
  Iteration = rep(1:iterations_fixed, num_subsets),
  MSE = unlist(mse_values_gbm_fixed)
)
write.csv(mse_df, file.path(results_folder, "mse_values.csv"), row.names = FALSE)


# Exporting ranking metrics (NDCG)
# Initialize data frames to store NDCG values across subsets and iterations for each threshold
ndcg_at_1_df <- data.frame()
ndcg_at_5_df <- data.frame()
ndcg_at_10_df <- data.frame()
ndcg_at_20_percent_df <- data.frame()
ndcg_at_10_percent_df <- data.frame()
ndcg_at_5_percent_df <- data.frame()
ndcg_range_1_5_df <- data.frame()
ndcg_range_1_10_df <- data.frame()
ndcg_range_1_20_df <- data.frame()

# Define cutoff thresholds as percentages of testing size (for our context)
cutoff_20_percent <- round(nrow(fixed_testing_set) * 0.20)
cutoff_10_percent <- round(nrow(fixed_testing_set) * 0.10)
cutoff_5_percent <- round(nrow(fixed_testing_set) * 0.05)

# Loop through subsets
for (fixed_subset_index in 1:num_subsets) {
  for (stab_i in 1:iterations_fixed) {
    
    # Extract NDCG values for each threshold and range for the current subset and iteration
    ndcg_at_1_val <- ndcg_results_gbm_fixed[[fixed_subset_index]][[stab_i]]$NDCG[ndcg_results_gbm_fixed[[fixed_subset_index]][[stab_i]]$Cutoff == 1]
    ndcg_at_5_val <- ndcg_results_gbm_fixed[[fixed_subset_index]][[stab_i]]$NDCG[ndcg_results_gbm_fixed[[fixed_subset_index]][[stab_i]]$Cutoff == 5]
    ndcg_at_10_val <- ndcg_results_gbm_fixed[[fixed_subset_index]][[stab_i]]$NDCG[ndcg_results_gbm_fixed[[fixed_subset_index]][[stab_i]]$Cutoff == 10]
    ndcg_at_20_percent_val <- ndcg_results_gbm_fixed[[fixed_subset_index]][[stab_i]]$NDCG[ndcg_results_gbm_fixed[[fixed_subset_index]][[stab_i]]$Cutoff == cutoff_20_percent]
    ndcg_at_10_percent_val <- ndcg_results_gbm_fixed[[fixed_subset_index]][[stab_i]]$NDCG[ndcg_results_gbm_fixed[[fixed_subset_index]][[stab_i]]$Cutoff == cutoff_10_percent]
    ndcg_at_5_percent_val <- ndcg_results_gbm_fixed[[fixed_subset_index]][[stab_i]]$NDCG[ndcg_results_gbm_fixed[[fixed_subset_index]][[stab_i]]$Cutoff == cutoff_5_percent]
    
    # Calculate the mean NDCG for each range
    mean_ndcg_range_1_5_val <- mean(head(ndcg_results_gbm_fixed[[fixed_subset_index]][[stab_i]]$NDCG, 5))
    mean_ndcg_range_1_10_val <- mean(head(ndcg_results_gbm_fixed[[fixed_subset_index]][[stab_i]]$NDCG, 10))
    mean_ndcg_range_1_20_val <- mean(head(ndcg_results_gbm_fixed[[fixed_subset_index]][[stab_i]]$NDCG, cutoff_20_percent))
    
    # Add the NDCG values for each iteration to the dataframes
    ndcg_at_1_df <- rbind(ndcg_at_1_df, data.frame(Subset = fixed_subset_index, Iteration = stab_i, NDCG = ndcg_at_1_val))
    ndcg_at_5_df <- rbind(ndcg_at_5_df, data.frame(Subset = fixed_subset_index, Iteration = stab_i, NDCG = ndcg_at_5_val))
    ndcg_at_10_df <- rbind(ndcg_at_10_df, data.frame(Subset = fixed_subset_index, Iteration = stab_i, NDCG = ndcg_at_10_val))
    ndcg_at_20_percent_df <- rbind(ndcg_at_20_percent_df, data.frame(Subset = fixed_subset_index, Iteration = stab_i, NDCG = ndcg_at_20_percent_val))
    ndcg_at_10_percent_df <- rbind(ndcg_at_10_percent_df, data.frame(Subset = fixed_subset_index, Iteration = stab_i, NDCG = ndcg_at_10_percent_val))
    ndcg_at_5_percent_df <- rbind(ndcg_at_5_percent_df, data.frame(Subset = fixed_subset_index, Iteration = stab_i, NDCG = ndcg_at_5_percent_val))
    
    ndcg_range_1_5_df <- rbind(ndcg_range_1_5_df, data.frame(Subset = fixed_subset_index, Iteration = stab_i, NDCG = mean_ndcg_range_1_5_val))
    ndcg_range_1_10_df <- rbind(ndcg_range_1_10_df, data.frame(Subset = fixed_subset_index, Iteration = stab_i, NDCG = mean_ndcg_range_1_10_val))
    ndcg_range_1_20_df <- rbind(ndcg_range_1_20_df, data.frame(Subset = fixed_subset_index, Iteration = stab_i, NDCG = mean_ndcg_range_1_20_val))
  }
}

# Export each dataframe for the different thresholds and ranges, now including subsets and iterations
write.csv(ndcg_at_1_df, file.path(results_folder, "ndcg_at_1.csv"), row.names = FALSE)
write.csv(ndcg_at_5_df, file.path(results_folder, "ndcg_at_5.csv"), row.names = FALSE)
write.csv(ndcg_at_10_df, file.path(results_folder, "ndcg_at_10.csv"), row.names = FALSE)
write.csv(ndcg_at_20_percent_df, file.path(results_folder, "ndcg_at_20_percent.csv"), row.names = FALSE)
write.csv(ndcg_at_10_percent_df, file.path(results_folder, "ndcg_at_10_percent.csv"), row.names = FALSE)
write.csv(ndcg_at_5_percent_df, file.path(results_folder, "ndcg_at_5_percent.csv"), row.names = FALSE)
write.csv(ndcg_range_1_5_df, file.path(results_folder, "ndcg_range_1_5.csv"), row.names = FALSE)
write.csv(ndcg_range_1_10_df, file.path(results_folder, "ndcg_range_1_10.csv"), row.names = FALSE)
write.csv(ndcg_range_1_20_df, file.path(results_folder, "ndcg_range_1_20.csv"), row.names = FALSE)

cat("All NDCG dataframes for subsets and iterations have been successfully exported.\n")


# Exporting Classification Metrics (Sensitivity, Specificity, Precision, F1, AUC)
classification_df <- data.frame(
  Subset = rep(1:num_subsets, each = iterations_fixed),
  Iteration = rep(1:iterations_fixed, num_subsets),
  Sensitivity = unlist(sensitivity_gbm_fixed),
  Specificity = unlist(specificity_gbm_fixed),
  Precision = unlist(precision_gbm_fixed),
  F1_Score = unlist(f1_score_gbm_fixed),
  AUC = unlist(auc_gbm_fixed)
)
write.csv(classification_df, file.path(results_folder, "classification_metrics.csv"), row.names = FALSE)

# Exporting Stability Metrics (ICC and Fleiss Kappa)
icc_df <- data.frame(
  Subset = 1:num_subsets,
 # ICC_Feature_Importance = unlist(icc_results_list_feature_matrix_1),
  ICC_Prediction_Stability = unlist(icc_pred_values_gbm_fixed)
)
write.csv(icc_df, file.path(results_folder, "icc_results.csv"), row.names = FALSE)

fleiss_df <- data.frame(
  Subset = 1:num_subsets,
 # Fleiss_Kappa_Feature_Importance = fleiss_feature_importance_results_df$Fleiss_Kappa,
  Fleiss_Kappa_Prediction_Stability = fleiss_results_df_prediction_stability$Fleiss_Kappa
)
write.csv(fleiss_df, file.path(results_folder, "fleiss_kappa_results.csv"), row.names = FALSE)

# Confirming the export of files
cat("Exported regression, classification, and stability metrics to the 'Results' folder.\n")


# Save the Workspace
workspace_file <- file.path(results_folder, paste0("Workspace_trait_", phenotype_of_interest_config_filename, "_", "_GBM_package.RData"))

# if workspace already exists, concatenatew timestamp 
if (file.exists(workspace_file)) {
  timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
  workspace_file <- file.path(
    results_folder,
    paste0("Workspace_trait_", phenotype_of_interest_config_filename, "_", "_GBM_package_", timestamp, ".RData")
  )
}

# Save the workspace from the current (local) environment
save(list = ls(all.names = TRUE), file = workspace_file, envir = environment())
cat("Workspace saved as", basename(workspace_file),
    "for trait", phenotype_of_interest,
    "in the 'GBM_package_final_model_training_Environment_01.02.2025' folder.\n")



#  Super Dataframe 

# Add an Iteration column with value 0 to dataframes that don't have it
icc_df_super <- icc_df %>% mutate(Iteration = 0)
fleiss_df_super <- fleiss_df %>% mutate(Iteration = 0)

# Function to convert a wide dataframe to long format
convert_to_long <- function(df) {
  pivot_longer(df, cols = -c(Subset, Iteration), names_to = "Metric", values_to = "Value")
}

# List all dataframes that need to be reshaped
df_list <- list(
  spearman_accuracies_df = spearman_accuracies_df,
  pearson_accuracies_df = pearson_accuracies_df,
  r_squared_df = r_squared_df,
  mse_df = mse_df,
  # ndcg_at_1_df = ndcg_at_1_df,
  # ndcg_at_5_df = ndcg_at_5_df,
  # ndcg_at_10_df = ndcg_at_10_df,
  ndcg_at_20_percent_df = ndcg_at_20_percent_df,
  # ndcg_at_10_percent_df = ndcg_at_10_percent_df,
  # ndcg_at_5_percent_df = ndcg_at_5_percent_df,
  # ndcg_range_1_5_df = ndcg_range_1_5_df,
  # ndcg_range_1_10_df = ndcg_range_1_10_df,
  # ndcg_range_1_20_df = ndcg_range_1_20_df,
  classification_df = classification_df,
  icc_df_super = icc_df_super,
  fleiss_df_super = fleiss_df_super
)

# Apply the conversion to each dataframe
long_df_list <- lapply(df_list, convert_to_long)

# Merge all long-format dataframes into one master dataframe.
super_df <- bind_rows(long_df_list, .id = "Source")

# Add constant parameter columns: 
super_df <- super_df %>%
  mutate(
    ParameterLabel = param_label,
    LearningRate   = learning_rate,
    BoostRounds    = boost_rounds, 
    Phenotype      = phenotype_of_interest
  ) %>%
  # Rearrange columns in the desired order.
  select(ParameterLabel,Phenotype, Subset, Iteration, LearningRate, BoostRounds, Metric, Value)

# Save the super dataframe to a CSV file
write.csv(super_df, file.path( results_folder, "super_dataframe.csv"), row.names = FALSE)


#### ADDITION TO ENSURE I DONT GET LOSTN IN WHAT WAS RUN OR WHAT OUTPUT IS WHAT ####

# --- Flexible Script Path Detection ---
get_script_path <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    current_path <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(current_path)) return(current_path)
  }
  # Fallback for Rscript
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  fileArg <- "--file="
  scriptPath <- sub(fileArg, "", cmdArgs[grep(fileArg, cmdArgs)])
  if (length(scriptPath) == 0) scriptPath <- ""
  return(scriptPath)
}

# Automatically determine the current script's path
script_path <- get_script_path()

# --- Save a copy of the current script to the results folder ---

script_copy_path <- file.path(results_folder, paste0("script_copy_",phenotype_of_interest_config_filename, "_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".R"))
file.copy(from = script_path, to = script_copy_path)
cat("Script copied to", script_copy_path, "\n")

# --- Save session information for reproducibility ---
session_info_file <- file.path(results_folder, paste0("session_info_",phenotype_of_interest_config_filename, "_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".txt"))
capture.output(sessionInfo(), file = session_info_file)
cat("Session information saved to", session_info_file, "\n")

# --- Save key configuration values to a text file ---
config_file <- file.path(results_folder, paste0("config_TRAIT_", phenotype_of_interest_config_filename, "_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".txt"))
config_info <- c(
  paste("gbm_model_fixed call:", model_call_text),
  paste("tune_grid_fixed:","\n", paste(capture.output(print(tune_grid_arg)), collapse = "\n")),
  paste("all_traits:", paste(all_traits, collapse = ", ")),
  paste("phenotype_of_interest (TRAINING):", phenotype_of_interest),
  paste("phenotype_of_interest_year2 (TESTING):", phenotype_of_interest_year2),
  paste("iterations_fixed:", iterations_fixed),
  paste("num_subsets:", num_subsets),
  paste("FEATURE_SELECTION_THRESHOLD (proportion):", FEATURE_SELECTION_THRESHOLD),
  paste("PHENOTYPE_SELECTION_THRESHOLD (proportion):", PHENOTYPE_SELECTION_THRESHOLD),
  paste("Elapsed_Time:", elapsed_time)
)
writeLines(config_info, con = config_file)
cat("Configuration settings saved to", config_file, "\n")
}

