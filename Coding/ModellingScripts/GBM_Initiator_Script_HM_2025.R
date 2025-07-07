#############################################
## GBM_Initiator_Script_HM_2025.R
#############################################
# Former r file name, "run_all_model_pieces_Across_all_traits.R"




####1 Install / load packages 
install_if_not_present <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}


# list of packages 
packages <- c( "gbm", "dplyr", "tidyr", "pROC", "doParallel","rstudioapi"  , "foreach", "irr")
  


# Install + load them in the main session
invisible(lapply(packages, install_if_not_present))
invisible(lapply(packages, function(pkg) library(pkg, character.only = TRUE)))


#### 2. source the worker 
#ORIGINAL NAME: 
#source("DEFINITIVE_MT_gbm.fit_gbm()_TASK1_GradientBoost_DIRECT_rmFeat_P_matrices_parallel_MODEL_GRID_3.0_MasterRun_1.4.R")

#RENAMED: 
source("GBM_Worker_Script_HM_2025.R")

#### 3. Build Model Grid ####
shrinkage_seq <- seq(0.0025, 0.2, length.out = 6)
ntrees_seq    <- seq(100, 5000, length.out = 6)

model_grid <- expand.grid(shrinkage = shrinkage_seq,
                          ntrees    = ntrees_seq)

# sett all parameters to default, set based on r documentation: 
model_grid$interaction.depth <-1
model_grid$n.minobsinnode  <- 10
model_grid$bag.fraction <- 0.5
model_grid$train.fraction  <- 1
model_grid$cv.folds          <- 0  # achieved via manual subsets, taking pairs (i.e testing_set_1.csv and training_set_1.csv from folder, "AllTraitsMasterLoopTrainingTestingMissRangerImputed" etc... )
model_grid$keep.data    <- TRUE
model_grid$verbose  <- TRUE
model_grid$distribution     <-"gaussian"

# Split into a list of data frames, one per row
my_36_grids<- split(model_grid, seq(nrow(model_grid)))
total_grids <- length(my_36_grids)
cat("Number of parameter combos:", total_grids, "\n")


#### 4. setup Parallel Workers with doParallel
num_workers <- 4

library(doParallel)
cl <- makeCluster(num_workers)
registerDoParallel(cl)
cat("Started cluster with", num_workers, "parallel workers.\n")


## 5. parallel Loop with foreach 
start_index <- 1   
end_index  <- total_grids  


#### TRAITS ####

### GY ####

results <- foreach(i = start_index:end_index,
                   .packages = c("gbm", "dplyr", "tidyr", "irr", "pROC", "rstudioapi", "doParallel", "foreach"),
                   .export   = c("my_36_grids","run_experiment"),  # may get warning: "In e$fun(obj, substitute(ex), parent.frame(), e$data) : already exporting variable(s): my_36_grids, run_experiment" not an issue. 
                   .errorhandling = "pass") %dopar% {

                     # Pull out the parameter row
                     this_grid   <- my_36_grids[[i]]
                     # Bookkeeping for labeling
                     model_label <- paste0("GBM_Model_", i,
                                           "_eta_", this_grid$shrinkage,
                                           "_ntrees_", this_grid$ntrees)

                     cat("\n***** Worker running combo #", i, ":", model_label, "\n")

                     
                     run_experiment(
                       tune_grid_arg          = this_grid,
                       iterations_fixed       = 10,
                       num_subsets            = 5,
                       param_label            = model_label,
                       phenotype_of_interest  = "GY_1",
                       phenotype_of_interest_year2 = "GY_2",
                       learning_rate          = this_grid$shrinkage,
                       boost_rounds           = this_grid$ntrees,
                       base_path              = "COMPLETE_imputed_data_run/GBM/GY"  # overrides the default worker script base path
                     )
                     full_save_path <- file.path(base_path, paste0("my_results_", param_label, ".csv"))
                     write.csv(my_results, full_save_path)
                     return(paste("Done with combo #", i, "on worker", Sys.getpid()))
                   }


#### GPC ####
results <- foreach(i = start_index:end_index,
                   .packages = c("gbm", "dplyr", "tidyr", "irr", "pROC", "rstudioapi", "doParallel", "foreach"),
                   .export   = c("my_36_grids","run_experiment"), # .export sends local objects to each worker 
                   .errorhandling = "pass") %dopar% {

                     # Extract the parameter row
                     this_grid   <- my_36_grids[[i]]
                     # Bookkeeping for labeling
                     model_label <- paste0("GBM_Model_", i,
                                           "_eta_", this_grid$shrinkage,
                                           "_ntrees_", this_grid$ntrees)

                     cat("\n***** Worker running combo #", i, ":", model_label, "\n")

                     
                     run_experiment(
                       tune_grid_arg          = this_grid,
                       iterations_fixed       = 10,
                       num_subsets            = 5,  # worker script set default is 10
                       param_label            = model_label,
                       phenotype_of_interest  = "GPC_1",
                       phenotype_of_interest_year2 = "GPC_2",
                       learning_rate          = this_grid$shrinkage,
                       boost_rounds           = this_grid$ntrees,
                       base_path              = "COMPLETE_imputed_data_run/GBM/GPC"
                     )
                     full_save_path <- file.path(base_path, paste0("my_results_", param_label, ".csv"))
                     write.csv(my_results, full_save_path)
                     return(paste("Done with combo #", i, "on worker", Sys.getpid()))
                   }

#### TGW ####
results <- foreach(i = start_index:end_index, 
                   .packages = c("gbm", "dplyr", "tidyr", "irr", "pROC", "rstudioapi", "doParallel", "foreach"),
                   .export   = c("my_36_grids","run_experiment"),
                   .errorhandling = "pass") %dopar% {
                     
                     # Pull out the parameter row
                     this_grid   <- my_36_grids[[i]]
                     # Bookkeeping for labeling
                     model_label <- paste0("GBM_Model_", i, 
                                           "_eta_", this_grid$shrinkage, 
                                           "_ntrees_", this_grid$ntrees)
                     
                     cat("\n***** Worker running combo #", i, ":", model_label, "\n")
                     
                     
                     run_experiment(
                       tune_grid_arg          = this_grid,
                       iterations_fixed       = 10,
                       num_subsets            = 5,
                       param_label            = model_label,
                       phenotype_of_interest  = "TGW_1",
                       phenotype_of_interest_year2 = "TGW_2",
                       learning_rate          = this_grid$shrinkage,
                       boost_rounds           = this_grid$ntrees,
                       base_path              = "COMPLETE_imputed_data_run/GBM/TGW"
                     )
                     full_save_path <- file.path(base_path, paste0("my_results_", param_label, ".csv"))
                     write.csv(my_results, full_save_path)
                     return(paste("Done with combo #", i, "on worker", Sys.getpid()))
                   }
#### HET ####
results <- foreach(i = start_index:end_index, 
                   .packages = c("gbm", "dplyr", "tidyr", "irr", "pROC", "rstudioapi", "doParallel", "foreach"),
                   .export   = c("my_36_grids","run_experiment"),
                   .errorhandling = "pass") %dopar% {
                     
                     # Pull out the parameter row
                     this_grid   <- my_36_grids[[i]]
                     # Bookkeeping for labeling
                     model_label <- paste0("GBM_Model_", i, 
                                           "_eta_", this_grid$shrinkage, 
                                           "_ntrees_", this_grid$ntrees)
                     
                     cat("\n***** Worker running combo #", i, ":", model_label, "\n")
                     
                     
                     run_experiment(
                       tune_grid_arg          = this_grid,
                       iterations_fixed       = 10,
                       num_subsets            = 5,
                       param_label            = model_label,
                       phenotype_of_interest  = "HET_1",
                       phenotype_of_interest_year2 = "HET_2",
                       learning_rate          = this_grid$shrinkage,
                       boost_rounds           = this_grid$ntrees,
                       base_path              = "COMPLETE_imputed_data_run/GBM/HET"
                     )
                     full_save_path <- file.path(base_path, paste0("my_results_", param_label, ".csv"))
                     write.csv(my_results, full_save_path)

                     return(paste("Done with combo #", i, "on worker", Sys.getpid()))
                   }


### 6 Shut Down Cluster 
stopCluster(cl)

# Done. 

