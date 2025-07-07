
# Dataset preparation and Imputation 
# See Scott et al. 2021 for datasets
# Scott, M. F., Fradgley, N., Bentley, A. R., Brabbs, T., Corke, F., Gardner, K. A., … Cockram, J. (2021). Genome Biology, 22(1). doi:10.1186/s13059-021-02354-7
# Dataset can be found via the study above and : http://mtweb.cs.ucl.ac.uk/mus/www/MAGICdiverse/


# First, read in the phenotypic data and reduce to the trait you want to analyse
Y = read.table("PhenotypicData.csv", sep=";", header=T)

# read in the genomic data and merge with the phenotypic data
X = read.table("GenomicData.csv", header=T, sep=";", check.names=F)


#Verify that 'line_name' column exists in both datasets
if (!("line_name" %in% colnames(X)) | !("line_name" %in% colnames(Y))) {
  stop("The 'line_name' column must exist in both datasets")
}

# Identify non-matches
X_not_in_Y <- setdiff(X$line_name, Y$line_name)
Y_not_in_X <- setdiff(Y$line_name, X$line_name)

# Count and report non-matches
cat("Number of 'line_name' entries in X not in Y:", length(X_not_in_Y), "\n")
cat("Number of 'line_name' entries in Y not in X:", length(Y_not_in_X), "\n")

# Print non-matching entries for review
print(X_not_in_Y)
print(Y_not_in_X)

# Merge the datasets based on 'line_name'
# Merge the datasets based on 'line_name'
raw_data <- merge(Y, X, by="line_name")
nrow(raw_data)  # number of rows in the merged dataset = 500
ncol(Y)
ncol(X)
ncol(raw_data)


# replace hyphens with underscores or another non-operator character.
names(raw_data) <- gsub("-", "_", names(raw_data))

# Create the raw_data_phenotype subset from the first column to the 72nd
raw_data_phenotype <- raw_data[, 1:72]
write.csv(raw_data_phenotype, "RawPhenotypicData.csv", row.names = FALSE)

# Create the raw_data_genotype subset from the 73rd column to the last
SNP_data <- raw_data[, 73:ncol(raw_data)]
SNP_data[1:5,1:5]
# Add the 'line_name' column to the data_genotype subset
data_genotype <- data.frame(line_name = raw_data$line_name, SNP_data)


# Clean the SNP DATA ####
# Ensure missing values are properly encoded
data_genotype[data_genotype == "na"] <- NA
summary(data_genotype[1:5])



# Identify columns with only one unique value excluding NAs (and excluding the first column)
ncol(data_genotype)
columns_with_ONE_unique_values <- sapply(data_genotype[, -1], function(x) length(unique(x[!is.na(x)])) == 1)

# how many columns have the same value
sum(columns_with_ONE_unique_values)
data_genotype$AX_94399399
summary(data_genotype$AX_94399399) 

# Exclude these columns by keeping columns that do not meet the above condition
data_genotype_filtered <- data_genotype[, c(TRUE, !columns_with_ONE_unique_values)]

# Print the number of columns in the filtered data frame
print(ncol(data_genotype_filtered))

# Updating the original data_genotype variable with the filtered data frame
data_genotype <- data_genotype_filtered

df <-  data_genotype
# Count and sort NA values in each column
na_count_sorted <- sapply(df, function(x) sum(is.na(x))) %>% 
  sort(decreasing = TRUE)

# Print column number, name, and NA counts, sorted by counts
cat("Column Name, Counts of NA:\n")
names(na_count_sorted) <- paste0(names(na_count_sorted))
print(na_count_sorted)

# Calculate the threshold for NA percentage
na_threshold <- nrow(df) * 0.05

# Count NA values in each column
na_counts <- sapply(df, function(x) sum(is.na(x)))

# Identify columns with more than 5% NA values
columns_to_remove <- names(na_counts[na_counts > na_threshold])

# Print the names of columns to be removed
if (length(columns_to_remove) > 0) {
  cat("Total SNP loci above NA threshold:", paste(length(columns_to_remove)), "\n", "Columns removed due to more than 5% NA values:\n", paste(columns_to_remove, collapse = ", "), "\n")
} else {
  cat("No columns removed based on the 5% NA values criterion.\n")
}
ncol(df)  # number of columns in the dataframe before removal

# Remove the identified columns from the dataframe
df <- df[, !(names(df) %in% columns_to_remove)]
ncol(df) 

df <- write.csv(df , "FilteredSNPData_col_missNA_need_impute.csv", header = TRUE, stringsAsFactors = FALSE)

#### missRanger Imputation for SNP Data ####

# overview:
# script used to impute missing values in SNP data using missRanger
# SNP data  treated as categorical by converting to factors for the Random Forest algorithm
# mtry for the number of SNPs/variables/features randomly sampled at each split, using the default value: floor(sqrt(number of features)).



### Function to check and install packages  ####
install_if_not_present <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
  }
}

# packages to be installed and loaded
packages <- c( "dplyr",  "kernlab",  
               "foreach", "caret", "missRanger")

# install packages 
sapply(packages, install_if_not_present)

# Load
sapply(packages, library, character.only = TRUE)


# read in 
df <- read.csv("FilteredSNPData_col_missNA_need_impute.csv", header = TRUE, stringsAsFactors = FALSE)

# Exclude 'line_name' column before imputation
snp_data_only <- df[-1]  #remove 'line_name' in the first column
snp_data_only[1:5, 1:5]

#  snp_data_only 
snp_data_only_factors <- snp_data_only

# Convert all columns to factors, treated as categorical by randomForest
for(i in seq_along(snp_data_only_factors)) {
  snp_data_only_factors[[i]] <- as.factor(snp_data_only_factors[[i]])
}


# https://cran.r-project.org/web/packages/missRanger/vignettes/missRanger.html

# missRanger ####
# Using missRanger for imputation
library(missRanger)

# Perform Random Forest imputation on SNP data only using missRanger
system.time({
  imputing_snp_data_ranger <- missRanger(snp_data_only_factors, 
                                         num.trees = 100,  
                                         mtry = floor(sqrt(ncol(snp_data_only_factors))), # default mtry, manually set, 
                                         verbose = TRUE,
                                         num.threads = 70, 
                                         data_only = FALSE,  # more info can look into it if needed after the imputation
                                         pmm.k = 5) # pmm.k (Predictive Mean Matching with k neighbors), helps with reliability, the number of nearest neighbors to use for the PMM imputation method, finding the k closest predicted probabilities 
})


# Save entire output 
saveRDS(imputing_snp_data_ranger, "FINAL_imputing_snp_data_ranger_output.rds")

# missRanger directly returns the imputed data
imputed_snp_data_Ranger <- imputing_snp_data_ranger$data

# bind  'line_name' column to the imputed SNP data
df_imputed_reattached_line_name <- cbind(df[1], imputed_snp_data_Ranger)

# save
write.csv(df_imputed_reattached_line_name, "Final_Lab_imputed_snp_data_missRanger_IMPORTANT_nt100_mtrSqrt.csv")

# SNP imputed, ready for joining back with traits, and the next steps: create subsets, run predictions 




### outlier cleanup ####

# First, read in the phenotypic data and reduce to the trait you want to analyse
Y = read.table("PhenotypicData.csv", sep=";", header=T)

Y[1:10,1:15]

# Selecting traits of interest for this thesis
Y_ThesisTraits <-  data.frame(Y$line_name, Y$GY_1, Y$GY_2, Y$TGW_1, Y$TGW_2, Y$YR_1,Y$YR_2, Y$GPC_1, Y$GPC_2, Y$HET_1, Y$HET_2)

#Makes columns named as original
names(Y_ThesisTraits) <- c("line_name", "GY_1", "GY_2", "TGW_1", "TGW_2", "YR_1", "YR_2", "GPC_1", "GPC_2", "HET_1", "HET_2")


# Calculate the number of NA values per column
na_count <- sapply(Y_ThesisTraits, function(x) sum(is.na(x)))

# Sort columns by the number of NA values in descending order
sorted_na_count <- sort(na_count, decreasing = TRUE)

# Display the sorted list of columns by their number of NA values
sorted_na_count


# Remove rows with any missing values (NA or empty)
Y_clean <- na.omit(Y_ThesisTraits)

nrow(Y)
nrow(Y_clean)

Y_clean[1:10,1:10]

# For each trait, 2nd column onwards, remove all entries with outliers for any of the traits. 
# Outlier removal using IQR method for each trait column
for (i in 2:ncol(Y_clean)) {  # start loop from 2 to skip 'line_name'
  Q1 <- quantile(Y_clean[[i]], 0.25, na.rm = TRUE)
  Q3 <- quantile(Y_clean[[i]], 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 3 * IQR
  upper_bound <- Q3 + 3 * IQR
  
  # Subset data frame to exclude outliers
  Y_clean <- Y_clean[Y_clean[[i]] >= lower_bound & Y_clean[[i]] <= upper_bound, ]
}

# Number of rows after cleaning (entries with na and outliers removed)
print(nrow(Y_clean))


# see how the line_name matches are between pheotype and genotype datasets

# read in the genomic data and merge with the phenotypic data
X = read.table("Final_Lab_imputed_snp_data_missRanger_IMPORTANT_nt100_mtrSqrt.csv", sep = ",", header=T,  check.names=F)


XRawImpute <-  X

dim(XRawImpute)

XRawImpute[1:10,1:10]
# Remove the first column (1:no. col, redundant)
X <-  XRawImpute[,-1]
X[1:5,1:5]

# Step 1: Verify that 'line_name' column exists in both datasets
if (!("line_name" %in% colnames(X)) | !("line_name" %in% colnames(Y_clean))) {
  stop("The 'line_name' column must exist in both datasets")
}



#  Identify non-matches
X_not_in_Y <- setdiff(X$line_name, Y_clean$line_name)
Y_not_in_X <- setdiff(Y_clean$line_name, X$line_name)

#  Count and report non-matches
cat("Number of 'line_name' entries in X not in Y_clean:", length(X_not_in_Y), "\n")
cat("Number of 'line_name' entries in Y_clean not in X:", length(Y_not_in_X), "\n")

# Optional: Print non-matching entries for review
print(X_not_in_Y)
print(Y_not_in_X)

# Step 5: Merge the datasets based on 'line_name'
# Merge the datasets based on 'line_name'
Master_dataset_clean_imputed <- merge(Y_clean, X, by="line_name")

nrow(Master_dataset_clean_imputed)
Master_dataset_clean_imputed[1,10:20]

# save MASTER as csv file 
write.csv(Master_dataset_clean_imputed, "Master_dataset_clean_imputed_13.05.2024.csv", row.names = FALSE)


# SUBSETTING: WORKING DONT RUN UNNECESSARILY ####
iterations <- 50

set.seed(123)

# Load the dataset
wheat_dat <- read.csv("Master_dataset_clean_imputed_13.05.2024.csv", header = TRUE)
# Renaming the column GY_1 to PHENOTYPE
# names(wheat_dat)[names(wheat_dat) == "GY_1"] <- "PHENOTYPE"

# determine sizes for training and testing
total_entries <- nrow(wheat_dat)
training_size <- round(total_entries * 0.80)

# Generate indices
training_indices <- vector("list", iterations)
testing_indices <- vector("list", iterations)

for (i in 1:iterations) {
  inds <- sample(total_entries)
  training_indices[[i]] <- inds[1:training_size]
  testing_indices[[i]] <- inds[(training_size + 1):total_entries]
}

# Check if the folder exists, if not create it
folder_name <- "AllTraitsMasterLoopTrainingTestingMissRangerImputed"
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
}

# Save indices to CSV files
lapply(1:iterations, function(i) {
  write.csv(wheat_dat[training_indices[[i]], ], file = paste0(folder_name,"/training_set_", i, ".csv"), row.names = FALSE)
  write.csv(wheat_dat[testing_indices[[i]], ], file = paste0(folder_name,"/testing_set_", i, ".csv"), row.names = FALSE)
})

