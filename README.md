# Multi-objective-Evaluation-and-Optimization-of-Stochastic-Gradient-Boosting-Machine)  - Publishing

To run this analysis, first download the dataset from Scott et al. (2021) here:  http://mtweb.cs.ucl.ac.uk/mus/www/MAGICdiverse/
Scott, M. F., Fradgley, N., Bentley, A. R., Brabbs, T., Corke, F., Gardner, K. A., … Cockram, J. (2021). Genome Biology, 22(1). doi:10.1186/s13059-021-02354-7

### Workflow

#### Step 1: Data Preparation and Partitioning

Run the `Dataset preparation and Imputation.R` script. This script will:
1.  Process the raw phenotypic and genomic datasets.
2.  Impute missing values using the `missRanger` package.
3.  Generate five sets of 80/20 training and testing data splits, exporting them as CSV files.

#### Step 2: Model Execution

Before running the models, create the following directory structure to store the output files:
*   Coding/ModellingScripts/COMPLETE_imputed_data_run/GBM/GY/
*   Coding/ModellingScripts/COMPLETE_imputed_data_run/GBM/GPC/
*   Coding/ModellingScripts/COMPLETE_imputed_data_run/GBM/TGW/
*   Coding/ModellingScripts/COMPLETE_imputed_data_run/GBM/HET/

To begin the analysis, run the main initiator script: `GBM_Initiator_Script_HM_2025.R`.

**Important:** The worker script, `GBM_Worker_Script_HM_2025.R`, must be located in the same directory as the initiator script for the models to run correctly.

#### Step 3: Generating Figures

To generate the publication figures, first ensure the summary data is prepared. The plotting scripts require a merged `super_dataframe.csv` file for each model combination.

*   **Option A (Recommended):** Use the pre-merged CSV files provided in the `Coding/EvaluationPlotterScript_and_Evaluation_Data/` directory.
*   **Option B:** Manually merge the individual CSV files that were exported to the trait-specific folders during Step 2.

Once the data is ready, run the plotting script: `Super_Plotter_PUBLISHING_V5_CB_25.06.R`.

#### Step 4: Statistical Analysis (ANOVA)

To perform the ANOVA, pairwise comparisons, and `emmeans` analysis, run the script located in the `Coding/ANOVA_PAIRWISE_Script_and_Evaluation_Data/` directory. The necessary dataframes are also provided in this location.
