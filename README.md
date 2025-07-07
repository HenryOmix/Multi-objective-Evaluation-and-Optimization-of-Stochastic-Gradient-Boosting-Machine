# Multi-objective-Evaluation-and-Optimization-of-Stochastic-Gradient-Boosting-Machine)  - Publishing

To run this analysis, first download the dataset from Scott et al. (2021) here:  http://mtweb.cs.ucl.ac.uk/mus/www/MAGICdiverse/
Scott, M. F., Fradgley, N., Bentley, A. R., Brabbs, T., Corke, F., Gardner, K. A., … Cockram, J. (2021). Genome Biology, 22(1). doi:10.1186/s13059-021-02354-7

1) Use the "Dataset preparation and Imputation.R" script to wrangle the Phenotypic and Genomic Datasets, impute missing NA via missRanger, export 80-20 splits of the dataset into testing and training splits   (default is 50 splits, we used 5)

2) Use the Intitiator script, "GBM_Initiator_Script_HM_2025.R" to begin modelling and evaluation metric estimation.
NOTE: "GBM_Worker_Script_HM_2025.R" must be in the same location as the intiator script, Before modelling intiation,  set up folders for the outputs as shown in the repository: 
Grain yield: "Coding/ModellingScripts/COMPLETE_imputed_data_run/GBM/GY"
Grain Protein Content: "Coding/ModellingScripts/COMPLETE_imputed_data_run/GBM/GPC"
and so on.

3) Plotting, merge "super_dataframe.csv" for each model combination exported in their respective folders for  each trait. (merged csv files provided: "Coding/EvaluationPlotterScript_and_Evaluation_Data"
Plotter: run the plotting scripts to produce plots: "Super_Plotter_PUBLISHING_V5_CB_25.06.R" 

4) ANOVA, pairwise tests, emmeans via: "Coding/ANOVA_PAIRWISE_Script_and_Evaluation_Data" script and dataframes provided. 
