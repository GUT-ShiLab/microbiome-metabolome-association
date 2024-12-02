## Data_preprocess_for_comparison_methods

1.This section of code is provided in Jupyter Notebook format and facilitates data preprocessing for the comparison methods in Task 2. The specific preprocessing steps for each model are summarized in the table below. While some preprocessing tasks are already integrated into the model source code, this section is critical for bridging the input data with each model's requirements.

| Method    | Data Preprocessing Steps                                     |      |
| --------- | ------------------------------------------------------------ | ---- |
| Melonnpan | (1) Sample normalization: Normalize to relative abundance (2) Feature  selection: Retain features with average relative abundance > 0.01% in at  least 10% of samples (3) Feature transformation: Apply rntransform function  for microbial features and arcsine square root transformation for metabolite  features |      |
| ENVIM     | (1) Sample normalization: Normalize to relative abundance (2) Feature  selection: Retain microbial features with average relative abundance >  5×10^-5 and present in >90% of samples; retain metabolite features with  average relative abundance >10^-4 and present in >90% of samples (3)  Feature transformation: Apply rntransform function for microbial features and  Box-Cox transformation for metabolite features |      |
| MMINP     | (1) Sample normalization: Normalize to relative abundance (2) Feature  selection: Retain features with average relative abundance > 0.01% in at  least 10% of samples (3) Feature transformation: Apply Box-Cox transformation  and scaling for microbial and metabolite features, ensuring each feature has  a mean of 0 and variance of 1 |      |
| MiMeNet   | (1) Feature selection: Remove features present in less than 10% of  samples (2) Sample normalization: Apply centered log-ratio (CLR)  transformation to each sample |      |
| BiomeNED  | (1) Feature selection: Remove features with zero measurements in over  50% of samples (2) Sample normalization: Apply centered log-ratio (CLR)  transformation to each sample |      |
| mNODE     | Sample normalization: Apply centered log-ratio (CLR) transformation to  each sample |      |
| LOCATE    | (1) Sample normalization: Apply log normalization to data samples (2)  Feature transformation: Apply Z-score standardization to features, ensuring  each feature has a mean of 0 and variance of 1 |      |

2.The code processes three datasets: ERAWIJANTARI_GASTRIC_CANCER_2020, FRANZOSA_IBD_2019, and WANG_ESRD_2020. For each dataset, dedicated processing pipelines and results are provided for each method. Users can download the method source code, input the processed data, and perform model validation. Notably, the FRANZOSA_IBD_2019 dataset also includes preprocessing steps for some methods from Task 1.

## Experimental_results_analysis

1.This R-based code is designed for the method comparison pipeline in Task 1, covering feature selection for microbes and metabolites, sample filtering (PRISM cohort), CLR transformation, correlation analysis (Spearman, Sparse CCA, and O2PLS), differential metabolite testing, and result visualization. Users can follow the code step-by-step to replicate the results presented in the paper.

2.The code also includes the Friedman test and result visualization for Task 2 (friedman_test_8.R and Task2_training_plot_9.R), with all plotting data available for download.
