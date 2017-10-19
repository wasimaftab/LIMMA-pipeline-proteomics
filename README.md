# LIMMA pipeline

![limma_pipeline_finite](https://user-images.githubusercontent.com/29901809/31774275-917afc3c-b4e5-11e7-900b-9ee2f36c8fb0.gif)

1. Implementation of LIMMA (Linear Models for Microarray Data) in a pipeline for two group comparision in proteomic experiment using an empirical Bayes method as described in the paper http://www.sciencedirect.com/science/article/pii/S2212968515000069 (Kammers, Kai, et al. "Detecting significant changes in protein abundance." EuPA open proteomics 7 (2015): 11-19.)
2. The pipeline is implemented in R programming language and you need following R libraries:
dplyr, stringr, MASS, plotly, limma, qvalue, htmlwidgets
3. Just run the lima_main.R code and select a MaxQuant outputted proteingroups.txt file (see inside Example folder). The helper functions must be in the same directory as the main. After successful run it will create two volcano plots and a txt file with data in the present working directory. One plot is using limma moderated pval and the other one using ordinary t-test.
4. The code will force you to provide as input correct names for treatment and control as it appears in the proteingroups file. Part of the name (case insensitive) should also be fine.
5. If columns 'Potential contaminant', 'Reverse' and 'Only.identified by site' is not present in the proteingroups file the code terminates with a message.
5. Tested on win7, Professional, 64 bit using Rstudio 1.1.383
