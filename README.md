# LIMMA pipeline

![limma_pipeline_finite](https://user-images.githubusercontent.com/29901809/31774275-917afc3c-b4e5-11e7-900b-9ee2f36c8fb0.gif)

* Implementation of LIMMA (Linear Models for Microarray Data), an empirical Bayes method for two group comparision in a proteomic experiment.[1]
* The pipeline is implemented in R programming language and you need following R libraries:
dplyr, stringr, MASS, plotly, limma, qvalue, htmlwidgets
* Run the lima_main.R code and select a MaxQuant outputted proteingroups.txt file with either iBAQ/LFQ values (see Example folder). 
* The helper functions must be in the same directory as the main. 
* After successful run it will create two volcano plots and a tsv file containing data in the present working directory. 
* One plot is using limma moderated pval and the other one using ordinary t-test.
* The code will force you to provide as input correct names for treatment and control as it appears in the proteingroups file. Part of the name (case insensitive) should also be fine.
* If the columns 'Potential contaminant', 'Reverse' and 'Only.identified by site' is not present in the proteingroups.txt file the code terminates with a message.
* Tested on win7, Professional, 64 bit using Rstudio 1.1.383

### Areas for Improvements
* A mechanism to install the needed R packages automatically.
* Add more checks for undesired input.

### Reference:
[1] Kammers, Kai, et al. "Detecting significant changes in protein abundance." EuPA open proteomics 7 (2015): 11-19.

