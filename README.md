1. # LIMMA-pipeline-proteiomics
Implementation of LIMMA pipeline in R programming language for two group comparision in proteomics data as described in the paper http://www.sciencedirect.com/science/article/pii/S2212968515000069
2. You need following R libraries and a proteingroups file (see inside Example forlder)
dplyr, stringr, MASS, plotly, limma, qvalue, htmlwidgets
3. Just run the lima_main.R code. The helper functions must be in the same directory as the main. After successful run it will create two volcano plots and a txt file with data in the present working directory. One plot is using limma moderated pval and the other one using ordinary t-test.
4. The code will force you to provide as input correct names for treatment and control as it appers in the proteingroups file. Part of the name (case insensitive) should also be fine.
5. If columns 'Potential contaminant', 'Reverse' and 'Only.identified by site' is not present in the proteingroups file the code terminates with a message.

