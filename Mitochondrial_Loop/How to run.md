## This is a pipeline to analyze proteiomic data for the paper Molecular wiring of a mitochondrial translational feedback loop (accepted in Molecular Cell)

* Implementation of LIMMA (Linear Models for Microarray Data), an empirical Bayes method for two group comparision in a proteomic experiment [1].
* The pipeline is implemented in R programming language and needed libraries are installed automatically
* Run the lima_main.R code and select a MaxQuant outputted proteingroups.txt file with either iBAQ/LFQ values (see Example folder). 
* Make sure you have mitochondrial proteins (identified with uniprot id) in your proteingroups otherwise code will exit midway
* The helper functions must be in the same directory as the main. 
* After successful run it will create volcano plots both in html and pdf formats and a tsv file containing final data inside a folder called  "Results" in the same directory where the main function is present. 
* One plot is using limma moderated statistics and the other one using ordinary t-test.
* In addition, the code will create files containng Outlier proteins present (if any) in control and treatment replicates with lfq/ibaq values, for example an outlier in treatment is defined as follows,
## Uniprot  Symbol  treat_1   treat_2   treat_3   contrl_1    contrl_2    contrl_3
## -------------------------------------------------------------------------------
## P25554   SGF29	  0	        0	      0        2810900	     0	        0
## -------------------------------------------------------------------------------
* The code will force you to provide correct names as input for treatment and control as it appears in the proteingroups file. Part of the name (case insensitive) should also be fine.
* If the columns 'Potential contaminant', 'Reverse' and 'Only.identified by site' is not present in the proteingroups.txt file the code terminates with a message.
* Tested on win7 Professional 64 bit, Rstudio 1.1.463, R version 3.6.1 and using the following packages
--------------------
Package		Version
--------------------
dplyr       0.8.3
stringr     0.4.0
MASS        7.3-5.4
plotly      4.9.0
htmlwidgets 0.3
limma       3.42.0
qvalue      2.8.0

### Areas for Improvements
* Add more checks for undesired inputs.

### Reference:
[1] Kammers, Kai, et al. "Detecting significant changes in protein abundance." EuPA open proteomics 7 (2015): 11-19.

