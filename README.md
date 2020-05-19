# LIMMA pipeline

![limma_pipeline_infinite](https://user-images.githubusercontent.com/29901809/32143726-4a82b77c-bcae-11e7-8435-b447cc2274c0.gif)

* Implementation of LIMMA (Linear Models for Microarray Data), an empirical Bayes method for two group comparision in a proteomic experiment [1].
* The pipeline is implemented in R programming language and all the required packages will auto install when the script is run.
* Run the limma_main.R code and select a MaxQuant outputted proteingroups.txt file (see Example folder). 
* There are two modes: either use full data or remove outliers before analysis. Proteins that have zero intensities in all the replicates of treatment or control are defined as outliers. They do not participate in two group comparision and wrote in seperate excel files in the Results folder.
* Code will ask you to provide treatment and control names. It will guide by printing them on your RStudio console. For example, in the proteinsgroups file in the Example folder there are three treatment replicates iBAQ CA_1, iBAQ CA_2, iBAQ CA_3 and three control replicates iBAQ FA_1, iBAQ FA_2, iBAQ FA_3. So, if you input partial names i.e. ca or CA for the treatment in this example, the code is able to recognise the desired columns.
* In other words, code will force you to provide as input correct names for treatment and control as it appears in the proteingroups file. Part of the name (case insensitive) should also be fine.
* The helper functions must be in the same directory as the main. 
* After successful run it will create volcano plots in html format and a tsv file containing final data inside a folder called  "Results" in the same directory where the main function is present. 
* One plot is using limma moderated statistics and the other one using ordinary t-test.
* If the columns 'Potential contaminant', 'Reverse' and 'Only.identified by site' is not present in the proteingroups.txt file the code terminates with a message.
* Tested on Ubuntu 18.04, Windows 10 Enterprise 64 bit, RStudio >= 1.1.383 and R version >= 3.3.1
* Ubuntu users need to install **libssl-dev** and **libcurl4-openssl-dev** before running the script.

### Areas for Improvements
* Add more sanity checks for undesired inputs.
* Facilitate plotting in pdf/tiff formats.
* Add multigroup comparision.

### Reference
[1] Kammers, Kai, et al. "Detecting significant changes in protein abundance." EuPA open proteomics 7 (2015): 11-19.

### Cite Repository

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3834635.svg)](https://doi.org/10.5281/zenodo.3834635)
