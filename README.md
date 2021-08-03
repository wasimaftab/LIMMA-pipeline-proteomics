# LIMMA pipeline

![limma_pipeline_infinite](https://user-images.githubusercontent.com/29901809/32143726-4a82b77c-bcae-11e7-8435-b447cc2274c0.gif)

* Implementation of LIMMA (Linear Models for Microarray Data), an empirical Bayes method for two group comparision in a proteomic experiment [1].
* The pipeline is implemented in R programming language and all the required packages will auto install when the script is run.
* Run the limma_main.R code by clicking **Source** in RStudio and select a MaxQuant outputted proteingroups.txt file (see Example folder). The helper functions must be in the same directory as the main. 
* There are two modes: either use full data or remove exclusive proteins before analysis. Proteins that have zero intensities in all the replicates of one group are defined as exclusive proteins. For example, an 'exclusively enriched' proteins in treatment(Cbp1) is defined as follows

 |Uniprot | Symbol  |Cbp1_1  |  Cbp1_2  |  Cbp1_3  |  contrl_1 |   contrl_2  |  contrl_3|
 |--------|---------|--------|----------|----------|-----------|-------------|----------|
 |P25554  |SGF29    |0       |     0    |    0     |   2810900 |    4903800  |     0    |
* Proteins that do not participate in two group comparison are saved in the results folder as tsv files with their iBAQ/LFQ intensities.
* Code further removes the proteins based on k out of N criteria: where **N** is the number of replicates in one group (treatment/control) and for each protein, **k** implies the desired number of non zero values out of **N** replicates. The code applies this criterion individually for each group and keeps only those proteins that satisfy it in both the groups simultaneously.
* Code will ask you to provide treatment and control names. It will guide by printing them on your RStudio console. For example, in the proteinsgroups file in the Example folder there are three treatment replicates iBAQ CA_1, iBAQ CA_2, iBAQ CA_3 and three control replicates iBAQ FA_1, iBAQ FA_2, iBAQ FA_3. So, if you input partial names i.e. ca or CA for the treatment in this example, the code is able to recognise the desired columns.
* In other words, code will force you to provide as input correct names for treatment and control as it appears in the proteingroups file. Part of the name (case insensitive) should also be fine.
* It will also ask users if they want to normalize their data prior to two-group comparison. There are two modes of normalization supported.

  1. **Normalise by subtracting median:** This method normalizes the protein intensities in each experiemnt by substrating the median 
      of the corresponding experiment.
      
  2. **Column wise median normalization of the data matrix:** This method calculates for each sample the median change (i.e. the difference    between the observed value and the row average) and subtracts it from each row. Missing values are ignored in the procedure. The method is based
on the assumption that a majority of the rows did not change. By default, all the rows are used for normalization but if we assume that the first 3 proteins are spike-ins then the call to median_normalization function needs to be modified as: median_normalization(data, spike_in_rows = 1:3)
     
* After successful run, it will create a volcano plots in html format and a tsv file containing final data inside a folder called "Results_timestamp" with the current system "timestamp" in the same directory where the limma_main.R file is present. 
* You can view the plot in any browser and save it as png by clicking camera icon in plot
* In addition, the code will save 'exclusively enriched' proteins (if any) in control and treatment/bait replicates with corresponding LFQ/iBAQ values in the "Results_timestamp" folder.
* One plot is using limma moderated statistics and the other one using ordinary t-test.
* Tested on Ubuntu 20.04, Windows 10 Enterprise 64 bit, RStudio >= 1.3.959 and R version >= 3.6.3 and using the following packages

| Package  | Version |
| ------------- | ------------- |
| dplyr         | >= 0.8.3  |
| stringr       | >= 0.4.0  |
| MASS       | >= 7.3-5.4 |
| plotly       | >= 4.9.0  |
| htmlwidgets       | >= 0.3 |
| limma       | >= 3.42.0  |
| qvalue       | >= 2.8.0  |

* We recommend using the same versions of the packages (atleast for plotly)

### Troubleshooting
* Ubuntu users might need to install **libssl-dev** and **libcurl4-openssl-dev** before running the script.

### Areas for Improvements
* Add more sanity checks for undesired inputs.
* Facilitate plotting in pdf/tiff formats.
* Add multigroup comparision.

### Reference
[1] Kammers, Kai, et al. "Detecting significant changes in protein abundance." EuPA open proteomics 7 (2015): 11-19.

### Cite Repository

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4050581.svg)](https://doi.org/10.5281/zenodo.4050581)

