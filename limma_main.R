#This is an omplementation of LIMMA statistics on two group comparision as described in the paper http://www.sciencedirect.com/science/article/pii/S2212968515000069
#Author:Wasim Aftab
library(dplyr)
library(stringr)
library(MASS)
library(plotly)
library(limma)
library(qvalue)
library(htmlwidgets)
cat('\014')
rm(list = ls())
##Chdir to source dir
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
source("limma_helper_functions.R")
#Load the proteingroups file
myFilePath <- file.choose()
temp <- unlist(strsplit(myFilePath, "\\", fixed = TRUE))
proteingroups <-
  as.data.frame(read.table(myFilePath, header = TRUE, sep = "\t"))
###Kill the code if proteingroups does not contain crap columns#######################################
temp <-
  select(
    proteingroups,
    matches("(Reverse|Potential.contaminant|Only.identified.by.site)")
  )
if (!nrow(temp) * ncol(temp)) {
  stop("File error, It does not contain crap...enter another file with crap")
}
######################################################################################################
##Display data to faciliate choice of treatment and control
temp <- select(proteingroups, matches("(ibaq|lfq)"))
print(names(temp))

# #remove "+" identifeid rows from proteingroups##################################################
idx <- NULL
# temp <- as.matrix(select(proteingroups, matches("(Only.identified.by.site|Reverse|Potential.contaminant)")))
temp <- select(proteingroups, matches("(Only.identified.by.site|Reverse|Potential.contaminant)"))
for (i in 1:ncol(temp)){
  index <- which(unlist(!is.na(match(temp[,i], "+"))))
  idx <- union(idx, index)
}
proteingroups <- proteingroups[-idx, ] # removing indexed rows
# ################################################################################################

#Extrat Uniprot and gene symbols
Uniprot <- character(length = nrow(proteingroups))
Symbol <- character(length = nrow(proteingroups))
for (i in 1:nrow(proteingroups)) {
  temp <- as.character(proteingroups$Fasta.headers[i])
  splits <- unlist(strsplit(temp, '\\|'))
  Uniprot[i] <- splits[2]
  splits <- unlist(str_match(splits[3], "GN=(.*?) PE="))
  Symbol[i] <- splits[2]
}

#Extract required data for Limma
treatment <-
  readline('Enter treatment name(case insensitive) as it appeared in the iBAQ/LFQ column= ')
control <-
  readline('Enter control name(case insensitive) as it appeared in the iBAQ/LFQ column= ')
ibaq <- readinteger_binary('Enter 1 for iBAQ or 0 for LFQ= ')

if (ibaq) {
  temp <-
    select(proteingroups, matches(paste('^.*', "ibaq", '.*$', sep = '')))
  treatment_reps <- data_sanity_check(temp, 'treatment', treatment)
  control_reps <- select(temp, matches(control))
  control_reps <- data_sanity_check(temp, 'control', control)
  data <-
    cbind(treatment_reps,
          control_reps,
          select(proteingroups, matches("^id$")),
          Uniprot,
          Symbol)
} else {
  temp <-
    select(proteingroups, matches(paste('^.*', "lfq", '.*$', sep = '')))
  treatment_reps <- select(temp, matches(treatment))
  treatment_reps <- data_sanity_check(temp, 'treatment', treatment)
  control_reps <- select(temp, matches(control))
  control_reps <- data_sanity_check(temp, 'control', control)
  data <-
    cbind(treatment_reps,
          control_reps,
          select(proteingroups, matches("^id$")),
          Uniprot,
          Symbol)
}

#Impute data
print(names(data))
rep_treats <-
  readinteger("Enter the number of treatment replicates=")
rep_conts <- readinteger("Enter the number of control replicates=")
FC_Cutoff <- readfloat("Enter the fold change cut off=")
data_limma <- log2(as.matrix(data[c(1:(rep_treats + rep_conts))]))
data_limma[is.infinite(data_limma)] <- NA
nan_idx <- which(is.na(data_limma))

fit <- fitdistr(c(na.exclude(data_limma)), "normal")
mu <- as.double(fit$estimate[1])
sigma <- as.double(fit$estimate[2])
sigma_cutoff <- 6
new_width_cutoff <- 0.3
downshift <- 1.8
width <- sigma_cutoff * sigma
new_width <- width * new_width_cutoff
new_sigma <- new_width / sigma_cutoff
new_mean <- mu - downshift * sigma
imputed_vals_my = rnorm(length(nan_idx), new_mean, new_sigma)
data_limma[nan_idx] <- imputed_vals_my

##Limma main code
design <-
  model.matrix( ~ factor(c(rep(2, rep_treats), rep(1, rep_conts))))
colnames(design) <- c("Intercept", "Diff")
res.eb <- eb.fit(data_limma, design, Symbol)
Sig_FC_idx <-
  union(which(res.eb$logFC < (-FC_Cutoff)), which(res.eb$logFC > FC_Cutoff))
Sig_Pval_mod_idx <- which(res.eb$p.mod < 0.05)
Sig_Pval_ord_idx <- which(res.eb$p.ord < 0.05)
Sig_mod_idx <-  intersect(Sig_FC_idx, Sig_Pval_mod_idx)
Sig_ord_idx <-  intersect(Sig_FC_idx, Sig_Pval_ord_idx)
categ_Ord <- rep(c("Not Significant"), times = length(Symbol))
categ_Mod <- categ_Ord
categ_Mod[Sig_mod_idx] <- "Significant"
categ_Ord[Sig_ord_idx] <- "Significant"
dat <-
  cbind(
    res.eb,
    categ_Ord,
    categ_Mod,
    NegLogPvalMod = (-log10(res.eb$p.mod)),
    NegLogPvalOrd = (-log10(res.eb$p.ord))
  )

##Create plotly object and save plot as html
filename_mod <-
  readline('Enter a filename for limma plot= ')
filename_ord <-
  readline('Enter a filename for ordinary t-test plot= ')

display_plotly_figs(dat, FC_Cutoff, filename_mod, filename_ord)

##Save the data file
final_data <-
  cbind(select(proteingroups, matches("^id$")),
        Uniprot,
        Symbol,
        data_limma,
        dat)
final_data <- select(final_data, -matches("^gene$"))
filename_final_data <-
  readline('Enter a filename for final data= ')
write.table(
  final_data,
  paste(filename_final_data, '.tsv', sep = ''),
  sep = '\t',
  row.names = FALSE,
  col.names = TRUE
)