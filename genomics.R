library(ggplot2)
library(tidyverse)
library(magrittr)
library(dplyr)
library(pheatmap)

# Bioconductor #
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Biobase")

###################################################################################################################
# Reproduction of the results found in                                                                            #
# Maurizio Callari, Vera Cappelletti Francesca D'Aiuto, Valeria Musella, Antonio Lembo, Fabien Petel, Thomas Karn,# 
# Takayuki Iwamoto, Paolo Provero, Maria Grazia Daidone, Luca Gianni, Giampaolo Bianchini;                        #
# Subtype-Specific Metagene-Based Prediction of Outcome after Neoadjuvant and Adjuvant Treatment in Breast Cancer.#
# Clin Cancer Res 15 January 2016; 22 (2): 337-345. https://doi.org/10.1158/1078-0432.CCR-15-0757                 #
###################################################################################################################

# Datas were downloaded from https://drive.switch.ch/index.php/s/YRxKQBSPOb8xonk #

# Path #
path <- rstudioapi::getActiveDocumentContext()$path
Encoding(path) <- "UTF-8"
setwd(dirname(path))
getwd()

rm(list=ls())

#Functions#
source("utility.R")
# Set Seed so that same sample can be reproduced in future also #
set.seed(101) 

# Data #
chemo <- load("data\\Chemo.rda")
generic <- load("data\\Generic.rda")
prognostic <- load("data\\Prognostic.rda")
tam <- load("data\\TAM.rda")

# Extracting samples IDs and dmfs, as well as expression levels of the genes (exp_mat) #
exp_mat_prog<- as.data.frame(assayData(prognostic_normalized)$exprs)
exp_mat_generic <- as.data.frame(assayData(generic_normalized)$exprs)
exp_mat_tam<- as.data.frame(assayData(tam_normalized)$exprs)
exp_mat_chemo <- as.data.frame(assayData(chemo_normalized)$exprs)

pheno_prog <- get_pheno_data(prognostic_normalized)
pheno_generic <- get_pheno_data(generic_normalized)
pheno_tam <- get_pheno_data(tam_normalized)
pheno_chemo <- get_pheno_data(chemo_normalized)

# Removing unused variables #
rm(generic,chemo,prognostic,tam)
gc()

#################################################################
# Identification of ER & HER status, based on paper's metagenes #
#################################################################

erherstatus_generic<-ER_HER_status(exp_mat_generic)
erherstatus_prog<-ER_HER_status(exp_mat_prog)
erherstatus_tam<- ER_HER_status(exp_mat_tam)
erherstatus_chemo<-ER_HER_status(exp_mat_chemo)

########################################################
# CLUSTERING                                           #
# This may take some time and saturate memory          #
########################################################

source("clustering.R")

# Computing the time of dmfs. Every value above 5 is put to 5
# Distant metastasis detection is marked as a 1 event if occured below 5 years.
# 
proliferation <- pheno_prog[rownames(erherstatus_prog),]
proliferation$event[proliferation$t.dmfs.norm.years >=5] <- 0
proliferation$event[proliferation$t.dmfs.norm.years <5] <- 1

proliferation$time[proliferation$t.dmfs.norm.years >=5] <- 5
proliferation <- proliferation%>% 
  mutate(time = coalesce(time,t.dmfs.norm.years))

pheno_prog$event[pheno_prog$t.dmfs.norm.years >=5] <- 0
pheno_prog$event[pheno_prog$t.dmfs.norm.years <5] <- 1
pheno_prog$time[pheno_prog$t.dmfs.norm.years >=5] <- 5

pheno_prog<- pheno_prog%>% 
  mutate(time = coalesce(time,t.dmfs.norm.years))

pheno_tam$event[pheno_tam$t.dmfs.norm.years >=5] <- 0
pheno_tam$time[pheno_tam$t.dmfs.norm.years >=5] <- 5
pheno_tam$event[pheno_tam$t.dmfs.norm.years <5] <- 1
pheno_tam<- pheno_tam%>% 
  mutate(time = coalesce(time,t.dmfs.norm.years))


pheno_chemo$ER <- expr_means_chemo$ER_status
pheno_chemo$HER2 <- expr_means_chemo$HER2_status

pheno_chemo$event[pheno_chemo$t.dmfs.norm.years >=5] <- 0
pheno_chemo$time[pheno_chemo$t.dmfs.norm.years >=5] <- 5
pheno_chemo$event[pheno_chemo$t.dmfs.norm.years <5] <- 1
pheno_chemo<- pheno_chemo%>% 
  mutate(time = coalesce(time,t.dmfs.norm.years))

