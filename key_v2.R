#### this is the second version of the center code key.R which processes affymtx data ####

## this script will serve as a base which will link to various functions each performing specific tasks
## required meta-data is stored in ../../microarray_metadata.txt
## file containing the raw version of the pipeline - affy_trials_with_limma.R (local drive of RIVM)
## explaination of the timestamps - A = T0; B = T2; C = T3; E = T5; F = T6; G = T7 (from Elske's email)
## PLEASE NOTE - Timestamp A,B and C are for Flu; D,E and F are for PCV13


#### importing relevant libraries ####

library(limma) #the central package for performing DE analysis - also mentioned in Christophe's WP2 presentation
library(affy) #reading the CEL files from affymetrics
#library(gcrma) #adjusting background intensities in the affy mdata which could include optical noise and non-specific binding
#gcrma uses information on sequences to identify the non-specific binding of transcripts to probes
#to use this as a factor for normalizing
library(oligo) #reading CEL files pt.2 - directly downloads required annotations
library(ggplot2) #for plotting
library(tidyr) # for manipulating data to making dataframes that can be plotted easily
library(dplyr) # for manipulating data to making dataframes that can be plotted easily - pt 2
library(plotly) #for plotting alongwith ggplot
library(biomaRt) #for converting between probe id and gene id
library(gplots) #for making heatmaps
library(gprofiler2) #for performing gene ontologies
library(msigdbr) #alternative to gene ontology - in the presentation from Christophe
library(clariomshumantranscriptcluster.db)
library(stringr) #for string manipulation
library(sva) #correcting batch effects
library(variancePartition) #fitting linear models with repeated measures
library(edgeR) #for filtering and normalizing counts
library(tmod) #for performing GO or similar analyses - tutorial - https://tmod.online/articles/user_manual/tmod_user_manual.html
library(cowplot) #for plotting
library(BloodGen3Module) #for plotting the blood transcriptome modules

#### sourcing relevant libraries ####
source("step1_dataProcess.R")
source("step2_processExpression.R")
source("step3_estimateDEG.R")
source("step4_makeVolcanoHeatplots.R")
source("step5_tmodDGEmods.R")

#### initial steps ####

#importing the metadata file
metainfo = read.delim(file = "../../microarray_metadata.txt", sep = "\t", header = TRUE)

#listing names of the raw affymtx files
affymtx_filenames = list.files(path = "../../sample_data/raw/", pattern = "*.CEL", full.names = TRUE, recursive = TRUE)

#listing the backup folder in which all the intermediate files will be saved
#backup_folder = "../../backup_young_old_timeC/"

#making metainfo dataframes seperately for the flu and pneumo data
metainfo_flu = metainfo[metainfo$timepoint %in% c("A", "B", "C"),]
metainfo_pneumo = metainfo[metainfo$timepoint %in% c("D", "E", "F"),]

#### primary data preparation ####

#giving a call to the function that takes care of processing the primary input data
#here the return is a list of two - this first is the list of all the affymtx file paths and
#the second is the batch-specific meta_info
affy_preprocess = step1_dataProcess(variable = "timepoint", values = c("A", "B"), meta_names = c("T0", "T2"), meta_info = metainfo_flu)

#seperating the two outputs from the previous function
affy_filepaths = affy_preprocess[[1]]
batch_meta_info = affy_preprocess[[2]]

#dustbin - position 1
rm(affy_preprocess, metainfo, affymtx_filenames)


#### reading in the affymtx files, normalizing and creating their expression matrices ####

#giving call to the function that takes care of reading and prcoessing the expression matrices
#here the return is exclusively the expression values normalized first and then corrected for batch effects using ComBat
expr_mtrx = step2_processExpression(affy_filepaths = affy_filepaths, batch_meta_info = batch_meta_info)


#### estimating the DGEs ####

#giving call to the function that takes care of estimating the DGEs
#here the return is the table containing the DGEs, the ebayes table (the second used for making heatmap plots)
#and the voom weights
dge_info = step3_estimateDEG(expr_mtrx = expr_mtrx, batch_meta_info = batch_meta_info, form = ~(T2 - T0))

#seperating the two outputs
dge_tophits = dge_info[[1]]
dge_dream_fit_ebayes = dge_info[[2]]
dge_voom_weights = dge_info[[3]]


#### constructing the volcano and heatmap plots ####
 
#giving a call to the function that takes care of reading the tophit DGE table and making the plots
#here the return the volcano plot - heatmap could not be stored as a variable for replotting
volc_plot = step4_makeVolcanoHeatplots(dge_tophits = dge_tophits, batch_meta_info = batch_meta_info, cutoff = c(0.05, 0.05), dge_dream_fit_ebayes = dge_dream_fit_ebayes,
                                   dge_voom_weights = dge_voom_weights)


#### extracting the gene modules using tmod ####

#giving a call to the function that calculates the gene modules
#no return - the function plots modules of DGE
step5_tmodDGEmods(dge_tophits = dge_tophits, cutoff = c(0.05,0.05))







