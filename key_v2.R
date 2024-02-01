##########################################################################################

##### Tra-Ana-Vit #####

#### this is the second version of the center code key.R which processes affymtx data ####

#meta-information wrt the analysis goes here - 

## this script will serve as a base which will link to various functions each performing specific tasks
## required meta-data is stored in ../../microarray_metadata.txt
## file containing the raw version of the pipeline - affy_trials_with_limma.R (local drive of RIVM)

##pipeline maintainance - 12122023 || 15012024

##########################################################################################

#@clock check
start_clock = Sys.time()

#setting working directory
setwd("/rivm/r/IMI vaccination of elderly/11. Databases en manuscripten/Manas/manas_workfolder/microarray_pipeline/pipeline_scripts/")

#setting path to the affymetrix assay data
##! major change - 06102023 
#changing from CEL to CHP files - read notes in the whitepaper for TAC for more - in concise -
#we change from RMA norm to TAC norm - the chp files provided are already processed through TAC hence using them now
affymtx_filenames = list.files(path = "../../input_files/preprocessed_data/", pattern = "*.chp", full.names = TRUE, recursive = TRUE)

#setting the stats-based cutoffs - I do not use the lfc cutoff for now (from Christophes analyses plots it seems he too did not use this) - 12122023
p_val_cutoff = 0.05 #actual one - 0.05

#### ! Section 0.1 - Declaring libraries that will be used in the pipeline and sourcing functions ! ####

#### standard R libraries ####
library(ggplot2) #for plotting
library(ggrepel) #for adding labels to volcano plots
library(cowplot) #for plotting
library(umap) #for plotting UMAP plots
library(tidyr) # for manipulating data to making dataframes that can be plotted easily
library(dplyr) # for manipulating data to making dataframes that can be plotted easily - pt 2
library(tidyverse) #because no escaping tidyverse
library(stringr) #for string manipulation
library(kableExtra) #for making tables
library(readxl) #for importing xlsx files
#library(biomaRt) #for converting between probe id and gene id

#### libraries specific for transcrioptomics analyses and making plots ####
library(limma) #the central package for performing DE analysis - also mentioned in Christophe's WP2 presentation
#library(affy) #reading the CEL files from affymetrics
library(affxparser) #reading the chp files
#library(gcrma) #adjusting background intensities in the affy mdata which could include optical noise and non-specific binding
#gcrma uses information on sequences to identify the non-specific binding of transcripts to probes
#to use this as a factor for normalizing
library(oligo) #reading CEL files pt.2 - directly downloads required annotations
library(plotly) #for plotting alongwith ggplot
library(gplots) #for making heatmaps
library(gprofiler2) #for performing gene ontologies
library(msigdbr) #alternative to gene ontology - in the presentation from Christophe
library(clariomshumantranscriptcluster.db) #affymtx annotation pt. 1
library(pd.clariom.s.human.ht) #affymtx annotation pt. 2
library(sva) #correcting batch effects
library(variancePartition) #fitting linear models with repeated measures
#library(edgeR) #for filtering and normalizing counts
library(tmod) #for performing GO or similar analyses - tutorial - https://tmod.online/articles/user_manual/tmod_user_manual.html
library(BloodGen3Module) #for plotting the blood transcriptome modules
library(enrichplot) #plotting enrichment analysis
library(DOSE) #making enrichment analysis
library(pheatmap) #for making heatmaps
library(rWikiPathways) #for performing pathway enrichment analysis
library(clusterProfiler) #for performing pathway enrichment analysis
library(pathfindR) #for making pathway enrichment analysis

#####
#### sourcing relevant libraries ####

source("pre_process_scripts/preprocess_metainfo.R")
source("pre_process_scripts/preprocess_immunotype.R")
source("step1_dataProcess.R")
source("step2_processExpression.R")
source("step2_1_processExpressionC.R")
source("step3_estimateDEG.R")
source("step4_makeVolcanoHeatplots.R")
source("step5_tmodDGEmods.R")
source("step6_diseaseAnnotate.R")
source("step7_pathwayAnnotate.R")
source("step8_gseaAnalysis.R")

#### ! Section 0.2 - Providing meta-data that filters data that is to be analyses - EDIT BEFORE EXECUTING ! ####
#### importing different types of meta data and creating a single dataframe for metadata ####

### - main metainfo - to run everytime - ###

#calling function that will directly read pieces of information from file paths and return a metainfo
#I include the responder data (267 X 6) and day sampling data (150 X 22) (Day 1 or Day 2) within the main metainfo (897 X 9)
#On merging the responder data - I lose two individuals (ID - 118 and 281) - hence for the purpose of this analyses these two IDs are absent
#final dimension of metainfo - (885 X 14)
metainfo = preprocess_metainfo(path_to_metainfo = "../../input_files/microarray_metadata.txt",
                               path_to_responder_data = "../../../microarray_data/dual response scores voor Manas.xlsx",
                               path_to_sampling_data = "../../../day_1_2_vaccination/VITAL_Participant information Microarray dataset_20231207.xlsx")

### OPTIONAL -  immunotype data -  to run *only when looking at immunotypes* - ###

## since 12 individuals do not have immunotype assignment - running this everytime would result in loss of data ##
#calling function that will merge the immunotype data to the main metainfo
# metainfo = preprocess_immunotypes(metainfo = metainfo,
#                                   path_to_immunotype = "../../../immunotype_data/VITAL_Participant information Microarray dataset_20231221_immunotypes.xlsx")

### done passing metainfo through functions - continuing with next steps ###

#making metainfo dataframes seperately for the flu and pneumo data - based on the timepoints
#remember - in the dataframe metainfo - every row is a single timepoint per individual
metainfo_flu = metainfo[metainfo$timepoint %in% c("A", "B", "C"),]
metainfo_pneumo = metainfo[metainfo$timepoint %in% c("F", "E", "G"),]

## !ACHTUNG - explaination of the timestamps - A = T0; B = T2; C = T3; E = T5; F = T6; G = T7 (from Elske's email)
## !ACHTUNG - PLEASE NOTE - Timestamp A,B and C are for Flu; D,E and F are for PCV13

#### formula to be tested ####

## this formula could be written as follows for different variable types:
## Responder category - ~(res3 - res1); Age group - ~(MA - YA); Timepoint - ~(T6 - T5)
## immunotype category - ~(imm_8 - imm_1)
form = ~(T6 - T5)

#### data subsetting - usually used when data complexity is two fold - Ex: res 4 vs res 1 for T2 vs T0 ####

## subsetting data to look at specific timepoints, age, responder category ##
#subsetting to keep timepoint
#metainfo_flu = metainfo_flu[metainfo_flu$timepoint == "A",]
#metainfo_pneumo = metainfo_pneumo[metainfo_pneumo$timepoint == "E",]

#subsetting to keep only individuals who were vaccinated on Day 1 or Day 2
#metainfo_flu_timepoint = metainfo_flu[metainfo_flu$sampling_day_flu == 2,]
metainfo_pneumo_timepoint = metainfo_pneumo[metainfo_pneumo$sampling_day_pcv13 == 2,]

#subsetting to keep only responders
# metainfo_flu = metainfo_flu[metainfo_flu$res_flu == 4,]
# metainfo_flu_timepoint = metainfo_flu_timepoint[metainfo_flu_timepoint$res_flu == 4,]
# metainfo_pneumo = metainfo_pneumo[metainfo_pneumo$res_pcv13 == 4,]
# metainfo_pneumo_timepoint = metainfo_pneumo_timepoint[metainfo_pneumo_timepoint$res_pcv13 == 4,]

#subsetting to keep age groups
# metainfo_flu_timepoint = metainfo_flu_timepoint[metainfo_flu_timepoint$vital_age_group == "Older adults",]
# metainfo_flu = metainfo_flu[metainfo_flu$vital_age_group == "Older adults",]
# metainfo_pneumo_timepoint = metainfo_pneumo_timepoint[metainfo_pneumo_timepoint$vital_age_group == "Young adults",]
# metainfo_pneumo = metainfo_pneumo[metainfo_pneumo$vital_age_group == "Young adults",]
# metainfo_flu_timepoint = metainfo_flu_timepoint[metainfo_flu_timepoint$age_group_t0 %in% c("80-84", "85+"),]
# metainfo_flu = metainfo_flu[metainfo_flu$age_group_t0 %in% c("80-84", "85+"),]
# metainfo_pneumo_timepoint = metainfo_pneumo_timepoint[metainfo_pneumo_timepoint$age_group_t0 %in% c("75-79", "80-84"),]
# metainfo_pneumo = metainfo_pneumo[metainfo_pneumo$age_group_t0 %in% c("75-79", "80-84"),]

#subsetting to keep only individuals for specific immunotypes
#metainfo_flu = metainfo_flu[metainfo_flu$immtype == "imm_1",]
#metainfo_pneumo = metainfo_pneumo[metainfo_pneumo$immtype == 2,]

#### declaring variable to be analysed ####

## the name of variable could be written as follows for different variable types:
## Responder category - "res_pcv13"; Age group - "vital_age_group"; Timepoint - "timepoint" 
## immunotype category - "immtype"
variable = "timepoint" ##### !! variable has to be what I give in as "values"

#### values of the variable for comparing ####

## the name of values could be written as follows for different variable types:
## Responder category - c("1","2","3","4"); Age group - c("Young adults","Middle aged adults", "Older adults"); 
## Timepoint - c("A","B","C","E","F","G"); immunotype - c(imm_1, imm_2, imm_3 .. imm_8)
values = c("F", "E") #mind the order in which this is supplied

#### metanames that are to be given to the variables ####

## the metaname of values could be written as follows for different variable types:
## Responder category - c("res1","res2","res3","res4"); Age group - c("YA","MA", "OA"); 
## Timepoint - c("T0","T2","T3","T5","T6","T7"); immunotype - c(imm_1, imm_2, imm_3 .. imm_8)
meta_names = c("T6", "T5") #mind the order in which this is supplied

#### ! Section 0.3 - Final preparations for executing the steps within the pipeline ! ####

#supplying metainfo dataframes for the two comparing variables as list
#meta_info = list(metainfo_flu, metainfo_flu) #mind the order in which this is supplied
meta_info = list(metainfo_pneumo_timepoint, metainfo_pneumo)

#adding meta-info to file name for responder category
category = "samplingday_2"

#vaccine that is being analysed - pneumo or flu
vaccine_type = "pcv13"

#type of analysis that I want to do - timepoint, vital_age_group, responder
analysis_type = "timepoint"

#listing the backup folder in which all the intermediate files will be saved
backup_file = paste("../pipeline_backup_files/", analysis_type, "/", vaccine_type, 
                    "/", category, "_", meta_names[1], "_", meta_names[2], sep = "")

#checking if the backup file is created - if not then creating it to store outputs
if(!dir.exists(paste(backup_file))){
  
  #creating the directory if it does not exist
  dir.create(backup_file, recursive = TRUE)
  
  #print check
  print(paste("Created backup directory : ",backup_file,sep = ""))
  
}else{print(paste("Found backup directory : ",backup_file,sep=""))}


#### ! Section 1 - primary data preparation ! ####

## explanation of the timestamps - A = T0; B = T2; C = T3; E = T5; F = T6; G = T7 (from Elske's email)
## PLEASE NOTE - Timestamp A,B and C are for Flu; E,F and G are for PCV13

#giving a call to the function that takes care of processing the primary input data
#here the return is a list of two - this first is the list of all the affymtx file paths and
#the second is the batch-specific meta_info

#checking if the information is already present
if(!file.exists(paste(backup_file,"/affy_preprocess", sep = ""))){
  
  #call point for function to get the processed data
  affy_preprocess = step1_dataProcess(variable = variable, values = values, meta_names = meta_names, meta_info = meta_info)
  
  #saving the file in a backup directory
  save(affy_preprocess, file = paste(backup_file,"/affy_preprocess", sep = ""))
}

#loading the preprocessed data here
load(paste(backup_file,"/affy_preprocess", sep = ""))


################################################################################
##### --- processing expression - choose one of the two ways --- #####


#### Section 2.1 - reading in the affymtx files, normalizing and creating their expression matrices ####

#checking if the information is already present
# if(!file.exists(paste(backup_file,"/expr_mtrx", sep = ""))){
# 
#   #giving call to the function that takes care of reading and prcoessing the expression matrices
#   #here the return is exclusively the expression values normalized first and then corrected for batch effects using ComBat
#   expr_mtx = step2_processExpression(affy_filepaths = affy_preprocess$affy_filepaths,
#                                      batch_meta_info = affy_preprocess$batch_meta_info)
# 
#   #saving the file in a backup directory
#   save(expr_mtx, file = paste(backup_file,"/expr_mtrx", sep = ""))
# }
# 
# #loading the preprocessed data here
# load(paste(backup_file,"/expr_mtrx", sep = ""))

#### Section 2.2 - reading in the affymtx files, normalizing and creating their expression matrices - using CC exp values ####

#checking if the information is already present
if(!file.exists(paste(backup_file,"/expr_mtrx", sep = ""))){

  #giving call to the function that takes care of reading and prcoessing the expression matrices
  #here the return is exclusively the expression values normalized first and then corrected for batch effects using ComBat
  expr_mtx = step_2_1_processExpressionC(affy_filepaths = affy_preprocess$affy_filepaths,
                                     batch_meta_info = affy_preprocess$batch_meta_info)

  #saving the file in a backup directory
  save(expr_mtx, file = paste(backup_file,"/expr_mtrx", sep = ""))
}

#loading the preprocessed data here
load(paste(backup_file,"/expr_mtrx", sep = ""))


##### --- processing expression done --- #####

################################################################################



#### Section 3 - estimating the DGEs ####

#checking if the information is already present
if(!file.exists(paste(backup_file,"/dge_info", sep = ""))){
  
  #giving call to the function that takes care of estimating the DGEs
  #here the return is the table containing the DGEs, the ebayes table (the second used for making heatmap plots)
  #and the voom weights
  dge_info = step3_estimateDEG(expr_mtx = expr_mtx, batch_meta_info = affy_preprocess$batch_meta_info, form = form)
  
  #saving the file in a backup directory
  save(dge_info, file = paste(backup_file,"/dge_info", sep = ""))
}

#loading the preprocessed data here
load(paste(backup_file,"/dge_info", sep = ""))

#### Section 4 - constructing the volcano, PCA and heatmap plots ####
 
#giving a call to the function that takes care of reading the tophit DGE table and making the plots
#here the return the volcano plot - heatmap could not be stored as a variable for replotting
step4_makeVolcanoHeatplots(dge_tophits = dge_info$dge_tophits, diff_exp_genes = dge_info$diff_exp_genes,
                           batch_meta_info = affy_preprocess$batch_meta_info, expr_mtx = expr_mtx)

#### Section 5 -  extracting the gene modules using tmod ####

#giving a call to the function that calculates the gene modules
#no return - the function plots modules of DGE
step5_tmodDGEmods(dge_tophits = dge_info$dge_tophits)

#### Section 6 - making gene annotations from the perspective of biological pathways and terms (GO) ####
##searching for term enrichment for up- and down-regulated genes

#giving a call to the function that performs gene ontology
#no return - the function makes ontology plots
step7_pathwayAnnotate(dge_tophits = dge_info$dge_tophits)

##### Section 7 - making GSEA analysis from the perspective of immunological and hallmark gene sets ####

#giving a call to the function that performs gsea
#no-return - the function makes gsea plots and tables
step8_gseaAnalysis(dge_tophits = dge_info$dge_tophits)

#@time check
stop_clock = Sys.time()

#@time check - taking difference of start and stop time
total_time = difftime(stop_clock, start_clock, units = "mins")

#@time taken to complete
total_time


################################################

#p-val and lfc-cutoff for disease and gene ontology analysis?

# #applying a lfc and p value cutoff
# dge_tophits = dge_info$dge_tophits[dge_info$dge_tophits$P.Value < 0.05 & abs(dge_info$dge_tophits$logFC) > 0.3,]

# #### following a tutorial from rWikiPathways - https://bioconductor.org/packages/devel/bioc/vignettes/rWikiPathways/inst/doc/
# 
# #extracting the gene names of the up and down regulated genes from dge_tophits
# upreg_tophits = dge_tophits[dge_tophits$logFC > 0.3, "gene_symbol"]
# downreg_tophits = dge_tophits[dge_tophits$logFC < -0.3, "gene_symbol"]
# 
# #extracting all the gene symbols from dge_tophits - this serves as the background gene pool
# bkgrd_genes = dge_info$dge_tophits[,"gene_symbol"]
# 
# #adding a conversion step to get entrez ids
# regulated_entrezids = mapIds(clariomshumantranscriptcluster.db, keys = c(upreg_tophits, downreg_tophits),
#                              column = 'ENTREZID', keytype = 'SYMBOL')
# allgenes_entrezids = mapIds(clariomshumantranscriptcluster.db, keys = bkgrd_genes,
#                              column = 'ENTREZID', keytype = 'SYMBOL')
# 
# #downstream processing the names to make them into a dataframe
# regulated_entrezids = as.data.frame(regulated_entrezids)
# regulated_entrezids$gene_symbol = rownames(regulated_entrezids)
# rownames(regulated_entrezids) = NULL
# 
# allgenes_entrezids = as.data.frame(allgenes_entrezids)
# allgenes_entrezids$gene_symbol = rownames(allgenes_entrezids)
# rownames(allgenes_entrezids) = NULL
# 
# #removing NA - rows with no entrez gene ids
# regulated_entrezids = regulated_entrezids[!is.na(regulated_entrezids$regulated_entrezids),]
# allgenes_entrezids = allgenes_entrezids[!is.na(allgenes_entrezids$allgenes_entrezids),]
# 
# #making gene ontology analysis for up-regulated genes
# #ont - types of ontologies (MF - molecular function, "BP" - biological process, "CC" - cellular component)
# upreg_go = enrichGO(gene = regulated_entrezids$regulated_entrezids[regulated_entrezids$gene_symbol %in% upreg_tophits],
#                     universe = allgenes_entrezids$allgenes_entrezids, OrgDb = org.Hs.eg.db, ont = "MF", pAdjustMethod = "fdr", 
#                     pvalueCutoff = p_val_cutoff, readable = TRUE)
# 
# #making a network of gene ontologies
# goplot(upreg_go)
# 
# 
# #making wikipathways analysis
# ewp.up = clusterProfiler::enrichWP(gene = regulated_entrezids$regulated_entrezids[regulated_entrezids$gene_symbol %in% upreg_tophits],
#   universe = allgenes_entrezids$allgenes_entrezids,
#   organism = "Homo sapiens",
#   pAdjustMethod = "fdr",
#   pvalueCutoff = 0.1, #p.adjust cutoff; relaxed for demo purposes
#   )
# 
# ewp.up <- DOSE::setReadable(ewp.up, org.Hs.eg.db, keyType = "ENTREZID")
# barplot(ewp.up)
# #### Step 6 - making gene annotations from the perspective of disease networks ####
# 
# #giving a call to the function that performs disease ontology
# #no return - the function makes ontology plots
# step6_diseaseAnnotate(dge_tophits = dge_info$dge_tophits)

#### commenting this section as I already cover this in section 7
