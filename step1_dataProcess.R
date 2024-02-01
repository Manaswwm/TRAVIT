#### this script is the first step in the microarray analysis pipeline ####
## this script serves the purpose of preparing the input files to feed into the consecutive steps

## the function will take in four values - variable (to be tested), values (of the variable),
## their names in the meta file and finally the metainfo file which is to be processed - flu or pneumo

step1_dataProcess = function(variable, values, meta_names, meta_info){
  
  #@print check
  print("Inside Step 1 - Data pre-processing .. ")
  
  ##shortlisting file names on the basis of the factors that are to be tested
  affymtx_files1 = meta_info[[1]]$microarray_id[meta_info[[1]][,variable] == values[1]]
  affymtx_files2 = meta_info[[2]]$microarray_id[meta_info[[2]][,variable] == values[2]]
  
  #extracting the file locations for those microarrayIDs
  affymtx_files1 = affymtx_filenames[grep(pattern = paste(affymtx_files1, collapse = "|"),
                                            x = affymtx_filenames)]
  affymtx_files2 = affymtx_filenames[grep(pattern = paste(affymtx_files2, collapse = "|"),
                                            x = affymtx_filenames)]

  #merging the file paths
  affy_filepaths = c(affymtx_files1, affymtx_files2)
  
  #@print check
  print("(Step 1) - Done making affy file names")
  
  #creating a batch-specific metainfo table
  batch_meta_info = data.frame(chp_files = affy_filepaths)
  batch_meta_info$group[batch_meta_info$chp_files %in% affymtx_files1] = meta_names[1]
  batch_meta_info$group[batch_meta_info$chp_files %in% affymtx_files2] = meta_names[2]
  batch_meta_info$id = str_match(batch_meta_info$chp_files, "0000-*(.*?).sst-rma-gene-full.chp")[,2]
  #batch_meta_info$sample_id = paste(batch_meta_info$id,"_",batch_meta_info$group, sep = "")
  
  #for the "id" column - some entries have a "_2" added in the entries - this information is captured in sample_id column already
  #removing this from "id" column as this hinders with adding microarray id and batch no
  batch_meta_info$microarray_id = str_remove(string = batch_meta_info$id, pattern = "_.*") #removing everything after "_" only if it exists, does not affect for whom this does not exist
  
  #### addition for making batch correction as per the microarray batches - suggestion by CC - 17012024 ####
  
  #adding the microarray id per row - here every row denotes a new microarray run (6 total for an individual)
  #the microarray ids are included in the chp_file names - the start from CRID and end with 4 numbers
  #the "id" column in batch_meta_info are a subset of the microarray id - using them to extract the microarray id from metainfo
  batch_meta_info$microarray_id = unlist(lapply(batch_meta_info$microarray_id, function(x){metainfo$microarray_id[grep(pattern = x, x = metainfo$microarray_id)]}))

  #adding the batch number now
  batch_meta_info = merge(batch_meta_info, metainfo[,c("batch_no", "microarray_id", "subj_id", "timepoint")], by = "microarray_id")

  #first adding leading zeros because CC's expr mtx has those
  batch_meta_info$subj_id = sprintf("%03d", batch_meta_info$subj_id)

  #substitution of sample id column - MJ17012024
  batch_meta_info$sample_id = paste(batch_meta_info$subj_id, "_", batch_meta_info$timepoint, sep = "")

  #removing sample_id that are mentioned in CC's presentation that do not meet the threshold criterion (slide 5)
  batch_meta_info = batch_meta_info[!batch_meta_info$sample_id %in% c("184_G","145_A","274_F","179_C","012_C","232_C","199_E","182_B","213_A"),]
  
  #making seperate columns for the two groups - binary logic - 1 is yes and 0 is no
  #first creating blank columns per variable value
  batch_meta_info[,meta_names[1]] = NA
  batch_meta_info[,meta_names[2]] = NA

  #second filling in the values as per their respective group identifiers
  batch_meta_info[,meta_names[1]][batch_meta_info$group == meta_names[1]] = 1
  batch_meta_info[,meta_names[1]][is.na(batch_meta_info[,meta_names[1]])] = 0
  batch_meta_info[,meta_names[2]][batch_meta_info$group == meta_names[2]] = 1
  batch_meta_info[,meta_names[2]][is.na(batch_meta_info[,meta_names[2]])] = 0
  
  #have to add the sample_id as rownames --> the column names of the dge_mtx_norm$counts should match the rownames of the batch_meta_info
  rownames(batch_meta_info) = batch_meta_info$sample_id
  
  #@ sanity check - removing any affyfilepaths if they are not in the batch_meta_info
  #this would happen if I have removed any affymtx files due to their deletion from CC's criteria (slide 5 from CC's presentation)
  affy_filepaths = affy_filepaths[affy_filepaths %in% batch_meta_info$chp_files]
  
  #### visualizing how much individuals am I analysing ####
  ind_distribution = ggplot(batch_meta_info, aes(x = group, fill = group))+geom_bar(width = 0.4)+
    xlab("Group")+ylab("Number of individuals")+
    scale_fill_manual("Groups", values = c("#0A0DAE", "#F57A00"))+
    theme_bw()+theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 16),
                     legend.title = element_text(size = 16), legend.text = element_text(size = 14),
                     axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 14))
  
  #saving to a file
  ggsave(plot = ind_distribution, filename = paste(backup_file, "/group_distribution.jpg", sep = ""), width = 8, height = 10, dpi = 300)
  
  #@print check
  print("(Step 1) - Done making batch-specific meta-info")
  
  #@print check
  print("(Step 1) - Returning the necessary information to the mother function")
  
  #returning back the information
  return(list(affy_filepaths = affy_filepaths, batch_meta_info = batch_meta_info))
  
}


# #randomly picking 20 from each
# affymtx_files1 = sample(affymtx_files1, size = 20, replace = FALSE)
# affymtx_files2 = sample(affymtx_files2, size = 20, replace = FALSE)
# 


# #adding microarray id so that I can add other metainfo to this
# batch_meta_info$microarray_id = paste("CRID-",str_match(batch_meta_info$chp_files, "CRID-*(.*?).sst-rma-gene-full.chp")[,2], sep = "")
# 
# #removing the numbers after "_" if they exist - these exist only in file names but not in metainfo dataframe
# batch_meta_info$microarray_id = unlist(lapply(batch_meta_info$microarray_id, function(x){strsplit(x, split = "_")[[1]][1]}))
# 
# #merging to add vital_age_group and gender -- ACHTUNG - directly using the main metainfo dataframe here
# batch_meta_info = merge(batch_meta_info, metainfo[,c("microarray_id", "vital_age_group", "gender")], by = "microarray_id")
  
  # }else if(fact == 2){
  #   
  #   ##shortlisting file names on the basis of the factors that are to be tested
  #   affymtx_files1 = meta_info[[1]]$microarray_id[meta_info[[1]][,variable[1]] == values[1] & meta_info[[1]][,variable[2]] == values[3]]
  #   affymtx_files2 = meta_info[[2]]$microarray_id[meta_info[[1]][,variable[1]] == values[2] & meta_info[[1]][,variable[2]] == values[3]]
  #   affymtx_files3 = meta_info[[1]]$microarray_id[meta_info[[1]][,variable[1]] == values[1] & meta_info[[1]][,variable[2]] == values[4]]
  #   affymtx_files4 = meta_info[[2]]$microarray_id[meta_info[[1]][,variable[1]] == values[2] & meta_info[[1]][,variable[2]] == values[4]]
  #   
  #   #extracting the file locations for those microarrayIDs
  #   affymtx_files1 = affymtx_filenames[grep(pattern = paste(affymtx_files1, collapse = "|"),
  #                                           x = affymtx_filenames)]
  #   affymtx_files2 = affymtx_filenames[grep(pattern = paste(affymtx_files2, collapse = "|"),
  #                                           x = affymtx_filenames)]
  #   affymtx_files3 = affymtx_filenames[grep(pattern = paste(affymtx_files3, collapse = "|"),
  #                                           x = affymtx_filenames)]
  #   affymtx_files4 = affymtx_filenames[grep(pattern = paste(affymtx_files4, collapse = "|"),
  #                                           x = affymtx_filenames)]
  #   
  #   # #randomly picking 20 from each
  #   # affymtx_files1 = sample(affymtx_files1, size = 20, replace = FALSE)
  #   # affymtx_files2 = sample(affymtx_files2, size = 20, replace = FALSE)
  #   # 
  #   #merging the file paths
  #   affy_filepaths = c(affymtx_files1, affymtx_files2, affymtx_files3, affymtx_files4)
  #   
  #   #@print check
  #   print("(Step 1) - Done making affy file names")
  #   
  #   #creating a batch-specific metainfo table
  #   batch_meta_info = data.frame(chp_files = affy_filepaths)
  #   
  #   #the first two affy file names represent a single comparison
  #   batch_meta_info$group[batch_meta_info$chp_files %in% affymtx_files1] = meta_names[1]
  #   batch_meta_info$group[batch_meta_info$chp_files %in% affymtx_files2] = meta_names[2]
  #   batch_meta_info$group[batch_meta_info$chp_files %in% affymtx_files3] = meta_names[3]
  #   batch_meta_info$group[batch_meta_info$chp_files %in% affymtx_files4] = meta_names[4]
  #   batch_meta_info$id = str_match(batch_meta_info$chp_files, "0000-*(.*?).sst-rma-gene-full.chp")[,2]
  #   batch_meta_info$sample_id = paste(batch_meta_info$id,"_",batch_meta_info$group, sep = "")
  #   
  #   # #adding microarray id so that I can add other metainfo to this
  #   # batch_meta_info$microarray_id = paste("CRID-",str_match(batch_meta_info$chp_files, "CRID-*(.*?).sst-rma-gene-full.chp")[,2], sep = "")
  #   # 
  #   # #removing the numbers after "_" if they exist - these exist only in file names but not in metainfo dataframe
  #   # batch_meta_info$microarray_id = unlist(lapply(batch_meta_info$microarray_id, function(x){strsplit(x, split = "_")[[1]][1]}))
  #   # 
  #   # #merging to add vital_age_group and gender -- ACHTUNG - directly using the main metainfo dataframe here
  #   # batch_meta_info = merge(batch_meta_info, metainfo[,c("microarray_id", "vital_age_group", "gender")], by = "microarray_id")
  #   
  #   #making seperate columns for the two groups - binary logic - 1 is yes and 0 is no
  #   #first creating blank columns per variable value
  #   batch_meta_info[,meta_names[1]] = NA
  #   batch_meta_info[,meta_names[2]] = NA
  #   
  #   #second filling in the values as per their respective group identifiers
  #   batch_meta_info[,meta_names[1]][batch_meta_info$group == meta_names[1]] = 1
  #   batch_meta_info[,meta_names[1]][is.na(batch_meta_info[,meta_names[1]])] = 0
  #   batch_meta_info[,meta_names[2]][batch_meta_info$group == meta_names[2]] = 1
  #   batch_meta_info[,meta_names[2]][is.na(batch_meta_info[,meta_names[2]])] = 0
  #   
  #   #making seperate columns for the two groups - binary logic - 1 is yes and 0 is no
  #   #first creating blank columns per variable value
  #   batch_meta_info[,meta_names[1]] = NA
  #   batch_meta_info[,meta_names[2]] = NA
  #   batch_meta_info[,meta_names[3]] = NA
  #   batch_meta_info[,meta_names[4]] = NA
  #   
  #   
  #   #second filling in the values as per their respective group identifiers
  #   batch_meta_info[,meta_names[1]][batch_meta_info$group == meta_names[1]] = 1
  #   batch_meta_info[,meta_names[1]][is.na(batch_meta_info[,meta_names[1]])] = 0
  #   batch_meta_info[,meta_names[2]][batch_meta_info$group == meta_names[2]] = 1
  #   batch_meta_info[,meta_names[2]][is.na(batch_meta_info[,meta_names[2]])] = 0
  #   batch_meta_info[,meta_names[3]][batch_meta_info$group == meta_names[3]] = 1
  #   batch_meta_info[,meta_names[3]][is.na(batch_meta_info[,meta_names[3]])] = 0
  #   batch_meta_info[,meta_names[4]][batch_meta_info$group == meta_names[4]] = 1
  #   batch_meta_info[,meta_names[4]][is.na(batch_meta_info[,meta_names[4]])] = 0
  #   
  #   #@print check
  #   print("(Step 1) - Done making batch-specific meta-info")
  #   
  #   #@print check
  #   print("(Step 1) - Returning the necessary information to the mother function")
  #   
  #   #returning back the information
  #   return(list(affy_filepaths = affy_filepaths, batch_meta_info = batch_meta_info))
  # }