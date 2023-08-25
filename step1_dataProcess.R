#### this script is the first step in the microarray analysis pipeline ####
## this script serves the purpose of preparing the input files to feed into the consecutive steps

## the function will take in four values - variable (to be tested), values (of the variable),
## their names in the meta file and finally the metainfo file which is to be processed - flu or pneumo

step1_dataProcess = function(variable, values, meta_names, meta_info){
  
  #@print check
  print("Inside Step 1 - Data pre-processing .. ")
  
  ##shortlisting file names on the basis of the factors that are to be tested
  affymtx_files1 = meta_info$microarray_id[meta_info[,variable] == values[1]]
  affymtx_files2 = meta_info$microarray_id[meta_info[,variable] == values[2]]
  
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
  batch_meta_info = data.frame(cel_files = affy_filepaths)
  batch_meta_info$group[batch_meta_info$cel_files %in% affymtx_files1] = meta_names[1]
  batch_meta_info$group[batch_meta_info$cel_files %in% affymtx_files2] = meta_names[2]
  batch_meta_info$id = str_match(batch_meta_info$cel_files, "0000-*(.*?).CEL")[,2]
  batch_meta_info$sample_id = paste(batch_meta_info$id,"_",batch_meta_info$group, sep = "")
  
  #making seperate columns for the two groups - binary logic - 1 is yes and 0 is no
  #first creating blank columns per variable value
  batch_meta_info[,meta_names[1]] = NA
  batch_meta_info[,meta_names[2]] = NA
  
  #second filling in the values as per their respective group identifiers
  batch_meta_info[,meta_names[1]][batch_meta_info$group == meta_names[1]] = 1
  batch_meta_info[,meta_names[1]][is.na(batch_meta_info[,meta_names[1]])] = 0
  batch_meta_info[,meta_names[2]][batch_meta_info$group == meta_names[2]] = 1
  batch_meta_info[,meta_names[2]][is.na(batch_meta_info[,meta_names[2]])] = 0
  
  #@print check
  print("(Step 1) - Done making batch-specific meta-info")
  
  #@print check
  print("(Step 1) - Returning the necessary information to the mother function")
  
  #returning back the information
  return(list(affy_filepaths, batch_meta_info))
  
}