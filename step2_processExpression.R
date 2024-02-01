#### this script is the second step in the microarray analysis pipeline ####
## this script serves the purpose of taking in the affymetrix file paths, read expression data
## normalize it using rma and further correct for batch effects using combat from sva

## the function will take in two values - the filepaths for the affymetrix files
## and the batch-specific metainfo

step2_processExpression = function(affy_filepaths, batch_meta_info){
  
  #@print check
  print("Inside Step 2 - Reading and processing the expression data .. ")
  
  ## writing a short R function that will take as input a single chp file and generate the probe-specific expression values in rows
  read_chp_file = function(filename){
    
    #importing the file - inspired from https://www.reddit.com/r/bioinformatics/comments/aesx07/chp_files/
    dat = readChp(filename = filename)
    
    #extracting the element of interest
    unform_signal = dat[['QuantificationEntries']]
    
    #extracting sample id
    sampleid = str_match(filename, "0000-*(.*?).sst-rma-gene-full.chp")[,2]
    
    #extracting the expression values
    file_expr = data.frame(probename = unform_signal$ProbeSetName, sampleid = unform_signal$QuantificationValue)
    
    #changing the second column name
    colnames(file_expr)[2] = sampleid

    #returning
    return(file_expr)
  }
  
  #### adding a step here - converting the affymetrix ids to their gene IDs - thereby ommitting entries that do not have a gene ID ####
  
  #giving call to the function to extract the expression values per file and storing it is a list
  expr_mtx = lapply(affy_filepaths, function(x){read_chp_file(x)})
  
  #@print check
  print("(Step 2) - Done extracting expression values for the sample sets")
  
  #joining all by their probe names
  expr_mtx = expr_mtx %>% purrr::reduce(left_join, by = "probename")
  
  #adding probename to rownames
  rownames(expr_mtx) = expr_mtx$probename
  
  #removing probename
  expr_mtx = subset(expr_mtx, select = -c(probename))
  
  #changing the column name to also include the timestamp
  colnames(expr_mtx) = unlist(lapply(colnames(expr_mtx), function(x){batch_meta_info$sample_id[batch_meta_info$id == x]}))
  #colnames(expr_mtx) = batch_meta_info$sample_id[batch_meta_info$id == colnames(expr_mtx)]
  
  #extracting the gene symbols
  genename_info = mapIds(clariomshumantranscriptcluster.db, keys = rownames(expr_mtx),
                               column = 'SYMBOL', keytype = 'PROBEID')
  
  #downstream processing the names to make them into a dataframe
  genename_info = as.data.frame(genename_info)
  genename_info$probe_id = rownames(genename_info)
  
  #removing the rows that do not have any gene names
  genename_info = genename_info[!is.na(genename_info$genename_info),]
  
  ## I have to remove any duplicates in the genename column - keeping only a single representative probe ID in case of duplication
  ## probes are proxies for genes, so context dependent upregulation of genes should have upregulation of all duplicate probes (?)
  ## number of genes with duplicate probe IDs - 169/19287 (0.8%)
  genename_info = genename_info[!duplicated(genename_info$genename_info),] #- the first entry is kept per gene name
  
  #removing expression matrix entries now that do not have gene names
  expr_mtx = expr_mtx[rownames(expr_mtx) %in% genename_info$probe_id, ]
  
  #giving back the gene name to the expr mtx
  rownames(expr_mtx) = genename_info$genename_info[genename_info$probe_id == rownames(expr_mtx)]
  
  #have to add the sample_id as rownames --> the column names of the dge_mtx_norm$counts should match the rownames of the batch_meta_info
  rownames(batch_meta_info) = batch_meta_info$sample_id
  
  #### end of geneID conversion ####
  
  #@print check
  print("(Step 2) - Done making expression matrix")
  
  #correcting for batches
  #the batches previously were based on groups that I am studying - T2 v T0 (for example)
  #after correspondence with CC - changing this to the microarray batches (10 in total) - this is based on the design of microarray experiments
  #MJ-17012024
  expr_mtx = ComBat(dat = expr_mtx, batch = batch_meta_info$batch_no)
  
  #@print check
  print("(Step 2) - Done processing the expression matrix through ComBat (sva - correcting batch effect)")
  
  #@print check
  print("(Step 2) - Returning the necessary information to the mother function")

  #returning the expression matrix back that is processed by combat
  return(expr_mtx)
  
}

# ##### section to check the actual change in expression #####
# 
# ## writing a function that can do this for me
# exp_change = function(df){
#   
#   #counting the mean of t0 and t2
#   t0_mean = mean(df[1:150])
#   t2_mean = mean(df[151:300])
#   
#   #taking exp change ratio here - (max/min) ratio
#   change = max(t0_mean, t2_mean)/min(t0_mean, t2_mean)
#   
#   #putting this together in dataframe
#   df_res = data.frame(t0 = t0_mean, t2 = t2_mean, change = change)
#   
#   #returning
#   return(df_res)
#   
# 
# }
# 
# #converting matrix to a dataframe
# expr_mtx_df = as.data.frame(expr_mtx)
# 
# #querying all the rows sequentially
# exp_change_df = apply(expr_mtx_df, 1, function(x){exp_change(x)})
# 
# #binding toegther
# exp_change_df = do.call("rbind", exp_change_df)
# 
# ##### section closed #####

#### correcting for batch effects with ComBat ####
#making a basic model
#mod0 = model.matrix(~1, data = batch_meta_info)

# #accessing the phenoData
# phenoData = pData(affy_read_rma)
# 
# #adding an id column in the pData so that I can link this to the batch meta info
# phenoData$id = str_match(rownames(phenoData), "0000-*(.*?).CEL")[,2]
# 
# #preserving the actual names which are the row names here
# phenoData$names = rownames(phenoData)
# 
# #merging with the batch_meta_info to add the "group" column
# phenoData = merge(phenoData[,c("index", "id", "names")], batch_meta_info[,c("id", "group")], by = "id")
# 
# #giving back the rownames
# rownames(phenoData) = phenoData$names
# 
# #keeping only columns of interest
# phenoData = subset(phenoData, select = c(index, group))

#running combat here- correcting for effects arising in different batchs (ex - T2 vs T0), using empirical bayesian framework,
#make plots to check the fit of the model - prior.plots = TRUE - if the correction fits for expected model
#expr_mtx = ComBat(dat = expr_mtx, mod = mod0, batch = batch_meta_info$group, par.prior = TRUE)


## part where I used cel files and rma normalization ##
# #reading the CEL files using oligo packages
# affy_read = read.celfiles(filenames = affy_filepaths)
# 
# #checking the signal intensities distribution per sample
# #hist(affy_read)
# 
# #@print check
# print("Inside Step 2 - RMA normalization .. ")
# 
# #normalizing the data using RMA method here
# ## Error here when running on Rstudio sever - have to disable multithreading - BiocManager::install("preprocessCore", configure.args="--disable-threading", force = TRUE)
# affy_read_rma = oligo::rma(affy_read) #rma only works with read.celfiles - oligo #I have to add oligo::rma as only rma comes from affy package and that gives errors here
# 
# #checking the sample intensities distribution per sample
# #hist(affy_read_rma)
# 
# #with rma expression check
# expr_mtx = exprs(affy_read_rma)
# 
# #adding probe names
# colnames(expr_mtx) = str_match(colnames(expr_mtx), "0000-*(.*?).CEL")[,2]
# colnames(expr_mtx) = batch_meta_info$sample_id[batch_meta_info$id == colnames(expr_mtx)]
