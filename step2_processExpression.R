#### this script is the second step in the microarray analysis pipeline ####
## this script serves the purpose of taking in the affymetrix file paths, read expression data
## normalize it using rma and further correct for batch effects using combat from sva

## the function will take in two values - the filepaths for the affymetrix files
## and the batch-specific metainfo

step2_processExpression = function(affy_filepaths, batch_meta_info){
  
  
  #@print check
  print("Inside Step 2 - Reading and processing the expression data .. ")
  
  #reading the CEL files using oligo packages
  affy_read = read.celfiles(filenames = affy_filepaths)
  
  #@print check
  print("Inside Step 2 - RMA normalization .. ")
  
  #normalizing the data using RMA method here
  affy_read_rma = rma(affy_read) #rma only works with read.celfiles - oligo
  
  #with rma expression check
  expr_mtx = exprs(affy_read_rma)
  
  #adding probe names
  colnames(expr_mtx) = str_match(colnames(expr_mtx), "0000-*(.*?).CEL")[,2]
  colnames(expr_mtx) = batch_meta_info$sample_id[batch_meta_info$id == colnames(expr_mtx)]
  
  #@print check
  print("(Step 2) - Done reading files and making expression matrix")
  
  #### correcting for batch effects with ComBat ####
  #making a basic model
  mod0 = model.matrix(~1, data = batch_meta_info)
  
  #accessing the phenoData
  phenoData = pData(affy_read_rma)
  
  #adding an id column in the pData so that I can link this to the batch meta info
  phenoData$id = str_match(rownames(phenoData), "0000-*(.*?).CEL")[,2]
  
  #preserving the actual names which are the row names here
  phenoData$names = rownames(phenoData)
  
  #merging with the batch_meta_info to add the "group" column
  phenoData = merge(phenoData[,c("index", "id", "names")], batch_meta_info[,c("id", "group")], by = "id")
  
  #giving back the rownames
  rownames(phenoData) = phenoData$names
  
  #keeping only columns of interest
  phenoData = subset(phenoData, select = c(index, group))
  
  #running combat here
  expr_mtx = ComBat(dat = expr_mtx, mod = mod0, batch = phenoData$group, par.prior = TRUE)
  
  #@print check
  print("(Step 2) - Done processing the expression matrix through ComBat (sva - correcting batch effect)")
  
  #@print check
  print("(Step 2) - Returning the necessary information to the mother function")
  
  #returning the expression matrix back that is processed by combat
  return(expr_mtx)
  
}