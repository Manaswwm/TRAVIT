#### this scirpt will directly extract the expression values from the expresison set provided by CC
## the expression values here are pre-normalized - I just have to extract them as per the requirement

## the function will take in two values - the filepaths for the affymetrix files
## and the batch-specific metainfo
step_2_1_processExpressionC = function(affy_filepaths, batch_meta_info){
  
  #first, importing the expression set from CC
  #loading the expression set that has been batch corrected
  eset = readRDS("../../../microarray_data/data_from_Christophe/eset.rds")
  
  #extracting the expression values
  expr_mtx_cc = exprs(eset)
  
  #extracting the columns of interest using the sample_id column from batch_meta_info
  expr_mtx_cc = expr_mtx_cc[,as.character(batch_meta_info$sample_id)]
  
  #extracting the gene symbols
  genename_info = mapIds(clariomshumantranscriptcluster.db, keys = rownames(expr_mtx_cc),
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
  expr_mtx_cc = expr_mtx_cc[rownames(expr_mtx_cc) %in% genename_info$probe_id, ]
  
  #giving back the gene name to the expr mtx
  rownames(expr_mtx_cc) = genename_info$genename_info[genename_info$probe_id == rownames(expr_mtx_cc)]
  
  #checking dim after the gene name additiong
  #dim(expr_mtx_cc) #19287X297
  
  #removing the expression values for the two individuals that do not have a responder data - MJ17012024
  #expr_mtx_cc = expr_mtx_cc[,!colnames(expr_mtx_cc) %in% c("118_A", "118_B", "281_A", "281_B")]
  
  #return
  return(expr_mtx = expr_mtx_cc)
  
  
      
}