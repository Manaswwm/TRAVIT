#### this script is the third step in the microarray analysis pipeline ####
## this script serves the purpose of taking in the reading the expression matrix, 
## using DREAM function to estimate the DE and constructing a DGE list

## the function will take in three values - the expression matrix, batch meta info
## and the formula to be used for calculating the DE

step3_estimateDEG = function(expr_mtrx, batch_meta_info, form){
  
  #@print check
  print("Inside Step 3 - reading expression data and estimating the DEGs .. ")
  
  #checking the distribution of the expression count per probe (please not2 this increases with increase in sample size - so every experiment has a customised threshold)
  #hist(rowSums(expr_mtx)) 
  
  # filter genes by number of counts
  probe_filter = rowSums(expr_mtrx) > 0 #-- optional
  
  #removing negative counts -- this happens when I use raw affy files with combat (without rma)
  #expr_mtx[expr_mtx < 0] = 0
  
  #making a DGEList object based on the filter
  dge_mtx = DGEList(expr_mtrx[probe_filter, ]) #edgeR
  
  #@print check
  print("(Step 3) - Normalizing the counts")
  
  #normalizing the counts
  dge_mtx_norm = calcNormFactors(dge_mtx) #edgeR
  
  #parameters for parallel processing
  param = SnowParam(4, "SOCK", progressbar = TRUE) #variancePartition - default
  
  #have to add the sample_id as rownames --> the column names of the dge_mtx_norm$counts should match the rownames of the batch_meta_info
  rownames(batch_meta_info) = batch_meta_info$sample_id
  
  #estimate weights using linear mixed model of dream
  dge_voom_weights = voomWithDreamWeights(counts = dge_mtx_norm, formula = form, data = batch_meta_info, BPPARAM = param) #variancePartition -- this gives same values as limma:voom - volcano plot
  
  #@print check
  print("(Step 3) - Running DREAM")
  
  #fitting the model using the dream function
  dge_dream_fit = dream(dge_voom_weights, form, batch_meta_info) #variancePartition
  
  #@print check
  print("(Step 3) - Computing statistics with eBayes")
  
  #computing the t, F and log-odds statistics from a given linear model fit from dream
  dge_dream_fit_ebayes = eBayes(dge_dream_fit) #limma
  
  #@print check
  print("(Step 3) - Extracting the top hits from the DGE")
  
  #extracting the DGEs - two ways of doing this
  #1. extract only the top 'n'genes (ranked on lfc and pval)
  dgelist_tophits = topTable(dge_dream_fit_ebayes, number = 10000) #variancePartition - similar to limma
  
  #2. extracting DGEs on the basis of customised pval and lfc cutoffs
  #dgelist_tophits = topTable(dge_dream_fit_ebayes, p.value = 0.5, lfc = 0.01) #making experiments here

  #adding back the probe IDs
  dgelist_tophits$probeID = rownames(dgelist_tophits)
  
  #@print check
  print("(Step 3) - Extracting the top hits from the DGE")
  
  #returning the tophits table back
  return(list(dgelist_tophits, dge_dream_fit_ebayes, dge_voom_weights))

  
}