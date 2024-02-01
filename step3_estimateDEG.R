#### this script is the third step in the microarray analysis pipeline ####
## this script serves the purpose of taking in the reading the expression matrix, 
## using DREAM function to estimate the DE and constructing a DGE list

## the function will take in three values - the expression matrix, batch meta info
## and the formula to be used for calculating the DE

step3_estimateDEG = function(expr_mtx, batch_meta_info, form){
  
  #@print check
  print("Inside Step 3 - reading expression data and estimating the DEGs .. ")
  
  #parameters for parallel processing
  param = SnowParam(workers = 8, "SOCK", progressbar = TRUE) #variancePartition - default
  
  #### adding an additional step here to reorder the colnames of the expression matrix according to the rownames of the batch_meta_info ####
  expr_mtx = expr_mtx[,match(rownames(batch_meta_info), colnames(expr_mtx))]
  
  #estimate weights using linear mixed model of dream -- Christophe suggests not doing this as this function is mainly for RNA-seq and not microarray
  #dge_voom_weights = voomWithDreamWeights(counts = dge_mtx_norm, formula = form, data = batch_meta_info, BPPARAM = param) #variancePartition -- this gives same values as limma:voom - volcano plot

  #@print check
  print("(Step 3) - Running DREAM")

  #fitting the model using the dream function - using parallel computing
  dge_dream_fit = dream(exprObj = expr_mtx, formula = form, data = batch_meta_info, BPPARAM = param) #variancePartition
  
  #@print check
  print("(Step 3) - Computing statistics with eBayes")

  #computing the t, F and log-odds statistics from a given linear model fit from dream
  dge_dream_fit_ebayes = eBayes(dge_dream_fit) #limma
  
  #extracting the DGEs - two ways of doing this
  #1. extract only the top 'n'genes (ranked on lfc and pval) - setting to 20000 to get all 19000+ genes
  dge_tophits = limma::topTable(dge_dream_fit_ebayes, number = 20000, sort.by = "logFC") #variancePartition - similar to limma - only runs with ebayes
  
  #2. extracting DGEs on the basis of customised pval and lfc cutoffs - if I do this I am left with <50 hits - depending on cutoff
  #dgelist_tophits = topTable(dge_dream_fit_ebayes, p.value = p_val_cutoff, lfc = lfc_cutoff) #making experiments here
  
  #adding back the probe IDs
  dge_tophits$gene_symbol = rownames(dge_tophits)
  
  #@print check
  print("(Step 3) - Accessing the differentially expressed genes")
  
  ###using decideTests to annotate genes to be either up or downreuglated based on the cutoffs that I mention
  #I have to supply the model from eb-fitting
  #setting other criteria - BH - multiple testing correction
  #hard cutoff for p-value and lfc (lfc of 2 is four fold change in expression)
  dge_results = decideTests(dge_dream_fit_ebayes, method = "global", adjust.method = "BH") ## changing method = seperate to method = hierarchical
  
  #extracting the names of genes that are differentially expressed - a dirty trick
  diff_exp_genenames = dge_results[,2] != 0 #here -1 is downregulated, 1 is upregulated and 0 is not significant
  diff_exp_genenames = as.data.frame(diff_exp_genenames)
  diff_exp_genenames = rownames(diff_exp_genenames)[diff_exp_genenames[,1] == TRUE] #since there is only one condition I can select only first column - cautious for multiple conditions!
  
  #subsetting the expression matrix
  diff_exp_genes = subset(expr_mtx, rownames(expr_mtx) %in% diff_exp_genenames)
  
  #giving back sample names - healthy or diseased
  colnames(diff_exp_genes) = batch_meta_info$sample_id
  
  #@print check
  print("(Step 3) - Tabulating the information on the differentially expressed genes")
  
  ## explaination of why I do not use the lfc and p-value cutoff here
  ## https://support.bioconductor.org/p/97271/#97281 - p-value cutoff is at 0.05 by-default, the fold change of 0.5 means a lfc 
  ## of 1.414 (2^0.5) - which is higher hence reduced number of up- or down-regulated genes

  #@print check
  print("(Step 3) - Extracting the top hits from the DGE")
  
  #returning the tophits table back
  return(list(dge_tophits = dge_tophits, diff_exp_genes = diff_exp_genes))

  
}


#### trying code from CC's shared code - 15012024 ####

# #model matrix
# design = model.matrix(~0 + batch_meta_info$group)
# colnames(design) = levels(as.factor(batch_meta_info$group))
# 
# #fitting through linear model
# fit_limma = lmFit(expr_mtx, design)
# 
# #making contrasts
# contrasts_limma = makeContrasts(T2vT0 = T2 - T0, levels = design)
# 
# #making fits to the contrast
# contrast_fit = contrasts.fit(fit_limma, contrasts_limma)
# 
# #making ebayes correction
# dge_dream_fit_ebayes = eBayes(contrast_fit, trend = TRUE, robust = FALSE)

#### resuming the code section ####


#checking the distribution of the expression count per probe (please not2 this increases with increase in sample size - so every experiment has a customised threshold)
#hist(rowSums(expr_mtx)) 

# filter genes by number of counts
#probe_filter = rowSums(expr_mtrx) > 0 #-- optional

#removing negative counts -- this happens when I use raw affy files with combat (without rma)
#expr_mtx[expr_mtx < 0] = 0

#making a DGEList object based on the filter
# dge_mtx = DGEList(expr_mtx) #edgeR
# 
# #@print check
# print("(Step 3) - Normalizing the counts")
# 
# #normalizing the counts
# dge_mtx_norm = calcNormFactors(dge_mtx) #edgeR - post normalization - expression values are still same
# expr_mtx = dge_mtx_norm


#trying a multicore parallel approach
#param = MulticoreParam(workers = 4, progressbar = TRUE, )
#param = SerialParam(progressbar = TRUE)


#fitting the model using the dream function - without parallel computing - 
#dge_dream_fit = dream(exprObj = expr_mtx, formula = form, data = batch_meta_info) #variancePartition

#fitting the model using the dream function
#dge_dream_fit = dream(dge_voom_weights, formula = form, data = batch_meta_info, BPPARAM = param) #variancePartition


#keeping only hits whose p-value is lower than the threshold
#dge_tophits = dge_tophits[dge_tophits$P.Value < p_val_cutoff,] -- trials for cerno plots


#summary of the up and downregulated genes
#summary(dge_results[,2])


#dge_results = decideTests(dge_dream_fit_ebayes, method = "global", adjust.method = "BH",
#                           p.value = 0.05, lfc = 0.05)

# #summary of the upregulated, downregulated and not-significant ones
# summary(dge_results) -- this can be used to get exact count of up and downregulated genes

# #using the data matrix from the differentially expressed genes matrix
# #modifying to keep only the expression values - do not want a tibble
# # I have to select the second column as this first is an intercept
# diff_exp_genes = dge_voom_weights$E[dge_results[,2] !=0, ]

