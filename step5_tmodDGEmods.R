#### this script is the fifth step in the microarray analysis pipeline ####
## this script serves the purpose of identifying the gene modules
## from the tophits from the DGE analysis in step 3

## the function will take in two values - the tophits DGE table
## and the p-value and lfc cutoff (in that order) for heatmap

step5_tmodDGEmods = function(dge_tophits){
  
  ##### formality step - checking and making folder for saving the tmod analysis results #####
  
  #making a folder to save the tmod analysis plots in the backup folder - checking if the folder is already existing
  if(!dir.exists(paste(backup_file, "/tmod_analysis/", sep = ""))){
    
    #creating the directory if it does not exist
    dir.create(paste(backup_file, "/tmod_analysis/", sep = ""))
    
    #print check
    print(paste("Created a new directory to save the tmod analysis results",sep = ""))
    
  }else{print(paste("Found backup directory to save the tmod analysis",sep=""))}
  
  
  #### script section ####
  
  #@print check
  print("Step 5 - Identifying DGE mods using tmod .. ")

  #@print check
  print("(Step 5) - performing cernotest")
  
  #reordering with the adjusted p values
  dge_tophits = dge_tophits[order(dge_tophits$adj.P.Val),]
  
  #manipulating the dataframe to fit for cerno test
  #colnames(dge_tophits)[1] = "gene_symbol"
  rownames(dge_tophits) = NULL
  
  #running the tmodcernotest
  cerno_results = tmodCERNOtest(dge_tophits$gene_symbol, order.by = "pval") ## ----- no data here yet - maybe the log fold change is not drastic enough -- investigate !! -- see line 343 comment - have to reorder p-values
  ## cerno test fails when we use combat to correct for batch effect
  
  # #### brief digression on how I can access the participating genes within each module ####
  # 
  # #extracting participating genes of modules
  # m5_1 = tmod::getGenes(gs = "DC.M5.1")
  # 
  # #finding how many genes are overlapping between two modules
  # overlap_m57_m51 = tmod::modOverlaps(modules = c("DC.M5.1", "DC.M5.7"))
  # 
  # ### digression done ####
    
  #checking if I have more than two rows here - else next steps will fail
  if(dim(cerno_results)[1] >= 2){
  
    #shortlisting the top 30 terms only - according to the adjusted p values
    cerno_results = cerno_results[order(cerno_results$adj.P.Val),]
    cerno_results = cerno_results[1:30,]
    
    #calculating number of significant genes per module - default lfc threshold - 0.5, default pval threshold - 0.05
    #old - sig_genes = tmodDecideTests(g = dge_tophits$gene_symbol, lfc.thr = 0.5, pval.thr = 0.05)
    sig_genes = tmodDecideTests(g = dge_tophits$gene_symbol, lfc = dge_tophits$logFC, pval = dge_tophits$P.Value,
                                lfc.thr = 0.01, pval.thr = 0.05)
    names(sig_genes) = "DGE_mods"
    
    #adding second panel plot for significant genes
    tmod_plots = ggPanelplot(res = list(DGE_mods = cerno_results),sgenes = sig_genes)
    
    #merging the two plots
    tmod_plots = tmod_plots+ggtitle(paste("BTM analysis for ", variable, " between ", meta_names[1], " & ", meta_names[2], sep = ""))+
      theme(plot.title = element_text(size = 18, hjust = 0.5), axis.text.y =element_text(size = 14), axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14),
            legend.title = element_text(size = 14), legend.text = element_text(size = 12))
      
    #saving the plot
    ggsave(plot = tmod_plots, filename = paste(backup_file,"/tmod_analysis/tmod_plots_DEG_modules.png", sep = ""), device = "png", 
           width = 14, height = 9, dpi = 300)
  
  
  }else{print("cerno test failed - less than 2 modules")}
  
}


# #@print check
# print("(Step 5) - Pre-processing to get the gene symbols .. ")

# #extracting the gene symbols of the tophit genes
# regulated_genenames = mapIds(clariomshumantranscriptcluster.db, keys = dge_tophits$probeID,
#                              column = 'SYMBOL', keytype = 'PROBEID')
# 
# #downstream processing the names to make them into a dataframe
# regulated_genenames = as.data.frame(regulated_genenames)
# regulated_genenames$probe_id = rownames(regulated_genenames)
# 
# #manipulating the column and rownames
# rownames(regulated_genenames) = NULL
# colnames(regulated_genenames) = c("gene_symbol", "probe_id")
# 
# #extracting information on the upregulated and downregulated probes just after processing them from topTable - here line 182 (or somewhere around that)
# #regulated_toptable = dgelist_tophits[dgelist_tophits$probeID %in% regulated$probe_id,]
# regulated_genenames = merge(dge_tophits, regulated_genenames, by.x = "probeID", by.y = "probe_id")
# 
# #removing probes for which there are no gene symbols - these could potentially be noncoding elements or lncRNAs
# regulated_genenames = regulated_genenames[!is.na(regulated_genenames$gene_symbol),]
# 
# #ordering the results on the adjusted p values -- without this the test returns no values -- strange!
# regulated_genenames = regulated_genenames[order(regulated_genenames$adj.P.Val),]
# 
# #removing probeID
# regulated_genenames = subset(regulated_genenames, select = -c(probeID))


#making first panel plot
#tmod_plot1 = ggPanelplot(list(DGE_mods = cerno_results))

#adding two panel plots - from tmod tutorial -- optional - to add two lists
#ggPanelplot(list(T5_MiddleAdults_YoungAdults = cerno_results_t5_MA_YA, T5_OldAdults_YoungAdults = cerno_results_t5_OA_YA))