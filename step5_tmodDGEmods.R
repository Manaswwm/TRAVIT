#### this script is the fifth step in the microarray analysis pipeline ####
## this script serves the purpose of identifying the gene modules
## from the tophits from the DGE analysis in step 3

## the function will take in two values - the tophits DGE table
## and the p-value and lfc cutoff (in that order) for heatmap

step5_tmodDGEmods = function(dge_tophits, cutoff){
  
  #@print check
  print("Step 5 - Identifying DGE mods using tmod .. ")
  
  #@print check
  print("(Step 5) - Pre-processing to get the gene symbols .. ")
  
  #extracting the gene symbols of the tophit genes
  regulated_genenames = mapIds(clariomshumantranscriptcluster.db, keys = dge_tophits$probeID,
                               column = 'SYMBOL', keytype = 'PROBEID')
  
  #downstream processing the names to make them into a dataframe
  regulated_genenames = as.data.frame(regulated_genenames)
  regulated_genenames$probe_id = rownames(regulated_genenames)
  
  #manipulating the column and rownames
  rownames(regulated_genenames) = NULL
  colnames(regulated_genenames) = c("gene_symbol", "probe_id")
  
  #extracting information on the upregulated and downregulated probes just after processing them from topTable - here line 182 (or somewhere around that)
  #regulated_toptable = dgelist_tophits[dgelist_tophits$probeID %in% regulated$probe_id,]
  regulated_genenames = merge(dge_tophits, regulated_genenames, by.x = "probeID", by.y = "probe_id")
  
  #removing probes for which there are no gene symbols - these could potentially be noncoding elements or lncRNAs
  regulated_genenames = regulated_genenames[!is.na(regulated_genenames$gene_symbol),]
  
  #ordering the results on the adjusted p values -- without this the test returns no values -- strange!
  regulated_genenames = regulated_genenames[order(regulated_genenames$adj.P.Val),]
  
  #removing probeID
  regulated_genenames = subset(regulated_genenames, select = -c(probeID))
  
  #@print check
  print("(Step 5) - Performing cernotest .. ")
  
  #running the tmodcernotest
  cerno_results = tmodCERNOtest(regulated_genenames$gene_symbol) ## ----- no data here yet - maybe the log fold change is not drastic enough -- investigate !! -- see line 343 comment - have to reorder p-values
  ## cerno test fails when we use combat to correct for batch effect
  
  #checking if I have more than two rows here - else next steps will fail
  if(dim(cerno_results)[1] > 0){
  
    #making first panel plot
    tmod_plot1 = ggPanelplot(list(DGE_mods = cerno_results))
    
    #adding two panel plots - from tmod tutorial -- optional - to add two lists
    #ggPanelplot(list(T5_MiddleAdults_YoungAdults = cerno_results_t5_MA_YA, T5_OldAdults_YoungAdults = cerno_results_t5_OA_YA))
    
    #calculating number of significant genes per module
    sig_genes = tmodDecideTests(g = regulated_genenames$gene_symbol, lfc = regulated_genenames$logFC, 
                                pval = regulated_genenames$adj.P.Val, lfc.thr = cutoff[1], pval.thr = cutoff[2])
    names(sig_genes) = "DGE_mods"
    
    #adding second panel plot for significant genes
    tmod_plot2 = ggPanelplot(list(DGE_mods = cerno_results), sgenes = sig_genes)
    
    #merging the two plots
    plot_grid(tmod_plot1, tmod_plot2)
    
    #saving the plot
    #ggsave(filename = "../../output_plots/comparative_plots/T6_vs_T5_cerno_sigGenes_lfc0.01_pval0.5.png", 
    # device = "png", width = 15, height = 12, dpi = 300) -- optional
  
  
  }else{print("cerno test failed - less than 2 modules")}
  
}