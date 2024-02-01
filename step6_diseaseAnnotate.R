## using DOSE and enrichplot libraries to create disease annotation plots
step6_diseaseAnnotate = function(dge_tophits){
  
  ##### formality step - checking and making folder for saving the disease ontology results #####
  
  #making a folder to save the disease ontology plots in the backup folder - checking if the folder is already existing
  if(!dir.exists(paste(backup_file, "/disease_annotate/", sep = ""))){
    
    #creating the directory if it does not exist
    dir.create(paste(backup_file, "/disease_annotate/", sep = ""))
    
    #print check
    print(paste("Created a new directory to save the disease annotation results",sep = ""))
    
  }else{print(paste("Found backup directory to save the disease annotation",sep=""))}
  
  
  #### script section ####
  
  #### section 1 - pre-processing the data ####
  #@print check
  print("Step 6 (optional) - performing analysis of enrichment of disease pathways")
  
  #@print check
  print("(Step 6) - Pre-processing step to get the entrez gene ids from gene symbols and filtering")
  
  ##importantly - to make disease-based annotation - ensuring that I have significance levels - p-value and lfc thresholds
  #subsetting the dge_tophits to keep only genes that are of significance
  #caution - I have to take abs so that I take both up and down regulated genes and not only upregulated
  dge_tophits = dge_tophits[dge_tophits$P.Value < 0.05 & abs(dge_tophits$logFC) > 0.3,]
  
  ##the package DOSE and enrichplot heavily use entrez_ids and not gene symbols
  #extracting the gene symbols of the tophit genes
  regulated_entrezids = mapIds(clariomshumantranscriptcluster.db, keys = dge_tophits$gene_symbol,
                               column = 'ENTREZID', keytype = 'SYMBOL')
  
  #downstream processing the names to make them into a dataframe
  regulated_entrezids = as.data.frame(regulated_entrezids)
  regulated_entrezids$gene_symbol = rownames(regulated_entrezids)
  rownames(regulated_entrezids) = NULL
  
  #removing NA - rows with no entrez gene ids
  regulated_entrezids = regulated_entrezids[!is.na(regulated_entrezids$regulated_entrezids),]
  
  #checking if I am at all left with any genes here (atleast 2)
  if(dim(regulated_entrezids)[1] > 1){
  
    #joining to dge_tophits
    dge_tophits = merge(dge_tophits, regulated_entrezids, by = "gene_symbol")
    
    #making a vector of the logFC - entrezIDs as names and logFC as values
    regulated_entrezids_logFC_vector = dge_tophits$logFC
    names(regulated_entrezids_logFC_vector) = dge_tophits$regulated_entrezids
    
    #sorting it on the logFC
    regulated_entrezids_logFC_vector = sort(regulated_entrezids_logFC_vector, decreasing = TRUE)
    
    #### Section 2 - making pathway enrichment analysis ####
    
    #@print check
    print("(Step 6) - Making barplot of the distribution of genes in the disease gene networks")
    
    #first performing enrichment - on human disease gene networks - barplot
    enrich_disgen = enrichDGN(dge_tophits$regulated_entrezids)
    
    #opening the barplot save code
    #jpeg(paste(backup_file, "barplot_disgen_network_barplot.jpg", sep = ""), width = 12, height = 10, units = "in", res = 300)
    
    #making enrichment plot
    barplot = barplot(enrich_disgen)
    
    #saving with ggsave
    ggsave(filename = paste(backup_file, "/disease_annotate/barplot_disgen_network_barplot.jpg", sep = ""), plot = barplot,
           dpi = 300, height = 10, width = 12, device = "jpeg", bg = "white")
    
    #closing the connection
    #dev.off()
    
    #@print check
    print("(Step 6) - Making the network of disease and the associated genes")
    
    #brief-digression - bringing back the gene names from the entrez ids
    enrich_design_genenames = setReadable(enrich_disgen, 'org.Hs.eg.db', 'ENTREZID')
    
    # #opening network save code
    # jpeg(paste(backup_file, "circ_plot_dis_gene_network.jpg", sep = ""), width = 12, height = 10, units = "in", res = 300)
    
    #plotting the disease and gene network
    cnet_plot = cnetplot(enrich_design_genenames, showCategory = 5, color.params = list(foldChange = regulated_entrezids_logFC_vector, edge = TRUE), 
             circular = TRUE, layout = "kk")+ scale_color_gradient2(name='Fold Change', low='blue', high='red')
    
    #saving with ggsave
    ggsave(filename = paste(backup_file, "/disease_annotate/circ_plot_dis_gene_network.jpg", sep = ""), plot = cnet_plot, 
           dpi = 300, height = 10, width = 12, device = "jpeg", bg = "white")
    
    # #closing the connection
    # dev.off()
    
    #@print check
    print("(Step 6) - Making the treeplot to cluster the disease condition on the gsea")
    
    #getting the similarity matrix for the diseases
    enrich_pairwise_disease = pairwise_termsim(enrich_design_genenames)
    
    # #opening tree plot save code
    # jpeg(paste(backup_file, "treeplot_disease_cluster.jpg", sep = ""), width = 12, height = 10, units = "in", res = 300)
    
    #for tree formation - checking if there are atleast 3 terms found
    if(length(enrich_design_genenames$ID) > 2){
      
      #making the treeplot
      treeplot = treeplot(enrich_pairwise_disease, geneClusterPanel = "dotplot", cluster.params = list(method = "average"))
      
      #saving with ggsave
      ggsave(filename = paste(backup_file, "/disease_annotate/treeplot_disease_cluster.jpg", sep = ""), plot = treeplot, 
             dpi = 300, height = 10, width = 12, device = "jpeg", bg = "white")
    
    }
    # #closing the connection
    # dev.off()
    
    #@print check
    print("(Step 6) - Performing gene set enrichment analysis for disease ontology")
    
    #performing the gene set analysis
    gene_set_enrich = gseDO(regulated_entrezids_logFC_vector)
    
    #checking if there are any sets as output - only then proceeding
    if(length(gene_set_enrich$pvalue) > 0){
    
      # #opening the barplot save code
      # jpeg(paste(backup_file, "gse_overview_dotplot.jpg", sep = ""), width = 12, height = 10, units = "in", res = 300)
      
      #plotting the results from the gene set enrichment analysis
      dotplot = dotplot(gene_set_enrich)
      
      #saving with ggsave
      ggsave(filename = paste(backup_file, "/disease_annotate/gse_overview_dotplot.jpg", sep = ""), plot = dotplot, 
             dpi = 300, height = 10, width = 12, device = "jpeg", bg = "white")
      
      # #closing the connection
      # dev.off()
    
    }
  }
}

### for the dotplot I need to process this a bit

#extracting the info on logFC
#regulated_entrezids_logFC = merge(regulated_entrezids, dge_tophits_subset[,c("logFC", "gene_symbol")], by = "gene_symbol")

#deleting probeIDs
#regulated_entrezids_logFC = regulated_entrezids_logFC[,c("regulated_entrezids", "logFC")] #### please check - duplicated entrez ids - multiple probes matching to single entrez ids - how to decide

#removing one of the two or more duplicates for now - FIX THIS
#regulated_entrezids_logFC = regulated_entrezids_logFC[!(duplicated(regulated_entrezids_logFC) | duplicated(regulated_entrezids_logFC, fromLast = TRUE)),] #fall from 

# heatplot(enrich_design_genenames, foldChange =regulated_entrezids_logFC_vector)