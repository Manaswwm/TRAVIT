#this script makes use of the dge tophits that were made in step 3 to make GO and pathway enrichment
step7_pathwayAnnotate = function(dge_tophits){

  
  ##### formality step - checking and making folder for saving the pathway analysis results #####
  
  #making a folder to save the gene ontology plots in the backup folder - checking if the folder is already existing
  if(!dir.exists(paste(backup_file, "/pathway_annotate/", sep = ""))){
    
    #creating the directory if it does not exist
    dir.create(paste(backup_file, "/pathway_annotate/", sep = ""))
    
    #print check
    print(paste("Created a new directory to save the pathway annotation results",sep = ""))
    
  }else{print(paste("Found backup directory to save the pathway annotation",sep=""))}
  

  #### script section ####
  
  #@print check
  print("Step 7 - Performing pathway, biological process and molecular function enrichment analyses ... ")
  
  #subsetting columns in the total dge_tophits
  dge_tophits = dge_tophits[,c("gene_symbol", "logFC", "P.Value")]
  colnames(dge_tophits) = c("Gene.Symbol", "logFC", "P.Value")
  
  ## the analysis here will be on 3 factors - Molecular function, Biological Process and KEGG pathway
  
  #### identification of the biological processes that are enriched within the genes ####
  
  #@ print check
  print("(Step 7) - Performing enrichment analysis of the biological processes")
  
  #running pathfindR 
  #pathfindR automatically detects the number of cores to be used - this might be causing crashes
  #over-riding this to ask pathfindR to use only 4 cores
  dge_tophits_processed_bp = run_pathfindR(dge_tophits, gene_sets = "GO-BP", plot_enrichment_chart = FALSE, n_processes = 6, p_val_threshold = p_val_cutoff) #default is 0.05
  
  ### block of code to check if there are atleast 2 genes in the output - if there is one gene then the heatmap does not work ###
  #collecting all the listed genes
  num_genes = c(dge_tophits_processed_bp$Up_regulated, dge_tophits_processed_bp$Down_regulated)
  
  #removing blank spaces
  num_genes = num_genes[!num_genes == ""]
  
  #checking if the length of unique entries is more than 1
  if(length(unique(num_genes)) > 1){
  
    #@print check
    print("(Step 7) Making heatmap for Biological processes - more than 1 gene detected")
    
    #making heatmap visualization
    go_heatmap_bp = term_gene_heatmap(result_df = dge_tophits_processed_bp, genes_df = dge_tophits, use_description = TRUE, 
                                      legend_title = "Fold Change", sort_terms_by_p = TRUE, num_terms = 3)
    
    #combining with ggplot to add more features in the plot
    go_heatmap_bp = go_heatmap_bp + theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 20), 
                                          legend.title = element_text(size = 18), 
            legend.text = element_text(size = 15), plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) + 
      ggtitle(paste("Biological Process enrichment - ", meta_names[1], " vs ", meta_names[2], "; variable - ", variable, " (p-val = < 0.05)", sep = ""))
      
    #saving the go plot
    ggsave(filename = paste(backup_file, "/pathway_annotate/GO_term_enrichment_heatmap_biological_processes_3.png", sep = ""),
           plot = go_heatmap_bp, device = "png", width = 20, height = 10, dpi = 300)
   
    #clustering the enriched terms
    # dge_tophits_processed_clustered_bp = cluster_enriched_terms(dge_tophits_processed_bp, plot_clusters_graph = FALSE, use_description = TRUE,
    #                                                             plot_dend = FALSE, method = "hierarchical")
    
    #saving the top 10 clusters to a file as a pdf
    kable(dge_tophits_processed_bp[1:10,]) %>% kable_styling(full_width = FALSE, bootstrap_options = "striped") %>%
      save_kable(paste(backup_file, "/pathway_annotate/GO_term_enrichment_table_BP_3.html", sep = ""))
     
  }
  
  #### identification of the molecular functions that are enriched within the genes ####
  
  #@ print check
  print("(Step 7) - Performing enrichment analysis of the molecular functions")
  
  #running pathfindR
  #pathfindR automatically detects the number of cores to be used - this might be causing crashes
  #over-riding this to ask pathfindR to use only 4 cores
  dge_tophits_processed_mf = run_pathfindR(dge_tophits, gene_sets = "GO-MF", plot_enrichment_chart = FALSE, n_processes = 6, p_val_threshold = p_val_cutoff)
  
  ### block of code to check if there are atleast 2 genes in the output - if there is one gene then the heatmap does not work ###
  #collecting all the listed genes
  num_genes = c(dge_tophits_processed_mf$Up_regulated, dge_tophits_processed_mf$Down_regulated)
  
  #removing blank spaces
  num_genes = num_genes[!num_genes == ""]
  
  #checking if the length of unique entries is more than 1
  if(length(unique(num_genes)) > 1){
  
    #@print check
    print("(Step 7) -Making heatmap for Molecular function - more than 1 gene detected")
    
    #making heatmap visualization
    go_heatmap_mf = term_gene_heatmap(result_df = dge_tophits_processed_mf, genes_df = dge_tophits, use_description = TRUE, 
                                      legend_title = "Fold Change", sort_terms_by_p = TRUE, num_terms = 3)
    
    #combining with ggplot to add more features in the plot
    go_heatmap_mf = go_heatmap_mf + theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 20), 
                                          legend.title = element_text(size = 18), 
                                          legend.text = element_text(size = 15), plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) + 
      ggtitle(paste("Molecular Functions enrichment - ", meta_names[1], " vs ", meta_names[2], "; variable - ", variable, " (p-val = < 0.05)", sep = ""))
  
    #saving the go plot
    ggsave(filename = paste(backup_file, "/pathway_annotate/GO_term_enrichment_heatmap_molecular_functions_3.png", sep = ""),
           plot = go_heatmap_mf, device = "png", width = 20, height = 10, dpi = 300)
    
    #clustering the enriched terms
    # dge_tophits_processed_clustered_mf = cluster_enriched_terms(dge_tophits_processed_mf, plot_clusters_graph = FALSE, use_description = TRUE,
    #                                                             plot_dend = FALSE, method = "hierarchical")
    
    #saving the top 10 clusters to a file as a pdf
    kable(dge_tophits_processed_mf[1:10,]) %>% kable_styling(full_width = FALSE, bootstrap_options = "striped") %>%
      save_kable(paste(backup_file, "/pathway_annotate/GO_term_enrichment_table_MF_3.html", sep = ""))
  }
  
  #### identification of the KEGG pathways that are enriched within the genes ####
  
  #@ print check
  print("(Step 7) - Performing enrichment analysis of the KEGG pathways")
  
  #running pathfindR
  #pathfindR automatically detects the number of cores to be used - this might be causing crashes
  #over-riding this to ask pathfindR to use only 4 cores
  dge_tophits_processed_kegg = run_pathfindR(dge_tophits, gene_sets = "KEGG", plot_enrichment_chart = FALSE, n_processes = 6, p_val_threshold = p_val_cutoff)
  
  ### block of code to check if there are atleast 2 genes in the output - if there is one gene then the heatmap does not work ###
  #collecting all the listed genes
  num_genes = c(dge_tophits_processed_kegg$Up_regulated, dge_tophits_processed_kegg$Down_regulated)
  
  #removing blank spaces
  num_genes = num_genes[!num_genes == ""]
  
  #checking if the length of unique entries is more than 1
  if(length(unique(num_genes)) > 1){
    
    #@print check
    print("(Step 7) - Making heatmap for KEGG pathways - more than 1 gene detected")
  
    #making heatmap visualization
    go_heatmap_kegg = term_gene_heatmap(result_df = dge_tophits_processed_kegg, genes_df = dge_tophits, use_description = TRUE, 
                                        legend_title = "Fold Change", sort_terms_by_p = TRUE, num_terms = 3)
    
    #combining with ggplot to add more features in the plot
    go_heatmap_kegg = go_heatmap_kegg + theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 20), 
                                              legend.title = element_text(size = 18), 
                                              legend.text = element_text(size = 15), plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) + 
      ggtitle(paste("KEGG pathway enrichment - ", meta_names[1], " vs ", meta_names[2], "; variable - ", variable, " (p-val = < 0.05)", sep = ""))
    
    #saving the go plot
    ggsave(filename = paste(backup_file, "/pathway_annotate/GO_term_enrichment_heatmap_KEGG_pathways_3.png", sep = ""),
           plot = go_heatmap_kegg, device = "png", width = 20, height = 10, dpi = 300)
  
    #clustering the enriched terms
    # dge_tophits_processed_clustered_kegg = cluster_enriched_terms(dge_tophits_processed_kegg, plot_clusters_graph = FALSE, use_description = TRUE,
    #                                                               plot_dend = FALSE, method = "hierarchical")
    
    #saving the top 10 clusters to a file as a pdf
    kable(dge_tophits_processed_kegg[1:10,]) %>% kable_styling(full_width = FALSE, bootstrap_options = "striped") %>%
      save_kable(paste(backup_file, "/pathway_annotate/GO_term_enrichment_table_KEGG_3.html", sep = ""))
    
  }
  
}