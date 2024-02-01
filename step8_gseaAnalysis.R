#### this script is the eigth step in the microarray analysis pipeline 
## this would be an additional/alternative to the already used pathways analysis using path_findR
## the tutorial is accessed here - https://yulab-smu.top/biomedical-knowledge-mining-book/universal-api.html

## the function will take in information on the dge_tophits
step8_gseaAnalysis = function(dge_tophits){

  #@print check
  print("Step 8 - Performing gene set enrichment analysis .. ")
  
  ##### formality step - checking and making folder for saving the plots made in this function #####
  
  #making a folder to save the gene ontology plots in the backup folder - checking if the folder is already existing
  if(!dir.exists(paste(backup_file, "/GSEA/", sep = ""))){
    
    #creating the directory if it does not exist
    dir.create(paste(backup_file, "/GSEA/", sep = ""))
    
    #print check
    print(paste("Created a new directory to save the GSEA analysis results",sep = ""))
    
  }else{print(paste("Found backup directory to save the GSEA analysis",sep=""))}
  
  #### script section ####
  
  #metainformation on all gene sets included in humans
  #hsap_msig_meta = msigdbr(species = "Homo sapiens") - loading the general msigdb for hsap
  
  #### accessing the different GSEA modules and modifying their dataframes ####
  
  #@print check
  print("(Step 8) - preaparing the gsea modules")
  
  #category C7 corresponds to immunological signature gene sets - subetting to keep columns of interest
  c7_sign = msigdbr(species = "Homo sapiens", category = "C7")
  
  #making the dataframe to shortlist the info relevant for gsea
  c7_sign_subset = c7_sign[,c("gs_name", "gene_symbol")]
  
  #category c7 and subcategory vax corresponds to immunological gene sets that are specific for vaccination
  c7_vax_sign = msigdbr(species = "Homo sapiens", category = "C7", subcategory = "VAX")
  
  #making the dataframe to shortlist the info relevant for gsea
  c7_vax_sign_subset = c7_vax_sign[,c("gs_name", "gene_symbol")]
  
  #category hallmark corresponds to hallmark gene sets - well-defined biological sets or processes
  hallmark_sign = msigdbr(species = "Homo sapiens", category = "H")
  
  #making the dataframe to shortlist the info relevant for gsea
  hallmark_sign_subset = hallmark_sign[,c("gs_name", "gene_symbol")]
  
  
  #### accessing the dge_tophits information and formulating the data for supplying it to GSEA ####

  #accessing the upregulated genes
  upreg_genes = dge_tophits$gene_symbol[dge_tophits$logFC > 0 & dge_tophits$P.Value < 0.05]
  downreg_genes = dge_tophits$gene_symbol[dge_tophits$logFC < 0 & dge_tophits$P.Value < 0.05]
  
  #also adding lof fold change
  de_genes_df = data.frame(genes = c(upreg_genes,downreg_genes), lfc = dge_tophits$logFC[dge_tophits$gene_symbol %in% c(upreg_genes, downreg_genes)])
  
  #converting this into a ranked numeric vector
  de_genes_info = de_genes_df$lfc
  names(de_genes_info) = de_genes_df$genes
  
  #sorting
  de_genes_info = sort(de_genes_info, decreasing = TRUE)
  
  #### making the analysis plots here - for each of the three modules I will make a dotplot, heatmap plot and table ####
  
  ### Part 1 - general C7 ###
  
  #@print check
  print("(Step 8) - performing gsea for c7 (general immunology module)")
  
  #making GSEA analysis
  c7_gsea = GSEA(de_genes_info, TERM2GENE = c7_sign_subset, pvalueCutoff = 0.05)
  
  #extracting the output table of results
  c7_gsea_result = c7_gsea@result
  
  #checking if the dimension of this dataframe is more than 2
  if(dim(c7_gsea_result)[1] > 2){
    
    #@print check
    print("(Step 8) - for c7 - found more than 2 gene sets hence making next analysis")
  
    #plotting the gsea results
    c7_gsea_plot = enrichplot::dotplot(object = c7_gsea, showCategory = 5)
    
    #beautifying the plot with ggplot
    c7_gsea_plot = c7_gsea_plot+ggtitle(paste("GSEA for Cluster 7: ", variable, " between ", meta_names[1], " & ", meta_names[2], sep = ""))+
      theme_bw()+theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
                       axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
                       legend.title = element_text(size = 16), legend.text = element_text(size = 14),
                       plot.title = element_text(size = 18, hjust = 0.5))
    
    #to the table, adding the description of the gene set term
    c7_gsea_result$description_entire = unlist(lapply(c7_gsea_result$ID, function(x){unique(c7_sign$gs_description[c7_sign$gs_name == x])}))
    
    #saving the plot
    ggsave(plot = c7_gsea_plot, filename = paste(backup_file, "/GSEA/c7_gsea_plot.jpeg", sep = ""),
           dpi = 300, height = 10, width = 10, device = "jpeg", bg = "white")
    
    #saving the tabulated results
    write.table(x = c7_gsea_result, file = paste(backup_file, "/GSEA/c7_gsea_table.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  }
  
  ### Part 2 - general C7 - VAX ###
  
  #@print check
  print("(Step 8) - performing gsea for c7 subcategory VAX")
  
  #making GSEA analysis
  c7_vax_gsea = GSEA(de_genes_info, TERM2GENE = c7_vax_sign_subset, pvalueCutoff = 0.05)
  
  #extracting the output table of results
  c7_vax_gsea_result = c7_vax_gsea@result
  
  #checking if the dimension of this dataframe is more than 2
  if(dim(c7_vax_gsea_result)[1] > 2){
    
    #@print check
    print("(Step 8) - for c7 VAX - found more than 2 gene sets hence making next analysis")
    
    #plotting the gsea results
    c7_vax_gsea_plot = enrichplot::dotplot(object = c7_vax_gsea, showCategory = 5)
    
    #beautifying the plot with ggplot
    c7_vax_gsea_plot = c7_vax_gsea_plot+ggtitle(paste("GSEA for Cluster 7 (VAX): ", variable, " between ", meta_names[1], " & ", meta_names[2], sep = ""))+
      theme_bw()+theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
                       axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
                       legend.title = element_text(size = 16), legend.text = element_text(size = 14),
                       plot.title = element_text(size = 18, hjust = 0.5))
    
    #to the table, adding the description of the gene set term
    c7_vax_gsea_result$description_entire = unlist(lapply(c7_vax_gsea_result$ID, function(x){unique(c7_vax_sign$gs_description[c7_vax_sign$gs_name == x])}))
    
    #saving the plot
    ggsave(plot = c7_vax_gsea_plot, filename = paste(backup_file, "/GSEA/c7_vax_gsea_plot.jpeg", sep = ""),
           dpi = 300, height = 12, width = 10, device = "jpeg", bg = "white")
    
    #saving the tabulated results
    write.table(x = c7_vax_gsea_result, file = paste(backup_file, "/GSEA/c7_vax_gsea_table.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  }
  
  ### Part 3 - Hallmark ###
  
  #@print check
  print("(Step 8) - performing gsea for Hallmark genes")
  
  #making GSEA analysis
  hallmark_gsea = GSEA(de_genes_info, TERM2GENE = hallmark_sign_subset, pvalueCutoff = 0.05)
  
  #extracting the output table of results
  hallmark_gsea_result = hallmark_gsea@result
  
  #checking if the dimension of this dataframe is more than 2
  if(dim(hallmark_gsea_result)[1] > 2){
  
    #@print check
    print("(Step 8) - for hallmark - found more than 2 gene sets hence making next analysis")
    
    #plotting the gsea results
    hallmark_gsea_plot = enrichplot::dotplot(object = hallmark_gsea, showCategory = 5)
    
    #beautifying the plot with ggplot
    hallmark_gsea_plot = hallmark_gsea_plot+ggtitle(paste("Hallmark: ", variable, " between ", meta_names[1], " & ", meta_names[2], sep = ""))+
      theme_bw()+theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
                       axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
                       legend.title = element_text(size = 16), legend.text = element_text(size = 14),
                       plot.title = element_text(size = 18, hjust = 0.5))
    
    #to the table, adding the description of the gene set term
    hallmark_gsea_result$description_entire = unlist(lapply(hallmark_gsea_result$ID, function(x){unique(hallmark_sign$gs_description[hallmark_sign$gs_name == x])}))
    
    #saving the plot
    ggsave(plot = hallmark_gsea_plot, filename = paste(backup_file, "/GSEA/hallmark_gsea_plot.jpeg", sep = ""),
           dpi = 300, height = 12, width = 10, device = "jpeg", bg = "white")
    
    #saving the tabulated results
    write.table(x = c7_vax_gsea_result, file = paste(backup_file, "/GSEA/hallmark_gsea_table.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  }
}
##### enricher plot #####

# ## reason why I only do GSEA and prefer this over enricher - https://www.biostars.org/p/9468150/#:~:text=enrichr()%20performs%20a%20hypergeometric%20test,list%20than%20expected%20by%20chance.
# 
# #using enricher to perform overrepresentation analysis
# vax_enricher = enricher(upreg_genes, TERM2GENE = c7_vax_sign_subset, pvalueCutoff = 0.05, pAdjustMethod = "BH")
# 
# #making enricher plots
# enrichplot::heatplot(vax_enricher, showCategory = 5