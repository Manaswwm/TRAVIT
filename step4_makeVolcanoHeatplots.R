#### this script is the fourth (optional) step in the microarray analysis pipeline ####
## this script serves the purpose of making the volcano plots and heatmaps 
## to visualize the DGEs

## the function will take in two values - the tophits DGE table
## and the p-value and lfc cutoff (in that order) for heatmap
step4_makeVolcanoHeatplots = function(dge_tophits, batch_meta_info, cutoff, dge_dream_fit_ebayes, dge_voom_weights){
  
  #@print check
  print("Step 4 (optional) - making volcano plots and heatmaps to visualize the DEGs .. ")
  
  
  #@print check
  print("(Step 4) - First making the volcano plot")
  
  #visualizing the dge using a volcano plot
  volcano_plot = ggplot(dge_tophits)+
    geom_point(size = 1, aes(y = -log10(P.Value), x = logFC))+
    #geom_text(aes(y = -log10(P.Value), x = logFC, label = gene_names), data = dgelist_tophits[dgelist_tophits$gene_names %in% chris_ids,])+ -- this is where I can annotate the gene names
    #scale_x_continuous(breaks = seq(-0.2, 0.2, by = 0.1), limits = c(-0.2,0.2))+ -- this is where I can customize the X and Y axis
    ggtitle(label = "Differential gene expression volcano plot")+
    theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
  
  #saving the ggplot to a file
  #ggsave(filename = "../../output_plots/comparative_plots/T2_T0_volcano_plot.jpeg", height = 12, width = 15, dpi = 300, device = "jpeg")
  
  #giving this plot to plotly
  #ggplotly(dfe_volcano) -- for interactive volcano plot visualization
  
  
  #@print check
  print("(Step 4) - Second making the heatmap plots")
  
  #using decideTests to annotate genes to be either up or downreuglated based on the cutoffs that I mention
  #I have to supply the model from eb-fitting
  #setting other criteria - BH - multiple testing correction
  #hard cutoff for p-value and lfc (lfc of 2 is four fold change in expression)
  dge_results = decideTests(dge_dream_fit_ebayes, method = "global", adjust.method = "BH", 
                            p.value = cutoff[1], lfc = cutoff[2])

  # dge_results = decideTests(dge_dream_fit_ebayes, method = "global", adjust.method = "BH",
  #                           p.value = 0.05, lfc = 0.05)


  # #summary of the upregulated, downregulated and not-significant ones
  # summary(dge_results) -- this can be used to get exact count of up and downregulated genes

  #using the data matrix from the differentially expressed genes matrix
  #modifying to keep only the expression values - do not want a tibble
  # I have to select the second column as this first is an intercept
  diff_exp_genes = dge_voom_weights$E[dge_results[,2] !=0, ]

  #giving back sample names - healthy or diseased
  colnames(diff_exp_genes) = batch_meta_info$sample_id
  
  #clustering the rows and columns
  #remember we cluster samples with spearman correlation (captures global patterns and less sensitive to outlier cases)
  #we cluster genes with pearson correlation (captures local patterns and is sensitive to outlier cases)
  #in case of samples - we do not want outliers influencing clustering; in genes - we do want to consider outliers for clustering
  
  #we use hierarchical cluster - unsupervised clustering algorithm
  
  #clustering columns
  #we calculate pairwise distances between genes (logically this is the strength of spearman)
  #calculating distance by 1-cor so that the resulting score is in the scale of 0-2
  clust_cols = hclust(as.dist(1-cor(diff_exp_genes, method = "spearman")), method = "complete")
  
  #clustering rows
  clust_rows = hclust(as.dist(1-cor(t(diff_exp_genes), method = "pearson")), method = "complete")
  
  #cutting the tree - stating how many samples we have in the study - clustering the rows accordingly
  #here k is the number of samples that we have (diseased and healthy)
  module_assign = cutree(clust_rows, k=2)
  
  #assigning colour to each module
  #here too - the colours are influenced by the number of samples (disease vs healthy)
  #getting the colours
  module_colour = rainbow(n = length(unique(module_assign)), start = 0.1, end = 0.9)
  
  #making copies of the colours on the basis of the group assignment
  module_colour = module_colour[as.vector(module_assign)]
  
  #defining a colour pallette
  col_pal = colorRampPalette(colors = c("green", "white", "blue"))(100)
  
  #plotting the heatmap
  #x indicates the exp values of DEGs, Rowv and Colv indicates the clustering of the rows and columns
  #RowSideColours indicates the colours to the row cluster, col indicates the colours for columns
  #by scaling for row I ensure that the comparison is normalised based and not magnitude based
  #for example difference of 100-10 is given the same importance as 10-1
  heatmap_plot = heatmap.2(x = diff_exp_genes, Rowv = as.dendrogram(clust_rows), Colv = as.dendrogram(clust_cols),
            RowSideColors = module_colour, col = col_pal, scale = "row",
            labRow = NA, density.info = "none", trace = "none", cexRow = 1, cexCol = 1, keysize = 1)
  
  #@print check
  print("(Step 4) - Done making the plots - returning them back to the mother function")
  
  #returning the plots back
  return(volcano_plot)
  

}

## code dump
#### @ chekcing probes and genes with highest log fold increase @####
# upreg_genes_info = dgelist_tophits[order(dgelist_tophits$logFC, decreasing = TRUE),]
# upreg_genes_info$logPval = -log10(upreg_genes_info$P.Value)
# upreg_genes_info = upreg_genes_info[upreg_genes_info$logPval>3 & upreg_genes_info$logFC > 0.0,]
# upreg_genes_info$gene_names = mapIds(clariomshumantranscriptcluster.db, keys = upreg_genes_info$probeID, column = "SYMBOL", keytype = "PROBEID")
# upreg_genes_info$gene_names = upreg_genes_names[!is.na(upreg_genes_names)]
# chris_ids = c("GBP2", "UBE2L6", "EPSTI1", "PR2RY14", "FCGR1BP", "GBP5", "GBP1", "FCGR1" , "TRIM22")
# 
# chris_ids[chris_ids %in% upreg_genes_names]

#### check done ####

