#### this script is the fourth (optional) step in the microarray analysis pipeline ####
## this script serves the purpose of making the volcano plots and heatmaps 
## to visualize the DGEs

## the function will take in two values - the tophits DGE table
## and the p-value and lfc cutoff (in that order) for heatmap
step4_makeVolcanoHeatplots = function(dge_tophits, diff_exp_genes, batch_meta_info, expr_mtx){
  
  #@print check
  print("Step 4 - Primary analysis of DGEs through various plots .. ")
  
  
  ##### formality step - checking and making folder for saving the plots made in this function #####
  
  #making a folder to save the gene ontology plots in the backup folder - checking if the folder is already existing
  if(!dir.exists(paste(backup_file, "/dge_plots/", sep = ""))){
    
    #creating the directory if it does not exist
    dir.create(paste(backup_file, "/dge_plots/", sep = ""))
    
    #print check
    print(paste("Created a new directory to save the DGE analysis results",sep = ""))
    
  }else{print(paste("Found backup directory to save the DGE analysis",sep=""))}
  
  #### script section ####
  
  #### @@ Part 1 - raw number of DEGs - lfc-cutoff > 0.1 and p-value cutoff <0.05 @@ ####
  
  #@print check
  print("(Step 4) - making a plot of the number of up and downregulated genes with p-val (0.05) and lfc cutoff (0.1)")
  
  #extracting the number of upregulated genes - pval cutoff - 0.05 and lfc cutoff - 0.1
  upreg_genes_num = dim(dge_tophits[dge_tophits$logFC > 0.001 & dge_tophits$P.Value < 0.05,])[1]
  downreg_genes_num = dim(dge_tophits[dge_tophits$logFC < -0.001 & dge_tophits$P.Value < 0.05,])[1]
  
  #making a dataframe of upreg and downreg
  de_genes = data.frame(genes = c(upreg_genes_num, downreg_genes_num), type = c("Upregulated", "Downregulated"))
  
  #making a barplot capturing the number of DE up and down regulated genes
  de_genes_plot = ggplot(de_genes, aes(y = type, x = as.numeric(genes), fill = type))+
    geom_bar(position = "dodge", stat = "identity", width = 0.4)+
    ylab("Type")+xlab("Number of genes")+geom_text(aes(label = genes), position = position_dodge(width = 0.9), hjust = -0.25, size = 6)+
    scale_fill_manual("Groups", values = c("#0A0DAE", "#F57A00"))+ ggtitle("Number of DEGs (pval < 0.05, lfc > 0.1)")+
    theme_bw()+theme(plot.title = element_text(size = 16, hjust = 0.5), axis.title.y = element_blank(), axis.title.x = element_text(size = 16),
                     legend.title = element_text(size = 16), legend.text = element_text(size = 14),
                     axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_text(size = 14))
  
  #saving the pca plot
  ggsave(plot = de_genes_plot, filename = paste(backup_file, "/dge_plots/de_genes.jpeg", sep = ""),
         dpi = 300, height = 6, width = 15, device = "jpeg", bg = "white")
  
  
  #### @@ Part 2 - clustering of samples through PCA @@ ####
  
  #@print check
  print("(Step 4) - performing clustering of samples using PCA")
  
  #establishing the pca components first
  pc_comp = prcomp(t(expr_mtx)) #using transpose so that clustering is on samples and not genes
  
  #extracting the eigen values per pc
  pc_var = pc_comp$sdev^2
  
  #calculating the percentage variance explained by each component
  pc_per = round(pc_var/sum(pc_comp$sdev^2)*100, 1)
  
  #making a tibble with the components
  pc_comp_df = as_tibble(pc_comp$x)
  
  #extracting the groups from batch_meta_info
  groups = batch_meta_info$group

  #making pca plot
  pca_plot = ggplot(pc_comp_df, aes(x=PC1, y=PC2, colour = groups))+
  geom_point(size=4, alpha = 0.7)+stat_ellipse(alpha = 0.5)+ggtitle(paste("PCA for ", variable, " between ", meta_names[1], " & ", meta_names[2], sep = ""))+
  xlab(paste0("PC1 (",pc_per[1],"%",")"))+ylab(paste0("PC2 (",pc_per[2],"%",")"))+
  scale_colour_manual("Groups", values = c("#EC0B43", "#58355E"))+
  theme_bw()+theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
                   legend.title = element_text(size = 16), legend.text = element_text(size = 14),
                   plot.title = element_text(size = 18, hjust = 0.5))
  
  #saving the pca plot
  ggsave(plot = pca_plot, filename = paste(backup_file, "/dge_plots/pca_plot.jpeg", sep = ""),
         dpi = 300, height = 10, width = 12, device = "jpeg", bg = "white")
   
  
  #### @@ Part 3 - clustering samples through UMAP @ ####
  
  #@print check
  print("(Step 4) - performing clustering of samples using UMAP")
  
  #performing umap analysis on the expr_mtx directly (similar to PCA)
  umap_comp = umap::umap(t(expr_mtx))
  
  #extracting the coordinates of the samples in the umap
  umap_layout = umap_comp[["layout"]]
  umap_layout = data.frame(umap_layout)
  
  #adding information on the sampleIDs as a column
  umap_layout$sample_id = rownames(umap_layout)
  
  #the sampleIDs are stored in the rownames - using this to add information on the groups
  umap_layout = merge(umap_layout, batch_meta_info[,c("sample_id", "group")], by = "sample_id")
  
  #removing sampleIDs now
  umap_layout = subset(umap_layout, select = -c(sample_id))
  
  #plotting the umap now
  umap_plot = ggplot(umap_layout, aes(x = X1, y = X2, color = group)) + geom_point(size = 4, alpha = 0.7)+stat_ellipse(alpha = 0.5)+
    xlab("X1")+ylab("X2")+ggtitle(paste("UMAP for ", variable, " between ", meta_names[1], " & ", meta_names[2], sep = ""))+
    scale_colour_manual("Groups", values = c("#EC0B43", "#58355E"))+
    theme_bw()+theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
                     legend.title = element_text(size = 16), legend.text = element_text(size = 14),
                     plot.title = element_text(size = 18, hjust = 0.5))
  
  #saving the plot
  ggsave(plot = umap_plot, filename = paste(backup_file, "/dge_plots/umap_plot.jpeg", sep = ""),
         dpi = 300, height = 10, width = 12, device = "jpeg", bg = "white")

    
  #### @@ Part 4 - clustering samples through UMAP @ ####
  
  #@print check
  print("(Step 4) - performing volcano plot to visualize the DEGs")
  
  #filtering out genes that are not significant hits
  dge_tophits_significant = dge_tophits[dge_tophits$P.Value < 0.05,]
  
  #extracting upreg and downreg info
  upreg_genes = dge_tophits_significant[dge_tophits_significant$logFC > 0,]
  downreg_genes = dge_tophits_significant[dge_tophits_significant$logFC < 0,]
  
  #sorting and keeping the top 5 entries only
  upreg_genes = upreg_genes[order(upreg_genes$logFC, decreasing = TRUE),][1:5,]
  downreg_genes = downreg_genes[order(downreg_genes$logFC, decreasing = FALSE),][1:5,]
  
  #adding labels
  upreg_genes$trend = rep("Up-regulated")
  downreg_genes$trend = rep("Down-regulated")
  
  #merging
  dge_info = rbind(upreg_genes, downreg_genes)
  
  #visualizing the dge using a volcano plot
  volcano_plot = ggplot(dge_tophits, aes(y = -log10(P.Value), x = logFC, colour = logFC > 0))+
    geom_point(size = 1, alpha = 0.6)+ geom_hline(yintercept = 1.30, linetype = "dashed", color = "black", linewidth = 1)+
    scale_color_manual(values = setNames(c("#F57A00", "#0A0DAE"), c(TRUE, FALSE)))+
    ggtitle(paste("Volcano plot for ", variable, " between ", meta_names[1], " & ", meta_names[2], sep = ""))+ theme_cowplot()+ylab("-log10(p-value)")+xlab("log(Fold Change)")+
    theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), legend.position = "none", 
          plot.title = element_text(size = 20, hjust = 0.5), axis.text.x =element_text(size = 16), axis.text.y = element_text(size = 16))
  
  #adding labels of top 5 up and downregulated genes
  volcano_plot = volcano_plot + geom_label_repel(data = dge_info, 
                                  mapping = aes(x = logFC, y = -log10(P.Value), label = gene_symbol, color = "black"), 
                                  size = 8)
  
  
  
  #saving the ggplot to a file - I have to add a white backgroun - strange
  ggsave(plot = volcano_plot, filename = paste(backup_file, "/dge_plots/volcano_plot_DEG_labels.jpeg", sep = ""),
         dpi = 300, height = 10, width = 12, device = "jpeg", bg = "white")
  
  
  #### @@ Part 4 - plotting DEGs through heatmap @ ####
  
  #@print check
  print("(Step 4) - performing heatmap plot to visualize the DEGs")
  
  ##upregulated genes - logfc is positive
  upreg_geneinfo = dge_tophits[dge_tophits$logFC > 0 & dge_tophits$P.Value < 0.05,] #setting p-val cutoff

  #ordering on increase in logfc
  upreg_geneinfo = upreg_geneinfo[order(upreg_geneinfo$logFC, decreasing = TRUE),]

  #taking the top 20 genes
  upreg_geneinfo = upreg_geneinfo$gene_symbol[1:5]

  ##downregulated genes - logfc is positive
  downreg_geneinfo = dge_tophits[dge_tophits$logFC < 0 & dge_tophits$P.Value < 0.05,]  #setting p-val cutoff

  #ordering on increase in logfc
  downreg_geneinfo = downreg_geneinfo[order(downreg_geneinfo$logFC, decreasing = FALSE),]

  #taking the top 20 genes
  downreg_geneinfo = downreg_geneinfo$gene_symbol[1:5]
  
  #subsetting the diff_exp_genes expression matrix to contain only the top 10 up and down regulated genes
  diff_exp_genes_subset = expr_mtx[rownames(expr_mtx) %in% c(upreg_geneinfo, downreg_geneinfo),]
  #diff_exp_genes_subset = diff_exp_genes[rownames(diff_exp_genes) %in% c(upreg_geneinfo, downreg_geneinfo),]
  
  #clustering columns
  #we calculate pairwise distances between genes (logically this is the strength of spearman)
  #calculating distance by 1-cor so that the resulting score is in the scale of 0-2
  clust_cols = hclust(as.dist(1-cor(diff_exp_genes_subset, method = "spearman")), method = "complete")
  
  #clustering rows
  clust_rows = hclust(as.dist(1-cor(t(diff_exp_genes_subset), method = "pearson")), method = "complete")
  
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
  col_pal = colorRampPalette(colors = c("red", "white", "green"))(10)
  
  #making row annotation
  my_gene_col = data.frame(cluster = ifelse(test = module_assign == 1, yes = "upregulated", no = "downregulated"))
  
  # #making column annotation - messy section but automated
  # #in short - the last three columns of the batch meta info are always - sample_id, group 1 and group 2 - using advantage of this structure
  # my_sample_col = batch_meta_info[,colnames(batch_meta_info[tail(names(batch_meta_info), 5)])]
  # my_sample_col$group[my_sample_col[,colnames(batch_meta_info[tail(names(batch_meta_info), 5)])[4]] == 1] = colnames(batch_meta_info[tail(names(batch_meta_info), 5)])[4]
  # my_sample_col$group[my_sample_col[,colnames(batch_meta_info[tail(names(batch_meta_info), 5)])[5]] == 1] = colnames(batch_meta_info[tail(names(batch_meta_info), 5)])[5]
  # rownames(my_sample_col) = my_sample_col$sample_id
  # my_sample_col = subset(my_sample_col, select = c(sample_id, group))
  
  #making a dataframe giving metainformation on the samples
  my_sample_col = batch_meta_info[,c("id", "group", "sample_id")]
  
  #renaming the column name
  colnames(my_sample_col) = c("subj_id", "group", "sample_id")
  
  #adding additional columns here - still under construction
  #gender distribution of the samples - extracting the gender from the metainfo dataframe by using expression matching
  my_sample_col$gender = unlist(lapply(sub("\\_.*", "", my_sample_col$subj_id), 
                                   function(x) {metainfo$gender[grep(x = metainfo$microarray_id, pattern = x)]}))
  
  #removing sample ID column
  my_sample_col = subset(my_sample_col, select = c(group, gender))

  #not clustering columns
  heatmap = pheatmap(mat = diff_exp_genes_subset, legend = TRUE,
           RowSideColors = module_colour, col = col_pal, scale = "row", cluster_cols = FALSE, cluster_rows = TRUE,
           labRow = NA, density.info = "none", trace = "none", cexRow = 1, keysize = 2, 
           annotation_row = my_gene_col, annotation_col = my_sample_col, labels_col = rep("", dim(my_sample_col)[1]), cutree_rows = 2)

  #saving with ggsave
  ggsave(filename = paste(backup_file, "/dge_plots/heatmap_plot_DEG.jpg", sep = ""), plot = heatmap,
         dpi = 300, height = 10, width = 12, device = "jpeg", bg = "white")
  
  #closing the dev
  dev.off()
  
  #@print check
  print("(Step 4) - Done making the plots - saved in the backup folder")
  
  #returning the plots back
  #return(volcano_plot)
  

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

#we use hierarchical cluster - unsupervised clustering algorithm

# #clustering columns
# #we calculate pairwise distances between genes (logically this is the strength of spearman)
# #calculating distance by 1-cor so that the resulting score is in the scale of 0-2
# clust_cols = hclust(as.dist(1-cor(diff_exp_genes, method = "spearman")), method = "complete")
# 
# #clustering rows
# clust_rows = hclust(as.dist(1-cor(t(diff_exp_genes), method = "pearson")), method = "complete")
# 
# #cutting the tree - stating how many samples we have in the study - clustering the rows accordingly
# #here k is the number of samples that we have (diseased and healthy)
# module_assign = cutree(clust_rows, k=2)
# 
# #assigning colour to each module
# #here too - the colours are influenced by the number of samples (disease vs healthy)
# #getting the colours
# module_colour = rainbow(n = length(unique(module_assign)), start = 0.1, end = 0.9)
# 
# #making copies of the colours on the basis of the group assignment
# module_colour = module_colour[as.vector(module_assign)]
# 
# #defining a colour pallette
# col_pal = colorRampPalette(colors = c("red", "white", "green"))(100)
# 
# #opening the heatmap save code
# png(paste(backup_file, "heatmap_plot_DEG.png"), width = 12, height = 10, units = "in", res = 300)
# 
# #plotting the heatmap
# #x indicates the exp values of DEGs, Rowv and Colv indicates the clustering of the rows and columns
# #RowSideColours indicates the colours to the row cluster, col indicates the colours for columns
# #by scaling for row I ensure that the comparison is normalised based and not magnitude based
# #for example difference of 100-10 is given the same importance as 10-1
# heatmap.2(x = diff_exp_genes, Rowv = as.dendrogram(clust_rows), Colv = as.dendrogram(clust_cols),
#           RowSideColors = module_colour, col = bluered(80), scale = "row",
#           labRow = NA, density.info = "none", trace = "none", cexRow = 1, cexCol = 1, keysize = 2, 
#           dendrogram = "column")

# #adding labels - T2 vs T0
# dge_tophits$genelabels = ifelse(dge_tophits$gene_symbol == "FCGR1A" |
#                                 dge_tophits$gene_symbol == "GBP5"|
#                                 dge_tophits$gene_symbol == "GBP1"|
#                                 dge_tophits$gene_symbol == "FCGR1BP"|
#                                 dge_tophits$gene_symbol == "GBP2"|
#                                 dge_tophits$gene_symbol == "UBE2L6"|
#                                 dge_tophits$gene_symbol == "EPSTI1"|
#                                 dge_tophits$gene_symbol == "STAT1"|
#                                 dge_tophits$gene_symbol == "TRIM22"|
#                                 dge_tophits$gene_symbol == "P2RY14", TRUE, FALSE)

#### extracting labels of top 5 up and downregulated genes ####

#checking if dev is off
#dev.off()

# #opening the heatmap save code
# png(paste(backup_file, "heatmap_plot_DEG.png", sep = ""), width = 12, height = 10, units = "in", res = 300)

#plotting the heatmap
#x indicates the exp values of DEGs, Rowv and Colv indicates the clustering of the rows and columns
#RowSideColours indicates the colours to the row cluster, col indicates the colours for columns
#by scaling for row I ensure that the comparison is normalised based and not magnitude based
#for example difference of 100-10 is given the same importance as 10-1
# pheatmap(mat = diff_exp_genes_subset, Rowv = as.dendrogram(clust_rows), Colv = as.dendrogram(clust_cols),
#           RowSideColors = module_colour, col = col_pal, scale = "row",
#           labRow = NA, density.info = "none", trace = "none", cexRow = 1, cexCol = 1, keysize = 2, 
#           dendrogram = "column", annotation_row = my_gene_col, annotation_col = my_sample_col, labels_col = rep("", dim(my_sample_col)[1]))

#clustering columns
# pheatmap(mat = diff_exp_genes_subset, legend = TRUE,
#          RowSideColors = module_colour, col = col_pal, scale = "row", cluster_cols = TRUE, cluster_rows = TRUE,
#          labRow = NA, density.info = "none", trace = "none", cexRow = 1, keysize = 2, 
#         annotation_row = my_gene_col, annotation_col = my_sample_col, labels_col = rep("", dim(my_sample_col)[1]), cutree_rows = 2)
# 


# #vital age distribution of the sample - extracting VAG from the metainfor dataframe by using expression matching
# my_sample_col$timepoint = unlist(lapply(sub("\\_.*", "", my_sample_col$sample_id), 
#                                      function(x) {metainfo$timepoint[grep(x = metainfo$microarray_id, pattern = x)]}))