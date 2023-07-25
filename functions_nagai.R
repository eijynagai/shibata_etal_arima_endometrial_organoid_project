# multipurpose functions

###################
# Quality control
###################
# How many genes have no expression in all the experiment?
count_zero_rows <- function(seurat_object) {
  # Extract the raw count matrix from the list
  my_data <- as.matrix(seurat_object@assays$RNA@counts)

  # Apply the all() function to each row in the matrix to check if all elements are equal to zero
  all_zeros <- apply(my_data, 1, function(row) all(row == 0))

  # Count how many rows have only zeros and return the result
  return(sum(all_zeros))
}
                     
# How many genes have no expression in all the experiment?
count_value_in_rows <- function(seurat_object, value) {
  # Extract the raw count matrix from the list
  my_data <- as.matrix(seurat_object@assays$RNA@counts)

  # Apply the all() function to each row in the matrix to check if all elements are equal to zero
  all_values <- apply(my_data, 1, function(row) all(row == value | row == 0))

  # Count how many rows have only zeros and return the result
  return(sum(all_values))
}
                     

do_QC <- function(seurat_object, outdir, title) {
  suppressPackageStartupMessages({
    library(Seurat)
    library(tidyverse)
    library(patchwork)
    library(sctransform)  
  })
  seurat_object[["percent_mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT[-\\.]")
  seurat_object[["percent_rb"]] <- PercentageFeatureSet(seurat_object, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")

  # Read count distribution report
  #print(seurat_object)
  print(paste0("Number of genes with only zero values: ", count_zero_rows(seurat_object)))
  print(paste0("Number of genes with only 1 or 0 values: ", count_value_in_rows(seurat_object, 1)))
    
    
  plot.width = 16
  plot.height = 4
  options(repr.plot.width = plot.width, repr.plot.height = plot.height)
  plot1 <- VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_rb"), ncol = 4)
  print(plot1)
  ggsave(file=paste0(outdir,"/QC_plots1_",title,".pdf"), width = plot.width, height = plot.height)

  plot.width = 16
  plot.height = 6
  options(repr.plot.width = plot.width, repr.plot.height = plot.height)
  plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent_mt")
  plot3 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(plot2 + plot3)
  ggsave(file=paste0(outdir,"/QC_plots2_",title,".pdf"), width = plot.width, height = plot.height)
    
  plot.width = 16
  plot.height = 6
  options(repr.plot.width = plot.width, repr.plot.height = plot.height)
  plot4 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent_rb")
  plot5 <- FeatureScatter(seurat_object, feature1 = "percent_mt", feature2 = "percent_rb")
  print(plot4 + plot5)
  ggsave(file=paste0(outdir,"/QC_plots3_",title,".pdf"), width = plot.width, height = plot.height)

  # Plots related to cUMI. The distribution of cells with low UMI can also helps defining the threshold 
  plot.width = 16
  plot.height = 6
  options(repr.plot.width = plot.width, repr.plot.height = plot.height)
  seurat_object <- seurat_object %>% AddMetaData(metadata = colSums(seurat_object@assays$RNA@counts), col.name = "nUMI")
  plot6 <- ggplot(seurat_object@meta.data, aes(x = nUMI)) + 
        geom_histogram(binwidth = 100) + 
        labs(x = "UMI Counts", y = "Number of Cells", title = "Histogram of UMI Counts") +
        scale_x_continuous(breaks = seq(min(seurat_object@meta.data$nUMI), max(seurat_object@meta.data$nUMI), by = 1000)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
        #geom_vline(xintercept = threshold_value, color = "red", linetype = "dashed", size = 1)  
  print(plot6)
  ggsave(file=paste0(outdir,"/QC_plots4_hist_",title,".pdf"), width = plot.width, height = plot.height)

  # Density plot
  plot.width = 16
  plot.height = 6
  options(repr.plot.width = plot.width, repr.plot.height = plot.height)
  plot7 <- ggplot(seurat_object@meta.data, aes(x = nUMI)) + 
        geom_density() + 
        labs(x = "UMI Counts", y = "Density", title = "Density Plot of UMI Counts")
  print(plot7)
  ggsave(file=paste0(outdir,"/QC_plots5_density_",title,".pdf"), width = plot.width, height = plot.height)

  # Elbow plot
  plot.width = 16
  plot.height = 6
  options(repr.plot.width = plot.width, repr.plot.height = plot.height)
  umi_counts_sorted <- sort(seurat_object@meta.data$nUMI, decreasing = TRUE)
  cell_rank <- 1:length(umi_counts_sorted)
  elbow_plot_data <- data.frame(Rank = cell_rank, UMI_Counts = umi_counts_sorted)
  plot8 <- ggplot(elbow_plot_data, aes(x = Rank, y = UMI_Counts)) + 
        geom_point() + scale_x_log10() + 
        scale_y_log10() + 
        labs(x = "Cell Rank (log10)", y = "UMI Counts (log10)", title = "Elbow Plot of UMI Counts")
  print(plot8)
  ggsave(file=paste0(outdir,"/QC_plots6_elbow_",title,".pdf"), width = plot.width, height = plot.height)
  
  return(seurat_object)
}



###################
# Analysis
###################

do_SCT_cluster <- function(seurat_object){
  seurat_object <- SCTransform(seurat_object, vars.to.regress = c("percent_mt"), verbose = FALSE)
  seurat_object <- RunPCA(seurat_object, verbose = FALSE)
  seurat_object <- RunUMAP(seurat_object, dims = 1:30, verbose = FALSE)
  seurat_object <- FindNeighbors(seurat_object, dims = 1:30, verbose = FALSE)
  seurat_object <- FindClusters(seurat_object, verbose = FALSE, resolution = c(0.1, 0.3, 0.5, 0.8))
  return(seurat_object)
}
    
    

write2gzip = function(input, outpath) {
    write.table(input, outpath, sep = "\t", row.names = FALSE, quote = FALSE)
    system(paste("gzip --force", outpath))
}


symbol2ensembl = function(dat, keys, from = "SYMBOL", to = "ENSEMBL", db = NULL) {
    if (is.null(db) | !"OrgDb" %in% class(db))
        stop("argument `db` must be an object of class OrgDb")
    if (is.null(from) | !from %in% keytypes(db))
        stop("invalid `from` argument: ", from)
    if (is.null(to) | !to %in% keytypes(db))
        stop("invalid `to` argument: ", to)
    map = suppressMessages(AnnotationDbi::select(db, keys = keys, columns = to,
                                               keytype = from))
    values = map[[to]][match(keys, map[[from]])]
    aggregate(dat, by = list(values), FUN = sum)
}

# Convert the symbols to ENTREZID (necessary for clusterprofiler)
convertID <- function(list_symbols){
    convertedIDs <- bitr(list_symbols,
                         fromType = "SYMBOL",
                         toType = "ENTREZID",
                         OrgDb = "org.Hs.eg.db",
                         drop = TRUE)
    return(convertedIDs$ENTREZID)
}
                      
# GO enrichment analysis 
# ENTREZ gene id is required
go_func <- function(gene_list_entrez, title = "BP: GO enrichment"){

    ego <- enrichGO(gene = gene_list_entrez, 
             OrgDb = org.Hs.eg.db,
             ont ="BP", 
             pAdjustMethod = "BH",
             qvalueCutoff  = 0.01,
             pvalueCutoff = 0.01, 
            )
    dotplot(ego) + 
        theme(text = element_text(size=12)) + 
        ggtitle(title)
}

go_func2 <- function(gene_list_entrez, title = "BP: GO enrichment"){
 
    ego <- enrichGO(gene = gene_list_entrez, 
             OrgDb = org.Hs.eg.db,
             ont ="BP", 
             pAdjustMethod = "BH",
             qvalueCutoff  = 0.01,
             pvalueCutoff = 0.01, 
            )
    mutate(ego, qscore = -log(p.adjust, base=10)) %>% 
            barplot(x = "qscore") + 
            ggtitle(title)
}