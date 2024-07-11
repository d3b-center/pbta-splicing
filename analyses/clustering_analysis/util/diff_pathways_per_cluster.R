# function to identify differential pathways using GSVA
suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(GSVA)
  library(pheatmap)
  library(Biobase)
})

# number of pathways per cluster to display
n <- 2

diff_pathways_per_cluster <- function(input_mat, input_clin, cluster_output, n_cluster, gene_set, prefix, output_dir){
  
  # output directory
  dir.create(output_dir, showWarnings = F, recursive = T)
  
  # read data
  input_mat <- readRDS(input_mat)
  input_clin <- readr::read_tsv(input_clin)
  cluster_output <- readRDS(cluster_output)
  gene_set <- readRDS(gene_set)
  
  # cluster tree
  cluster_tree <- cluster_output[[n_cluster]]$consensusTree
  
  # cluster class
  cluster_output <- cluster_output[[n_cluster]]$consensusClass %>%
    as.data.frame() %>%
    rownames_to_column('sample_id') %>%
    dplyr::rename('cluster_class' = '.')
  
  # match up
  expr_mat <- input_mat %>%
    as.data.frame() %>%
    dplyr::select(cluster_output$sample_id)
  cluster_output <- cluster_output %>%
    mutate(tmp = sample_id) %>%
    column_to_rownames('tmp')
  
  # convert into an eset with cluster info
  eset <- Biobase::ExpressionSet(assayData = as.matrix(expr_mat), 
                                 phenoData = new("AnnotatedDataFrame", data = cluster_output))
  
  # GSVA analysis
<<<<<<< HEAD
  gsva_eset <- GSVA::gsva(expr = eset,
                          gset.idx.list = gene_set, 
                          method = "ssgsea", 
                          kcdf = "Gaussian") 
=======
  gsea_scores_param <- gsvaParam(eset,
                                 geneSets = gene_set,
                                 kcdf = "Gaussian",
                                 assay = NA_character_,
                                 annotation = NA_character_,
                                 tau = 1,
                                 minSize = 1, 
                                 maxSize = 1500, ## Arguments from K. Rathi
                                 maxDiff = TRUE) ## Setting this argument to TRUE computes Gaussian-distributed scores (bimodal score distribution if FALSE)
  
  gsva_eset <- gsva(gsea_scores_param, verbose = TRUE)
>>>>>>> main
  
  # convert to long format
  gsva_scores <- gsva_eset@assayData$exprs %>%
    as.data.frame() %>%
    rownames_to_column("geneset") %>%
    gather("sample_id", "score", -c("geneset")) %>%
    dplyr::select(sample_id, geneset, score)
  
  # save gsva output
  readr::write_tsv(x = gsva_scores, file = file.path(output_dir, paste0(prefix, "_gsva_output.tsv")))
  
  # differential expression across clusters
  
  # per cluster analysis
  n_clusters <- unique(gsva_eset[["cluster_class"]])
  all_pathways <- c()
  for(i in 1:length(n_clusters)){
    
    # cluster of interest vs others
    cluster_of_interest <- n_clusters[i]
    type <- ifelse(gsva_eset[["cluster_class"]] == cluster_of_interest, "ref", "others")
    
    # create design and fit
    mod <- model.matrix(~ factor(type))
    colnames(mod) <- c("others", "ref_vs_others")
    rownames(mod) <- gsva_eset[["sample_id"]]
    fit <- limma::lmFit(object = gsva_eset, design = mod)
    fit <- limma::eBayes(fit)
    toptable_output <- limma::topTable(fit, coef = 2, n = Inf)
    toptable_output <- toptable_output %>%
      rownames_to_column("Geneset")
    
    # save filtered output
    out_file <- paste0(prefix, "_cluster_", i, "_pathway.tsv")
    readr::write_tsv(x = toptable_output %>% 
                       filter(P.Value < 0.05), file = file.path(output_dir, out_file))
    
    #  pull N most significant pathways only
    DEpwys <- toptable_output %>%
      filter(adj.P.Val < 0.05) %>%
      arrange(adj.P.Val) %>%
      head(n) %>%
      pull(Geneset) 
    all_pathways <- c(all_pathways, DEpwys)
  }
  
    # plot top 5 pathways per cluster
  all_pathways <- unique(all_pathways)

  DEpwys_es <- Biobase::exprs(gsva_eset[all_pathways, ])
  rownames(DEpwys_es) <- gsub("KEGG_|HALLMARK_", "", rownames(DEpwys_es))
  rownames(DEpwys_es) <- gsub("_", " ", rownames(DEpwys_es))
  DEpwys_annot <- pData(gsva_eset[DEpwys,])
  DEpwys_annot['sample_id'] <- NULL
  DEpwys_annot$cluster_class <- as.character(DEpwys_annot$cluster_class)
  
  DEpwys_annot <- DEpwys_annot %>%
    rownames_to_column("Kids_First_Biospecimen_ID") %>%
    inner_join(input_clin, by = "Kids_First_Biospecimen_ID") %>%
    #inner_join(palette_file, by = "short_histology") %>%
    dplyr::select(Kids_First_Biospecimen_ID, cluster_class, plot_group, plot_group_hex) %>%
    column_to_rownames("Kids_First_Biospecimen_ID")
  
  # rename annotation columns
  DEpwys_annot <- DEpwys_annot %>%
    dplyr::rename("Histology" = "plot_group",
                  "Cluster" = "cluster_class")
  
  # colors for marking different clusters - to match first clustering heatmap
  thisPal <- c("#B2DF8A","#E31A1C","#33A02C","#A6CEE3","#FB9A99","#FDBF6F",
               "#CAB2D6","#FFFF99","#1F78B4","#B15928","#6A3D9A","#FF7F00",
               "#2ef4ca","#f4cced","#bd18ea","#f4cc03","#05188a","#e5a25a",
               "#06f106", #bright green
               "#85848f", #med gray
               "#000000", #black
               "#076f25", #dark green
               "#93cd7f",#lime green
               "#4d0776" #dark purple
  )
  
  l <- thisPal[1:n_cluster]
  names(l) <- as.character(1:n_cluster)
  mycolors <- list()
  mycolors[['Cluster']] <- l 
  
  # create annotation for plot group
  short_histology_palettes <- DEpwys_annot %>%
    dplyr::select(Histology, plot_group_hex) %>%
    unique()
  mycolors[['Histology']] <- short_histology_palettes$plot_group_hex
  names(mycolors[['Histology']]) <- short_histology_palettes$Histology
  
  # remove colors from annotation table
  DEpwys_annot$plot_group_hex <- NULL
  color_palette <- colorRampPalette(c("blue", "white", "darkorange"))
  heat_colors <- color_palette(6)

  pheatmap::pheatmap(DEpwys_es, scale = "row", 
                     fontsize = 8,
                     treeheight_row = 20, 
                     treeheight_col = 20,
                     show_colnames = F, 
                     annotation = DEpwys_annot, 
                     annotation_colors = mycolors, 
                     cluster_cols = cluster_tree, 
                     color = heat_colors,
                     name = "GSVA score",
                     filename = file.path(output_dir, paste0(prefix, "_top", n, "_pathways.pdf")), 
<<<<<<< HEAD
                     width = 12, height = 6)
=======
                     width = 12, height = 5.5)
>>>>>>> main
}
