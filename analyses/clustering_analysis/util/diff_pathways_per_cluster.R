# function to identify differential pathways using GSVA
suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(GSVA)
  library(pheatmap)
  library(Biobase)
})

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
  gsva_eset <- GSVA::gsva(expr = eset,
                          gset.idx.list = gene_set, 
                          method = "ssgsea", 
                          kcdf = "Gaussian") 
  
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
    readr::write_tsv(x = toptable_output %>% filter(P.Value < 0.05), file = file.path(output_dir, out_file))
    
    #  pull 20 most significant pathways only
    DEpwys <- toptable_output %>%
      arrange(adj.P.Val) %>%
      filter(adj.P.Val < 0.05) %>%
      head(20) %>%
      pull(Geneset) 
    all_pathways <- c(all_pathways, DEpwys)
  }
  
  # plot top 20 pathways
  all_pathways <- unique(all_pathways)
  DEpwys_es <- Biobase::exprs(gsva_eset[all_pathways, ])
  DEpwys_annot <- pData(gsva_eset[DEpwys,])
  DEpwys_annot['sample_id'] <- NULL
  DEpwys_annot$cluster_class <- as.character(DEpwys_annot$cluster_class)
  
  # create annotation for cluster class
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  l <- gg_color_hue(length(n_clusters))
  names(l) <- as.character(n_clusters)
  mycolors <- list()
  mycolors[['cluster_class']] <- l 
  
  # color palette for short histology
  palettes_dir <- "../../palettes/"
  palette_file <- file.path(palettes_dir, "short_histology_color_palette.tsv") %>% read_tsv()
  DEpwys_annot <- DEpwys_annot %>%
    rownames_to_column("Kids_First_Biospecimen_ID") %>%
    inner_join(input_clin, by = "Kids_First_Biospecimen_ID") %>%
    inner_join(palette_file, by = "short_histology") %>%
    dplyr::select(Kids_First_Biospecimen_ID, cluster_class, short_histology, hex_code) %>%
    column_to_rownames("Kids_First_Biospecimen_ID")
  
  # create annotation for short histology
  short_histology_palettes <- DEpwys_annot %>%
    dplyr::select(short_histology, hex_code) %>%
    unique()
  mycolors[['short_histology']] <- short_histology_palettes$hex_code
  names(mycolors[['short_histology']]) <- short_histology_palettes$short_histology
  
  # remove colors from annotation table
  DEpwys_annot$hex_code <- NULL
  
  pheatmap::pheatmap(DEpwys_es, scale = "row", 
                     treeheight_row = 10, 
                     treeheight_col = 10,
                     fontsize_row = 8,
                     fontsize = 8, 
                     show_colnames = F, 
                     cellwidth = 0.5, cellheight = 10,
                     annotation = DEpwys_annot, 
                     annotation_colors = mycolors, 
                     cluster_cols = cluster_tree, 
                     filename = file.path(output_dir, paste0(prefix, '_top20_pathways.tiff')), 
                     width = 12, height = 8)
}