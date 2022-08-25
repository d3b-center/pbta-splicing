# Author: Komal S. Rathi
# Function: Identify differential pathway analysis for each cluster using GSVA

# load libraries
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(limma)
  library(EnhancedVolcano)
  library(GSVA)
  library(Biobase)
  library(pheatmap)
  library(DESeq2)
})

# function to identify differential genes per cluster and do pre-ranked gsea
diff_pathways_per_cluster <- function(input_mat, cluster_output, n_cluster, gene_set, prefix, output_dir){
  
  # output directory
  dir.create(output_dir, showWarnings = F, recursive = T)
  
  # read data
  input_mat <- readRDS(input_mat)
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
  DEpwys_es <- exprs(gsva_eset[all_pathways, ])
  DEpwys_annot <- pData(gsva_eset[DEpwys,])
  DEpwys_annot['sample_id'] <- NULL
  DEpwys_annot$cluster_class <- as.character(DEpwys_annot$cluster_class)
  
  # annotation colors
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  l <- gg_color_hue(length(n_clusters))
  names(l) <- as.character(n_clusters)
  mycolors <- list()
  mycolors[['cluster_class']] <- l 
  
  pheatmap::pheatmap(DEpwys_es, scale = "row", 
                     show_colnames = F, 
                     cellwidth = 0.5, cellheight = 10,
                     annotation = DEpwys_annot, 
                     annotation_colors = mycolors, 
                     cluster_cols = cluster_tree, 
                     filename = file.path(output_dir, paste0(prefix, '_top20_pathways.pdf')), 
                     width = 12, height = 15)
}

# differential pathway expression per cluster on pbta splicing data
diff_pathways_per_cluster(input_mat =  'output/ccp_output/non_expr_pan_cancer_splice_subset_km_euclidean_0_matrix.rds',
                          cluster_output = 'output/ccp_output/non_expr_pan_cancer_splice_subset_km_euclidean_0_ccp.rds', 
                          n_cluster = 3, 
                          gene_set = 'input/kegg_geneset_mrna.rds', 
                          prefix = 'non_expr_pan_cancer_splice_subset_km_euclidean_0', 
                          output_dir = 'output/diff_pathways/')

# differential pathway expression per cluster on pbta expression data
diff_pathways_per_cluster(input_mat =  'output/ccp_output/raw_counts_pbta_subset_km_euclidean_0_matrix.rds',
                          cluster_output = 'output/ccp_output/raw_counts_pbta_subset_km_euclidean_0_ccp.rds', 
                          n_cluster = 3, 
                          gene_set = 'input/kegg_geneset_mrna.rds', 
                          prefix = 'raw_counts_pbta_subset_km_euclidean_0', 
                          output_dir = 'output/diff_pathways/')

