# Author: Komal S. Rathi
# Function: Identify differential gene expression and pre-ranked gsea per cluster

# load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(EnhancedVolcano)
  library(Biobase)
  library(ggpubr)
  library(ggplot2)
  library(reshape2)
  library(fgsea)
  library(DESeq2)
})

# function to identify differential genes per cluster and do pre-ranked gsea
diff_genes_per_cluster <- function(input_mat, cluster_output, n_cluster, gene_set, prefix, output_dir){
  
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
  
  # differential expression across clusters
  # per cluster analysis
  n_clusters <- unique(cluster_output$cluster_class)
  for(i in 1:length(n_clusters)){
    
    # cluster of interest vs others
    cluster_of_interest <- n_clusters[i]
    type <- ifelse(cluster_output[["cluster_class"]] == cluster_of_interest, "ref", "others")
    
    # create design and fit
    mod <- model.matrix(~ factor(type))
    colnames(mod) <- c("others", "ref_vs_others")
    rownames(mod) <- cluster_output[["sample_id"]]
    fit <- limma::lmFit(object = expr_mat, design = mod)
    fit <- limma::eBayes(fit)
    toptable_output <- limma::topTable(fit, coef = 2, n = Inf)
    toptable_output <- toptable_output %>%
      rownames_to_column("Gene")
    
    # save filtered output
    out_file <- paste0(prefix, "_cluster_", i, "_limma_output.tsv")
    readr::write_tsv(x = toptable_output %>% filter(adj.P.Val < 0.05), 
                     file = file.path(output_dir, out_file))
    
    # pathway analysis using significantly differential expressed genes
    DE_genes <- toptable_output %>%
      filter(adj.P.Val < 0.05)
    ranks <- DE_genes$logFC
    names(ranks) <- DE_genes$Gene
    ranks <- ranks[order(ranks, decreasing = TRUE)]
    set.seed(42)
    fgseaRes <- fgsea(gene_set, ranks, minSize = 1, maxSize = Inf)
    
    # significant pathways
    fgseaRes <- fgseaRes %>%
      filter(padj < 0.05) %>%
      mutate(direction = ifelse(ES > 0, "up", "down")) %>%
      arrange(padj) %>%
      mutate(genes = sapply(leadingEdge, paste, collapse=",")) %>%
      dplyr::select(pathway, pval, padj, direction, genes)
    
    # save filtered output
    out_file <- paste0(prefix, "_cluster_", i, "_limma_output_pathway_ranked.tsv")
    readr::write_tsv(x = fgseaRes, file = file.path(output_dir, out_file))
  }
}

# differential gene expression per cluster + pre-ranked gsea on pbta splicing data
diff_genes_per_cluster(input_mat =  'output/ccp_output/non_expr_pan_cancer_splice_subset_km_euclidean_0_matrix.rds',
                       cluster_output = 'output/ccp_output/non_expr_pan_cancer_splice_subset_km_euclidean_0_ccp.rds', 
                       n_cluster = 3, 
                       gene_set = 'input/kegg_geneset_mrna.rds', 
                       prefix = 'non_expr_pan_cancer_splice_subset_km_euclidean_0', 
                       output_dir = 'output/diff_genes/')

# differential gene expression per cluster + pre-ranked gsea on pbta expression data
diff_genes_per_cluster(input_mat =  'output/ccp_output/raw_counts_pbta_subset_km_euclidean_0_matrix.rds',
                       cluster_output = 'output/ccp_output/raw_counts_pbta_subset_km_euclidean_0_ccp.rds', 
                       n_cluster = 3, 
                       gene_set = 'input/kegg_geneset_mrna.rds', 
                       prefix = 'raw_counts_pbta_subset_km_euclidean_0', 
                       output_dir = 'output/diff_genes/')

