# test Ammar's dataset
# top 10 pathways per cluster - 1 heatmap

setwd('~/Projects/OpenPBTA-miRNA/')
suppressPackageStartupMessages({
  library(limma)
  library(GSVA)
  library(EnhancedVolcano)
  library(pheatmap)
  library(tidyverse)
})


source('R/perform_diptest.R')
expr_mat = readRDS('data/non_expr_pan_cancer_splice_subset.rds')
res.ccp = readRDS('www/non_expr_pan_cancer_splice_subset_km_euclidean_0.rds')
ccp_tree <- res.ccp[[n_cluster]]$consensusTree 
n_cluster = 3
cluster_output = res.ccp[[n_cluster]]$consensusClass %>%
  as.data.frame() %>%
  rownames_to_column('sample_id') %>%
  dplyr::rename('cluster_class' = '.')

# remove rows where all values are 0
mat <- expr_mat
mat <- mat[rowSums(mat[])>0,]

# remove columns and rows with sd == 0 
mat <- mat[apply(mat, MARGIN = 1,function(x) { sd(x) != 0} ),] 
mat <- mat[,apply(mat, MARGIN = 2,function(x) { sd(x) != 0} )] 

# var genes = 0 (perform diptest)
mat <- perform_diptest(count_matrix = mat)
expr_mat <- mat

# convert into an eset with cluster info
cluster_output <- cluster_output %>%
  mutate(tmp = sample_id) %>%
  column_to_rownames('tmp')
cluster_output <- cluster_output[colnames(expr_mat),]
eset <- Biobase::ExpressionSet(assayData = as.matrix(expr_mat), 
                               phenoData = new("AnnotatedDataFrame", data = cluster_output))

# GSVA analysis
geneset = 'data/kegg_geneset_mrna.rds'
geneset <- readRDS(geneset)
gsva_eset <- GSVA::gsva(expr = eset,
                        gset.idx.list = geneset, 
                        method = "ssgsea", 
                        kcdf = "Gaussian") 

# convert to long format
gsva_scores <- gsva_eset@assayData$exprs %>%
  as.data.frame() %>%
  rownames_to_column("geneset") %>%
  gather("sample_id", "score", -c("geneset")) %>%
  dplyr::select(sample_id, geneset, score)

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
  
  #  pull 20 most significant pathways only
  # first try adjusted p-values
  DEpwys <- toptable_output %>%
    arrange(adj.P.Val) %>%
    filter(adj.P.Val < 0.05) %>%
    head(20) %>%
    pull(Geneset) 
  all_pathways <- c(all_pathways, DEpwys)
}

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
                   cluster_cols = ccp_tree, 
                   filename = 'www/non_expr_pan_cancer_splice_subset_km_euclidean_0_top20_pathways.pdf', 
                   width = 12, height = 15)
