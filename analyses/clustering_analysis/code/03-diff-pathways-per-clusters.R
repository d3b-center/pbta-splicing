# Author: Komal S. Rathi
# Function: Identify differential pathway analysis for each cluster using GSVA

# load libraries
suppressPackageStartupMessages({
  library(optparse)
})

# source function to perform dip test
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "clustering_analysis")
source(file.path(analysis_dir, 'util', 'diff_pathways_per_cluster.R'))

# parse parameters     
option_list <- list(
  make_option(c("--input_mat"), type = "character",
              help = "Input matrix collapsed to unique features (.rds)"),
  make_option(c("--cluster_output"), type = "character",
              help = "cluster output obtained after running 01-get-clustering-output script"),
  make_option(c("--n_cluster"), type = "character",
              help = "number of clusters to investigate"),
  make_option(c("--gene_set"), type = "character",
              help = "gene set of interest (.rds)"),
  make_option(c("--prefix"), type = "character",
              help = "prefix of output files"),
  make_option(c("--output_dir"), type = "character",
              help = "output directory path")
)
opt <- parse_args(OptionParser(option_list = option_list))
input_mat <- opt$input_mat
cluster_output <- opt$cluster_output
n_cluster <- opt$n_cluster
n_cluster <- as.numeric(n_cluster)
gene_set <- opt$gene_set
prefix <- opt$prefix
output_dir <- opt$output_dir

# differential pathway expression per cluster using GSVA
diff_pathways_per_cluster(input_mat =  input_mat,
                          cluster_output = cluster_output, 
                          n_cluster = n_cluster, 
                          gene_set = gene_set, 
                          prefix = prefix, 
                          output_dir = output_dir)
