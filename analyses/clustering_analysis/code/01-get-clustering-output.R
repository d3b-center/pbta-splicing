# Author: Komal S. Rathi
# Function: Get ccp clustering output for a specific combination of distance + algorithm + % variable genes 

# load libraries
suppressPackageStartupMessages({
  library(optparse)
})

# source function to perform dip test
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "clustering_analysis")
source(file.path(analysis_dir, 'util', 'perform_diptest.R'))
source(file.path(analysis_dir, 'util', 'get_ccp_output.R'))
         
# parse parameters     
option_list <- list(
  make_option(c("--input_mat"), type = "character",
              help = "Input matrix collapsed to unique features (.rds)"),
  make_option(c("--data_type"), type = "character",
              help = "data type of the input matrix (raw_counts or non_expr)"),
  make_option(c("--var_genes"), type = "character",
              help = "% variable genes"),
  make_option(c("--cluster_algorithm"), type = "character",
              help = "CCP clustering algorithm"),
  make_option(c("--cluster_distance"), type = "character",
              help = "CCP clustering distance"),
  make_option(c("--prefix"), type = "character",
              help = "prefix of output files")
)
opt <- parse_args(OptionParser(option_list = option_list))
input_mat <- opt$input_mat
input_mat <- readRDS(input_mat)
data_type <- opt$data_type
var_genes <- opt$var_genes
var_genes <- as.numeric(var_genes)
cluster_algorithm <- opt$cluster_algorithm
cluster_distance <- opt$cluster_distance
prefix <- opt$prefix

# get ccp clustering output for a specific combination of distance + algorithm + % variable genes 
get_ccp_output(input_mat,
               data_type = data_type, 
               var_genes = var_genes, 
               cluster_algorithm = cluster_algorithm, 
               cluster_distance = cluster_distance, 
               prefix = prefix)
