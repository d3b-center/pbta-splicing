# Author: Komal S. Rathi
# Function: Get ccp clustering output for a specific combination of distance + algorithm + % variable genes 
# for the splicing dataset, we chose km, euclidean at 0% (based on Ammar's manual inspection)

# load libraries
suppressPackageStartupMessages({
  library(optparse)
  library(ConsensusClusterPlus)
  library(tidyverse)
  library(DGCA)
  library(DESeq2)
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

# pbta splicing data
# read input file
# input_mat <- readRDS('input/non_expr_pan_cancer_splice_subset.rds')
get_ccp_output(input_mat,
               data_type = data_type, 
               var_genes = var_genes, 
               cluster_algorithm = cluster_algorithm, 
               cluster_distance = cluster_distance, 
               prefix = prefix)

# run function with same clustering parameters using pbta expression data
# pbta expression data subsetted to samples in the splicing data
# read input file
# input_mat <- readRDS('input/raw_counts_pbta_subset.rds')
# get_ccp_output(input_mat,
#                data_type = "raw_counts", 
#                var_genes = 0, 
#                cluster_algorithm = "km", 
#                cluster_distance = "euclidean", 
#                prefix = "raw_counts_pbta_subset")
