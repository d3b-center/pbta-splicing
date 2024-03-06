# Author: Komal S. Rathi
# Function: Get optimal clustering output by evaluating all algorithms and distance combinations

# load libraries
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
})

# parse parameters     
option_list <- list(
  make_option(c("--input_mat"), type = "character",
              help = "Input matrix collapsed to unique features (.rds)"),
  make_option(c("--hist_file"), type = "character", default = NULL,
              help = "histology file (.tsv)"),
  make_option(c("--cluster_algorithm"), type = "character",
              help = "comma separated list of CCP clustering algorithms"),
  make_option(c("--cluster_distance"), type = "character",
              help = "comma separated list of CCP clustering distances"),
  make_option(c("--filter_expr"), type = "logical",
              help = "filter expression (TRUE or FALSE)"),
  make_option(c("--protein_coding_only"), type = "logical",
              help = "filter to protein coding genes only (TRUE or FALSE)"),
  make_option(c("--gencode_version"), type = "numeric", default = NULL,
              help = "gencode version used for protein coding genes filter (>=27)"),
  make_option(c("--feature_selection"), type = "character", 
              help = "feature selection type (dip.test or variance)"),
  make_option(c("--var_prop"), type = "numeric", default = NULL,
              help = "% of variable features (> 0)"),
  make_option(c("--transformation_type"), type = "character",
              help = "transformation type (none, tmm, vst, uq, log2, rank)"),
  make_option(c("--max_k"), type = "numeric",
              help = "maximum k-values to evaluate (>= 2)")
)
opt <- parse_args(OptionParser(option_list = option_list))
# input matrix
input_mat <- opt$input_mat
print(input_mat)

# histology file (optional)
hist_file <- opt$hist_file
print(hist_file)

# algorithms
cluster_algorithm <- opt$cluster_algorithm 
cluster_algorithm <- strsplit(cluster_algorithm, split = ", ") %>% unlist()
stopifnot(cluster_algorithm %in% c("hc", "pam", "km"))
print(cluster_algorithm)

# distances
cluster_distance <- opt$cluster_distance 
cluster_distance <- strsplit(cluster_distance, split = ", ") %>% unlist()
stopifnot(cluster_distance %in% c("pearson", "spearman", "euclidean", "manhattan", "binary", "maximum", "canberra", "minkowski"))
print(cluster_distance)

# filter by expression? 
filter_expr <- opt$filter_expr
print(filter_expr)

# filter to protein coding genes only?
protein_coding_only <- opt$protein_coding_only
print(protein_coding_only)

# gencode version for protein coding genes
gencode_version <- opt$gencode_version 
if(!is.null(gencode_version)){
  # if gencode version is not NULL, it should be at least at version 27 
  # gencode_version <- gencode_version %>% as.numeric()
  stopifnot(gencode_version >= 27) 
}
print(gencode_version)

# type of feature selection
feature_selection <- opt$feature_selection
stopifnot(feature_selection %in% c("variance", "dip.test"))
print(feature_selection)

# proportion of features sorted by variance
var_prop <- opt$var_prop
# if feature selection is variance, var_prop cannot be NULL and should be > 0
if(feature_selection == "variance"){
  stopifnot(!is.null(var_prop) & var_prop > 0)
  var_prop <- var_prop %>% as.numeric()
}
print(var_prop)

# transformation type
transformation_type <- opt$transformation_type
stopifnot(transformation_type %in% c("none", "tmm", "vst", "uq", "log2", "rank"))
print(transformation_type)

# max k-values to evaluate
max_k <- opt$max_k %>% as.numeric()
print(max_k)

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "clustering_analysis")
output_dir <- file.path(analysis_dir, "output", "optimal_clustering")
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(analysis_dir, 'util', 'lspline_clustering.R'))

# compute optimal clustering across all combinations of algorithms, distances and input k-values
  lspline_clustering(expr_mat = input_mat,
                     hist_file = hist_file,
                     algorithms = cluster_algorithm,
                     distances = cluster_distance,
                     filter_expr = filter_expr,
                     protein_coding_only = protein_coding_only,
                     gencode_version = gencode_version,
                     feature_selection = feature_selection,
                     var_prop = var_prop,
                     transformation_type = transformation_type,
                     max_k = max_k,
                     output_dir = output_dir)
