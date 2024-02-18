## convert_to_rds.R
## Converts matrix tsv into RDS file for downstream clustering
## author: Ammar S Naqvi
################################################################################

library("tidyverse")

## directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "clustering_analysis")
results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")
input_dir <- file.path(analysis_dir, "input")

matrix_file <- file.path(input_dir,"pan_cancer_splicing_SE.gene.tsv")
rds_out <- file.path(input_dir,"pan_cancer_splicing_SE.gene.rds")

## get matrix and save to RDS file
splice_df <- read_tsv(matrix_file)
saveRDS(splice_df, rds_out)