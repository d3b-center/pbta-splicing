## convert_to_rds.R
## Converts matrix tsv into RDS file for downstream clustering
## author: Jo Lynne Rokita
################################################################################

library("tidyverse")

## directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "clustering_analysis")
results_dir <- file.path(analysis_dir, "functional-sites", "input")
input_dir <- file.path(root_dir, "analyses", "create-functional-site-matrix", "results")

indep_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv")
matrix_file <- file.path(input_dir,"pan-cancer-SE-func.tsv.gz")
rds_out <- file.path(results_dir,"pan_cancer_splicing_SE_func.rds")
counts_rds_output <- file.path(input_dir,"raw_counts_pbta_subset.rds")

indep_specimens <- read_tsv(indep_file)

## get matrix and save to RDS file
splice_mat <- read_tsv(matrix_file) %>%
  column_to_rownames("Splice_ID") %>%
  # select independent specimens from primary tumors only
  dplyr::select(any_of(indep_specimens$Kids_First_Biospecimen_ID))
    
# save
saveRDS(splice_mat, rds_out)
