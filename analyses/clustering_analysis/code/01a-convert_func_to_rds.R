## convert_to_rds.R
## Converts matrix tsv into RDS file for downstream clustering
## author: Jo Lynne Rokita
################################################################################

library("tidyverse")

## directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "clustering_analysis")
results_dir <- file.path(analysis_dir, "output", "functional-sites")
input_dir <- file.path(root_dir, "analyses", "create-functional-site-matrix", "results")

matrix_file <- file.path(input_dir,"pan-cancer-SE-func.tsv.gz")
rds_out <- file.path(results_dir,"pan_cancer_splicing_SE.gene.rds")
counts_rds_output <- file.path(input_dir,"raw_counts_pbta_subset.rds")

## get matrix and save to RDS file
splice_mat <- read_tsv(matrix_file) %>%
  column_to_rownames("Splice_ID")

# save
saveRDS(splice_mat, rds_out)

# read pbta dataset, filter to samples and genes in pbta splice matrix and save
pbta_subset <- readRDS(file.path(data_dir,"gene-counts-rsem-expected_count-collapsed.rds")) %>%
  dplyr::select(any_of(colnames(splice_mat))) %>%
  rownames_to_column("gene_symbol") %>%
  filter(gene_symbol %in% rownames(splice_mat)) %>%
  column_to_rownames("gene_symbol")

# save
saveRDS(pbta_subset, file = counts_rds_output)
