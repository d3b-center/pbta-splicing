# Author: Komal S. Rathi
# Function: 
# 1) Remove histologies with <= 5 samples from input matrices
# 2) Create KEGG geneset

suppressPackageStartupMessages({
  library(tidyverse)
  library(msigdbr)
})


#Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "clustering_analysis")
data_dir <- file.path(root_dir, "data")

results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")
input_dir <- file.path(analysis_dir, "input")


# read histologies file
clin_file <- file.path(input_dir, "histologies-plot-group.tsv") %>% read_tsv() %>% 
  dplyr::select(-short_histology) %>% 
  dplyr::rename(short_histology=plot_group) %>% 
  filter(experimental_strategy == "RNA-Seq")

# read splice dataset
splice_mat <- readRDS(file.path(input_dir, "pan_cancer_splicing_SE.gene.rds")) %>%
  column_to_rownames("Splice_ID") 

# remove histologies with <= 5 samples
remove_sh <- clin_file %>%
  filter(Kids_First_Biospecimen_ID %in% colnames(splice_mat)) %>%
  group_by(short_histology) %>%
  tally() %>%
  filter(n <= 5) %>% 
  pull(short_histology)

clin_file <- clin_file %>%
  filter(!short_histology %in% remove_sh)

# 1) remove samples from splice matrix and save
splice_mat <- splice_mat %>%
  dplyr::select(any_of(clin_file$Kids_First_Biospecimen_ID))

## set rownames
saveRDS(splice_mat, file = file.path(input_dir, 'non_expr_pan_cancer_splice.rds'))


# 2) read pbta dataset, filter to samples and genes in pbta splice matrix and save
pbta_subset <- file.path(data_dir,"gene-counts-rsem-expected_count-collapsed.rds") %>% readRDS()
pbta_subset <- pbta_subset %>%
  dplyr::select(colnames(splice_mat)) %>%
  rownames_to_column("gene_symbol") %>%
  filter(gene_symbol %in% rownames(splice_mat)) %>%
  column_to_rownames("gene_symbol")

counts_rds_output <- file.path(input_dir,"raw_counts_pbta_subset.rds")
saveRDS(pbta_subset, file = counts_rds_output)

# 3) create KEGG input file
kegg_db <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
kegg_db <- kegg_db %>%
  dplyr::select(gene_symbol, gs_name) %>%
  unique() %>%
  unstack()

keggdb_rds_output <- file.path(input_dir,"kegg_geneset_mrna.rds")
saveRDS(kegg_db, file = keggdb_rds_output)
