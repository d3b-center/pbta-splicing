# Author: Komal S. Rathi, Jo Lynne Rokita
# Function: 
# 1) Create KEGG geneset

suppressPackageStartupMessages({
  library(tidyverse)
  library(msigdbr)
  library(vroom)
})


#Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "clustering_analysis")
data_dir <- file.path(root_dir, "data")

results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")
hist_dir <- file.path(root_dir, "analyses", "cohort_summary", "results")
input_dir <- file.path(analysis_dir, "input")


# read histologies file
clin_file <- read_tsv(file.path(hist_dir, "histologies-plot-group.tsv")) %>% 
  filter(experimental_strategy == "RNA-Seq")

# read splice dataset
splice_mat <- readRDS(file.path(input_dir, "pan_cancer_splicing_SE.gene.rds"))

# 2) read pbta dataset, filter to samples and genes in pbta splice matrix and save
pbta_subset <- readRDS(file.path(data_dir,"gene-counts-rsem-expected_count-collapsed.rds")) %>%
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
