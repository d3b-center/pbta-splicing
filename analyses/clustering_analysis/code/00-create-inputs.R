# Author: Komal S. Rathi
# Function: 
# 1) Remove histologies with <= 5 samples from input matrices
# 2) Create KEGG geneset

suppressPackageStartupMessages({
  library(tidyverse)
  library(msigdbr)
})

# read histologies file
clin_file <- "../../data/histologies.tsv" %>% read_tsv()
clin_file <- clin_file %>% 
  filter(experimental_strategy == "RNA-Seq")

# read splice dataset
splice_mat <- "input/pan_cancer_splicing_SE.gene.txt" %>% 
  read_tsv() %>%
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
saveRDS(splice_mat, file = 'input/non_expr_pan_cancer_splice_subset.rds')

# 2) read pbta dataset, filter to samples and genes in pbta splice matrix and save
pbta_subset <- "../../data/gene-counts-rsem-expected_count-collapsed.rds" %>% readRDS()
pbta_subset <- pbta_subset %>%
  dplyr::select(colnames(splice_mat)) %>%
  rownames_to_column("gene_symbol") %>%
  filter(gene_symbol %in% rownames(splice_mat)) %>%
  column_to_rownames("gene_symbol")
saveRDS(pbta_subset, file = 'input/raw_counts_pbta_subset.rds')

# 3) create KEGG input file
kegg_db <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
kegg_db <- kegg_db %>%
  dplyr::select(gene_symbol, gs_name) %>%
  unique() %>%
  unstack()
saveRDS(kegg_db, file = file.path('input/kegg_geneset_mrna.rds'))
