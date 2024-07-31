## convert_to_rds.R
## Converts matrix tsv into RDS file for downstream clustering
## author: Jo Lynne Rokita
################################################################################

library("tidyverse")

## directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "clustering_analysis")
results_dir <- file.path(analysis_dir, "input", "functional-sites")
dir.create(results_dir, showWarnings = F, recursive = T)
input_dir <- file.path(root_dir, "analyses", "create-functional-site-matrix", "results")

# select independent specimens from PBTA primary tumors only
indep_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv") %>%
  read_tsv()
indep_file <- indep_file %>% 
  filter(cohort == "PBTA",
         experimental_strategy == "RNA-Seq") 

## get PSI matrix with functional sites, filter to independent samples and save to RDS file
matrix_file <- file.path(input_dir,"pan-cancer-SE-func.tsv.gz")
splice_mat <- read_tsv(matrix_file) %>%
  column_to_rownames("Splice_ID") %>%
  dplyr::select(any_of(indep_file$Kids_First_Biospecimen_ID))
  
# filter to OncoKB Cancer Gene List
oncokb_genelist <- file.path(analysis_dir, "input", "cancerGeneList.tsv")
oncokb_genelist <- read_tsv(oncokb_genelist)
splice_mat <- splice_mat %>%
  rownames_to_column("Splice_ID") %>%
  mutate(gene = gsub(":.*", "", Splice_ID)) %>%
  filter(gene %in% oncokb_genelist$`Hugo Symbol`) %>%
  dplyr::select(-c(gene)) %>%
  column_to_rownames("Splice_ID")

# save
rds_out <- file.path(results_dir,"pan_cancer_splicing_SE.gene.rds")
saveRDS(splice_mat, rds_out)
