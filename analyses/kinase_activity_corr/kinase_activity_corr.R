# Author: Run Jin
#
# Plotting of heatmap of kinase activities with clustering rows and columns

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(pheatmap)
  library(readxl)
})


#### Define Directories --------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "kinase_activity_corr")

input_dir <- file.path(analysis_dir, "input")

plots_dir <- file.path(analysis_dir, "plots")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

#### Read in files necessary for analyses --------------------------------------
# histology file 
histology_df <- readr::read_tsv(file.path(data_dir, "pbta-histologies.tsv"), guess_max=10000)

# kinase activity score
kinase_activity <- readxl::read_excel(file.path(input_dir, "1-s2.0-S0092867420314513-mmc4.xlsx"), sheet =3, skip=2)

# PSI group 
psi_group <- readr::read_tsv(file.path(input_dir, "CLK1_PSI_grouping.txt"), skip=1)

# match BS_ID to sample ID
histology_matched <- histology_df %>% 
  dplyr::filter(Kids_First_Biospecimen_ID %in% psi_group$BS_id) %>%
  dplyr::select(Kids_First_Biospecimen_ID, sample_id, harmonized_diagnosis)

# note that 2 samples are not in the histologies file 
print(psi_group$BS_id[!psi_group$BS_id %in% histology_openpedcan$Kids_First_Biospecimen_ID])
# BS_G43E1QBM & BS_QV367NEX

# only 7 samples have kinase activity scores from the publication
common_samples <- histology_matched$sample_id[histology_matched$sample_id %in% colnames(kinase_activity)] %>%
  unique()

# generate an annotation file for the heatmap  
histology_matched_filtered <- histology_matched %>% 
  dplyr::filter(sample_id %in% common_samples) %>% 
  left_join( psi_group, by=c("Kids_First_Biospecimen_ID"= "BS_id")) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE) %>%
  dplyr::select(-"Kids_First_Biospecimen_ID") %>% 
  tibble::column_to_rownames("sample_id") %>%
  arrange(Group, harmonized_diagnosis)

# filter kinase activity
kinase_activity_filtered <- kinase_activity %>% 
  tibble::column_to_rownames("proID") 

# clear the NA in the numbers
fix_NA_function <- function(df){
  df_fixed <- apply(df, 2, function(x) as.numeric(gsub("NA", "", x)))
  rownames(df_fixed) <- rownames(df)
  return(df_fixed)
}

kinase_activity_fixed <-fix_NA_function(kinase_activity_filtered) %>% 
  as.data.frame() %>% 
  # order the columns for plotting
  dplyr::select(rownames(histology_matched_filtered)) 

# generate heatmap
pheatmap::pheatmap(as.matrix(kinase_activity_fixed),
                   annotation_col = histology_matched_filtered,
                   cluster_rows=FALSE, 
                   cluster_cols=FALSE,
                   width = 10, 
                   height = 8,
                   filename = file.path(plots_dir, "kinase_group_psi_hist.pdf"))






