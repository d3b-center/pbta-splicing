# Author: Run Jin
# Generate kinase activity map with PSI scores from various genes 
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(pheatmap)
  library(readxl)
  library(optparse)
  library(Hmisc)
})

#### Parse command line options ------------------------------------------------
option_list <- list(
  make_option(c("-h", "--histology"),type="character",
              help="histology file for all OpenPedCan samples (.tsv) "),
  make_option(c("-k","--kinase"),type="character",
              help="Kinase activity scores for all samples (.xlsx) "),
  make_option(c("-g","--group"),type="character",
              help="File containing the PSI groups based on CLK1 psi scores (.tsv) "),
  make_option(c("-d","--pki_input_dir"),type="character",
              help="directory to files of PSI scores")
)
opt <- parse_args(OptionParser(option_list=option_list,add_help_option = FALSE))
pki_input_dir <- opt$pki_input_dir

#### Define Directories 
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "kinase_activity_corr")

plots_dir <- file.path(analysis_dir, "plots", "individual_genes", "all")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

plots_sig_dir <- file.path(analysis_dir, "plots", "individual_genes", "significant")
if(!dir.exists(plots_sig_dir)){
  dir.create(plots_sig_dir, recursive=TRUE)
}

#### Read in files necessary for analyses
# histology file 
histology_df <- readr::read_tsv(opt$histology, guess_max=10000)

# PSI group by CLK1 file
psi_group_CLK1 <- readr::read_tsv(opt$group) %>%
  dplyr::rename(Kids_First_Biospecimen_ID = samples,
                PSI_group_by_CLK1 = Group) %>%
  dplyr::select(-PSI)
psi_group_CLK1$PSI_group_by_CLK1 <- factor(psi_group_CLK1$PSI_group_by_CLK1, 
                                           levels = c("Low","Medium","High"))

# kinase activity score
kinase_activity <- readxl::read_excel(opt$kinase, sheet =3, skip=2) %>%
  tibble::column_to_rownames("proID") %>% 
  dplyr::select(-geneID)

# read in the files with PSI scores
file_list <- list.files(pki_input_dir)

for (i in 1:length(file_list)){
 
  # get the file of interest
  file_of_interest <- file_list[i]
  # read in the file of interest
  psi_of_interest <- readr::read_tsv(file.path(pki_input_dir, file_of_interest), col_names =FALSE)
  colnames(psi_of_interest) <- c("Kids_First_Biospecimen_ID", "PSI")
  
  # get the name to output
  file_id <- gsub("-.*", "", file_of_interest)
  
  # match BS ID to sample ID
  histology_matched <- histology_df %>% 
    dplyr::filter(Kids_First_Biospecimen_ID %in% psi_of_interest$Kids_First_Biospecimen_ID) %>%
    # add the CLK1 PSI group here
    dplyr::left_join(psi_group_CLK1) %>% 
    dplyr::select(Kids_First_Biospecimen_ID, sample_id, harmonized_diagnosis, short_histology, composition, RNA_library, PSI_group_by_CLK1)
  
  # filter to only contain common samples
  common_samples <- histology_matched$sample_id[histology_matched$sample_id %in% colnames(kinase_activity)] %>%
    unique()
  histology_matched_filtered <- histology_matched %>% 
    dplyr::filter(sample_id %in% common_samples) %>%
    dplyr::left_join(psi_of_interest)
  
  # find duplicated samples
  duplicated_sample_id <- histology_matched_filtered$sample_id[duplicated(histology_matched_filtered$sample_id)]
  histology_matched_filtered_dup <-histology_matched_filtered %>% filter(sample_id %in% duplicated_sample_id) 
  
  # filter for distinct stranded, solid tissue
  histology_matched_dup_match <- histology_matched_filtered_dup %>%
    dplyr::filter(RNA_library == "stranded",
                  composition == "Solid Tissue") %>% 
    dplyr::distinct(sample_id, .keep_all = TRUE) 
  
  histology_matched_dup_non_match <- histology_matched_filtered_dup %>%
    dplyr::filter(!sample_id %in% histology_matched_dup_match$sample_id) %>% 
    dplyr::distinct(sample_id, .keep_all = TRUE) 
  
  histology_matched_dup_fixed <- bind_rows(histology_matched_dup_match,
                                           histology_matched_dup_non_match)
  
  # combine duplicated with non-duplicated to get a complete df 
  histology_matched_filtered <- histology_matched_filtered %>%
    dplyr::filter(!sample_id %in% duplicated_sample_id) %>% 
    bind_rows(histology_matched_dup_fixed) %>%
    tibble::column_to_rownames("sample_id") %>%
    dplyr::select(-c("RNA_library", "composition", "Kids_First_Biospecimen_ID"))
  
  # order by PSI and then short histology
  plot_matrix_1 <-histology_matched_filtered %>% 
    arrange(PSI, short_histology)
  
  # order the columns to plotting
  kinase_activity_plot1 <-kinase_activity %>% 
    # order the columns for plotting
    dplyr::select(rownames(plot_matrix_1)) 
  
  pheatmap::pheatmap(as.matrix(kinase_activity_plot1),
                     annotation_col = plot_matrix_1,
                     cluster_rows=TRUE,
                     cluster_cols=FALSE,
                     width = 10,
                     height = 8,
                     show_colnames = F,
                     filename = file.path(plots_dir, paste0(file_id, "_vs_kinase_by_psi_short.pdf")))
  
  # order by short histology and then PSI
  plot_matrix_2 <-histology_matched_filtered %>% 
    arrange(short_histology, PSI)
  
  # order the columns to plotting
  kinase_activity_plot2 <-kinase_activity %>% 
    # order the columns for plotting
    dplyr::select(rownames(plot_matrix_2)) 
  
  pheatmap::pheatmap(as.matrix(kinase_activity_plot2),
                     annotation_col = plot_matrix_2,
                     cluster_rows=TRUE,
                     cluster_cols=FALSE,
                     width = 10,
                     height = 8,
                     show_colnames = F,
                     filename = file.path(plots_dir, paste0(file_id, "_vs_kinase_by_short_psi.pdf")))
  
  ################### Additionally, plot out heatmaps with only sig. correlated genes
  histology_matched_filtered_calc <- histology_matched_filtered %>%
    tibble::rownames_to_column("sample_id") %>%
    dplyr::select(sample_id, PSI)
  
  kinase_activity_calc <- kinase_activity %>%
    dplyr::select(histology_matched_filtered_calc$sample_id) %>%
    t() %>% as.data.frame() %>%
    tibble::rownames_to_column("sample_id") %>%
    left_join(histology_matched_filtered_calc) %>%
    tibble::column_to_rownames("sample_id")
  
  # get the correlation results 
  cor_sig_genes <- rcorr(as.matrix(kinase_activity_calc))$P %>%
    as.data.frame() %>%
    dplyr::select(PSI) %>%
    dplyr::filter(PSI < 0.05) %>%
    rownames() 
  
  if(length(cor_sig_genes) >=1){
      
    # order by PSI and then short histology
    plot_matrix_3 <-histology_matched_filtered %>% 
      arrange(PSI, short_histology)
    
    # order the columns to plotting
    kinase_activity_plot3 <-kinase_activity %>% 
      tibble::rownames_to_column("proID") %>%
      dplyr::filter(proID %in% cor_sig_genes) %>%
      tibble::column_to_rownames("proID") %>%
      # order the columns for plotting
      dplyr::select(rownames(plot_matrix_1)) 
    
    pheatmap::pheatmap(as.matrix(kinase_activity_plot3),
                       annotation_col = plot_matrix_3,
                       cluster_rows=(length(cor_sig_genes) >=2),
                       cluster_cols=FALSE,
                       width = 10,
                       height = 8,
                       show_colnames = F,
                       filename = file.path(plots_sig_dir, paste0(file_id, "_vs_kinase_by_psi_short.pdf")))
    
    # order by short histology and then PSI
    plot_matrix_4 <-histology_matched_filtered %>% 
      arrange(short_histology, PSI)
    
    # order the columns to plotting
    kinase_activity_plot4 <-kinase_activity %>% 
      tibble::rownames_to_column("proID") %>%
      dplyr::filter(proID %in% cor_sig_genes) %>%
      tibble::column_to_rownames("proID") %>%
      # order the columns for plotting
      dplyr::select(rownames(plot_matrix_4)) 
    
    pheatmap::pheatmap(as.matrix(kinase_activity_plot4),
                       annotation_col = plot_matrix_4,
                       cluster_rows=(length(cor_sig_genes) >=2),
                       cluster_cols=FALSE,
                       width = 10,
                       height = 8,
                       show_colnames = F,
                       filename = file.path(plots_sig_dir, paste0(file_id, "_vs_kinase_by_short_psi.pdf")))
    
  }
}


