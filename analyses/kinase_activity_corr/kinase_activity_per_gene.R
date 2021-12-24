# Author: Run Jin
# Generate kinase activity map with PSI scores from various genes 
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(pheatmap)
  library(readxl)
  library(optparse)
  library(Hmisc)
  library(RColorBrewer)
  library(magrittr)
  library(corrplot)
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

plots_sig_dir <- file.path(analysis_dir, "plots", "individual_genes", "significant_heatmap")
if(!dir.exists(plots_sig_dir)){
  dir.create(plots_sig_dir, recursive=TRUE)
}

corr_plots_dir <- file.path(analysis_dir, "plots", "individual_genes", "corrplot")
if(!dir.exists(corr_plots_dir)){
  dir.create(corr_plots_dir, recursive=TRUE)
}

results_dir <- file.path(analysis_dir, "results")
if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive=TRUE)
}

#### Read in files necessary for analyses
# histology file 
histology_df <- readr::read_tsv(opt$histology, guess_max=10000)

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
    dplyr::select(Kids_First_Biospecimen_ID, sample_id, short_histology, composition, RNA_library)
  
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
  
  
  ################### Additionally, calculate correlation 
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
  res <- rcorr(as.matrix(kinase_activity_calc))
  cor_sig_genes <- res$P %>%
    as.data.frame() %>%
    dplyr::select(PSI) %>%
    dplyr::filter(PSI < 0.05) %>%
    rownames() 
           
  ##### plot out heatmaps with only sig. correlated genes for SRRM1 only
  if(length(cor_sig_genes) >=1 && grepl("SRRM1", file_id)){
    
    # order by PSI and then short histology
    plot_matrix <-histology_matched_filtered %>% 
      arrange(PSI, short_histology)
    
    plot_matrix$short_histology <- factor(plot_matrix$short_histology)
    
    # generate color list for heatmaps
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    # select n distinct colors short histology
    n_short_hist <- plot_matrix %>% pull(short_histology) %>% unique() %>% length()
    # generate a list of colors for each annotation 
    set.seed(1015)
    short_hist_color <- sample(col_vector, n_short_hist) %>% 
      set_names(unique(plot_matrix$short_histology))
    
    anno_colors <- list(short_hist_color)
    names(anno_colors) <- c("short_histology")
    
    # order the columns to plotting
    kinase_activity_plot <-kinase_activity %>% 
      tibble::rownames_to_column("proID") %>%
      dplyr::filter(proID %in% cor_sig_genes) %>%
      tibble::column_to_rownames("proID") %>%
      # order the columns for plotting
      dplyr::select(rownames(plot_matrix)) 
    
    pheatmap::pheatmap(as.matrix(kinase_activity_plot),
                       annotation_col = plot_matrix,
                       cluster_rows=(length(cor_sig_genes) >=2),
                       cluster_cols=FALSE,
                       annotation_colors = anno_colors,
                       width = 10,
                       height = 8,
                       show_colnames = F,
                       filename = file.path(plots_sig_dir, paste0(file_id, "_vs_kinase_by_psi_short.pdf")))
  }
  
  ####### plot out correlation graph
  pdf(file = file.path(corr_plots_dir, paste0(file_id, "_corrplot.pdf")))
  p <- corrplot(as.matrix(res[["r"]]), 
                type="lower", 
                order="hclust",   
                bg="transparent", 
                sig.level = 0.05,
                insig = "blank")
  print(p)
  dev.off()
  
  #### output the file with P and R value 
  # get r score
  r_scores <- res$r %>% 
    as.data.frame() %>%
    dplyr::select(PSI) %>% 
    tibble::rownames_to_column("Kinase_name") %>%
    dplyr::rename(r_scores = PSI)
  
  # get p score
  p_scores <- res$P %>% 
    as.data.frame() %>%
    dplyr::select(PSI) %>% 
    tibble::rownames_to_column("Kinase_name") %>%
    dplyr::rename(p_scores = PSI)
  
  # combined r and p and add the name of the gene and write out results
  left_join(r_scores, p_scores) %>% 
    dplyr::mutate(RBP_gene_name = file_id) %>% 
    readr::write_tsv(file.path(results_dir, paste0(file_id, "_correlation_results.tsv")))
}



