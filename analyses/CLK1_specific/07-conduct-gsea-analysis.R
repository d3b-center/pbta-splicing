################################################################################
# 07-conduct_gsea_analysis.R
# Author: Stephanie J. Spielman for CCDL ALSF, 2020 and modified by Ammar S 
# Naqvi
# Description: This script conducts gene set enrichment analysis, using the GSVA 
# method [1] for scoring hallmark human pathway enrichment from RNA-Seq results.
#
# The GSVA scores (i.e., enrichment scores) are calculated to produce a 
# **Gaussian distribution of scores** "under the null hypothesis of no change 
# in the pathway activity throughout the sample population."
#
# The authors claim a benefit to this approach:
#   + "Penalizes deviations that are large in *both* tails"
#   + "Provides a 'normalization' of the enrichment score by 
#      subtracting potential noise
#   + "Emphasizes genes in pathways that are concordantly activated in one 
#      direction only"
#   + "For pathways containing genes strongly acting in both directions, the 
#      deviations with cancel each other out and show little or no enrichment."
#
#
# usage: Rscript 07-conduct_gsea_analysis.R
#
# Reference:
# [1] Sonja Hänzelmann, Robert Castelo, and Justin Guinney. 2013. “GSVA: Gene Set 
# Variation Analysis for Microarray and RNA-Seq Data.” BMC Bioinformatics 14 
# (1): 7. https://doi.org/10.1186/1471-2105-14-7.
################################################################################

## Load libraries
suppressPackageStartupMessages({
  # Library to perform GSVA
  library("GSVA")
  # Library for data manipulation
  library("tidyverse")
  # Library to get genesets from MSigDB
  library("msigdbr")
})


## Set directories
# Input directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1_specific")
input_dir   <- file.path(analysis_dir, "input")
# Output directories
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

## Load input files
# Specify gene counts data path
expr_path <- file.path(data_dir, "gene-counts-rsem-expected_count-collapsed.rds")
# Load expression data
expression_data <- readr::read_rds(expr_path) %>%
  as.data.frame()
# Load Hallmark gene sets from MSigDB
human_hallmark <- msigdbr(species = "Homo sapiens", category = "H") # The loaded data is a tibble.
# Load histology data
histology_df <- readr::read_tsv( file.path(data_dir, "histologies.tsv")  , guess_max = 100000)

#### Prepare hallmark genes: Create a list of hallmarks, each of which is a list of genes -----------------------------------------------
human_hallmark_twocols <- human_hallmark %>% 
  dplyr::select(gs_name, human_gene_symbol)
human_hallmark_list <- base::split(human_hallmark_twocols$human_gene_symbol, list(human_hallmark_twocols$gs_name))

#### Perform gene set enrichment analysis --------------------------------------------------------------------
# Prepare expression data: log2 transform and re-cast as matrix

# Filter histology dataframe
histology_rna_df <- histology_df %>% 
  # Only include RNA-Seq samples 
  filter(experimental_strategy == "RNA-Seq") %>% 
  # Remove NAs
  filter(!is.na(RNA_library)) %>%
  # Only include samples from PBTA cohort
  filter(cohort == "PBTA") %>%
  # Include samples that are "stranded"
  filter(RNA_library == "stranded") %>% 
  # Include those samples originated from the Midline
  filter(CNS_region == 'Midline') %>% 
  # Include only HGAT sample histology
  filter(short_histology == 'HGAT') %>%
  # Filter for module-specific samples that was previously used in CLK1 splicing/comparisons
  filter((Kids_First_Biospecimen_ID   == 'BS_Q13FQ8FV') | 
         (Kids_First_Biospecimen_ID   == 'BS_ZV1P6W9C') |
         (Kids_First_Biospecimen_ID   == 'BS_WH8G4VFB') | 
         (Kids_First_Biospecimen_ID   == 'BS_NNPEC7W1') |
         (Kids_First_Biospecimen_ID   == 'BS_PZVHMSYN') | 
         (Kids_First_Biospecimen_ID   == 'BS_DRY58DTF') | 
         (Kids_First_Biospecimen_ID   == 'BS_GXTFW99H') | 
         (Kids_First_Biospecimen_ID   == 'BS_E60JZ9Z3') |
         (Kids_First_Biospecimen_ID   == 'BS_9CA93S6D') )


# Subset expression data to include only samples from filtered histology filter dataset
expression_data <- expression_data %>%
  select(histology_rna_df$Kids_First_Biospecimen_ID)
 
# For each type of the RNA library, we subset the expression matrix accordingly and run gsea scores for each RNA library 
rna_library_list <- histology_rna_df %>% 
  pull(RNA_library) %>% 
  unique()

# Further subset to each cohort to deal with size issues
cohort_list <- histology_rna_df %>% 
  pull(cohort) %>% 
  unique()

# Create empty dataframe
gsea_scores_df_tidy <- data.frame()

# Iterate through each cohort and RNA library type 
for(i in 1:length(rna_library_list)) {
  rna_library = rna_library_list[i]
  # Get biospecimen id for one particular RNA library type
  rna_library_type_bs_id <- histology_rna_df %>% 
    filter(RNA_library == rna_library) %>% 
    pull(Kids_First_Biospecimen_ID) %>%
    unique()
  
  # Filter the expression data to this RNA library type
  # Subset to the remaining samples 
  expression_data_each <- expression_data %>% 
    select(rna_library_type_bs_id)
  
  ### Rownames are genes and column names are samples
  expression_data_each_log2_matrix <- as.matrix(log2(expression_data_each + 1))
  
  # Calculate the Gaussian-distributed scores
  gsea_scores_each <- gsva(expression_data_each_log2_matrix,
                           human_hallmark_list,
                           method = "gsva",
                           min.sz=1, max.sz=1500,## Arguments from K. Rathi
                           parallel.sz = 8, # For the bigger dataset, this ensures this won't crash due to memory problems
                           mx.diff = TRUE) # Setting this argument to TRUE computes Gaussian-distributed scores (bimodal score distribution if FALSE)
  
  ### Clean scoring into tidy format
  gsea_scores_each_df <- as.data.frame(gsea_scores_each) %>%
    rownames_to_column(var = "hallmark_name")
  
  #first/last_bs needed for use in gather (we are not on tidyr1.0)
  first_bs <- head(colnames(gsea_scores_each), n=1)
  last_bs  <- tail(colnames(gsea_scores_each), n=1)
  
  rna_library <- gsub(" ", "_", rna_library)
  rna_library <- stringr::str_to_lower(gsub("-", "", rna_library))
  
  gsea_scores_each_df_tidy <- gsea_scores_each_df %>%
    gather(Kids_First_Biospecimen_ID, gsea_score, !!first_bs : !!last_bs) %>%
    select(Kids_First_Biospecimen_ID, hallmark_name, gsea_score) %>%
    mutate(data_type = rna_library)
  
  gsea_scores_df_tidy <-  bind_rows(gsea_scores_df_tidy , gsea_scores_each_df_tidy)
}

#### Export GSEA scores to TSV --------------------------------------------------------------------
# Specify output path
scores_output_file <- (file.path(results_dir, "gsea_out.tsv"))
# Create output file
write_tsv(gsea_scores_df_tidy, scores_output_file)


