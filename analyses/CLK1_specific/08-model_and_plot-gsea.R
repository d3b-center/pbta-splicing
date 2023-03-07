################################################################################
# 08-model_and_plot-gsea.R
# Author: Ammar Naqvi (script adapted OpenPBTA)
#
# usage: Rscript 08-model_and_plot-gsea.R
################################################################################

## Load libraries: 
## load libraries
suppressPackageStartupMessages({
  # Library for data manipulation
  library("tidyverse")
  library("vroom")
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


# This script contains functions used for modeling GSVA scores
source(file.path(analysis_dir, "util", "hallmark_models.R"))

# Source function for plots theme
source(file.path(root_dir, "figures", "theme_for_plots.R"))

# Significance testing universal threshold
SIGNIFICANCE_THRESHOLD <- 0.05


# Load input files
## Load metadata file path
metadata_file <- file.path(data_dir, "histologies.tsv")
# Specify GSEA scores file path
scores_file <- file.path(results_dir, "gsea_out.tsv")

# # Load metadata
# metadata <- readr::read_tsv(metadata_file, guess_max = 100000) %>% 
#   # Select only RNA-Seq samples
#   filter(experimental_strategy == "RNA-Seq") %>% 
#   # Filter for module-specific samples that was previously used in CLK1 splicing/comparisons
#   filter(Kids_First_Biospecimen_ID %in% c('BS_Q13FQ8FV', 'BS_ZV1P6W9C', 'BS_WH8G4VFB', 
#                                           'BS_NNPEC7W1', 'BS_PZVHMSYN', 'BS_DRY58DTF', 
#                                           'BS_GXTFW99H', 'BS_E60JZ9Z3', 'BS_9CA93S6D')) %>% 
#   ## Define CLK1 groups as "High" or "Low"
#   mutate(CLK1_group = ifelse((Kids_First_Biospecimen_ID   == 'BS_Q13FQ8FV') |
#                              (Kids_First_Biospecimen_ID   == 'BS_ZV1P6W9C') |
#                              (Kids_First_Biospecimen_ID   == 'BS_WH8G4VFB') | 
#                              (Kids_First_Biospecimen_ID   == 'BS_NNPEC7W1') |
#                              (Kids_First_Biospecimen_ID   == 'BS_PZVHMSYN'), 
#                              "High", "Low"))





######## Load input files
## histology file
metadata <- readr::read_tsv(metadata_file, guess_max = 100000) %>% 
  # Only include RNA-Seq, PBTA cohorts, stranded, Midline HGGs samples
  dplyr::filter(experimental_strategy == "RNA-Seq", 
                cohort == "PBTA",
                CNS_region == 'Midline', 
                short_histology == 'HGAT') 

## Merged rmats output files
rmats_file <-  file.path(data_dir,"rMATS_merged.single.SE.tsv.gz")
rmats_df <-  vroom(rmats_file, comment = "#",delim="\t") %>%
  # select specific samples and extract CLK1 exon 4  
  dplyr::filter(geneSymbol=="CLK1") %>% 
  dplyr::filter(exonStart_0base=="200860124", exonEnd=="200860215") %>%  
  # join on metadata from above
  inner_join(metadata, by=c('sample'='Kids_First_Biospecimen_ID')) %>%  
  dplyr::select(sample, geneSymbol, IncLevel1) 

## compute quantiles to define high vs low Exon 4 PSI groups
quartiles_psi <- quantile(rmats_df$IncLevel1, probs=c(.25, .75), na.rm = FALSE)
IQR_psi <- IQR(rmats_df$IncLevel1)
lower_psi <- quartiles_psi[1] 
upper_psi <- quartiles_psi[2] 

rmats_subset_samples_lowEx4  <- rmats_df %>% dplyr::filter(rmats_df$IncLevel1 < lower_psi) %>% rename(Kids_First_Biospecimen_ID = sample)
rmats_subset_samples_highEx4 <- rmats_df %>% dplyr::filter(rmats_df$IncLevel1 > upper_psi) %>%  rename(Kids_First_Biospecimen_ID = sample)
rmats_subset_samples_Ex4 <- rbind(rmats_subset_samples_lowEx4,rmats_subset_samples_highEx4)

metadata <- inner_join(metadata,rmats_subset_samples_Ex4, by="Kids_First_Biospecimen_ID") %>%   
            mutate(CLK1_group = if_else(IncLevel1 < lower_psi, "Low", "High"))
write.csv(metadata, "CLK1_group_status.txt")

# Load GSEA scores
scores_file <- readr::read_tsv(scores_file) 

######## Find out unique RNA library types
rna_library_list <- scores_file %>% 
  pull(data_type) %>% 
  unique()

CLK1_group_tukey_outpaths <- lapply(rna_library_list, function(x){
  x<-gsub(" ", "_", x)
  x<-stringr::str_to_lower(gsub("-", "", x))
  file.path(results_dir, paste0("gsva_tukey_hm_", x, "_CLK1_group"))
})

CLK1_group_anova_outpaths <- lapply(rna_library_list, function(x){
  x<-gsub(" ", "_", x)
  x<-stringr::str_to_lower(gsub("-", "", x))
  file.path(results_dir, paste0("gsva_anova_hm_", x, "_CLK1_group.tsv"))
})

### ANOVA and Tukey analysis of GSVA scores

### Merge histology metadata with each set of gsea scores
## First, prepare the data for modeling:
metadata_with_gsva <- metadata %>%
  inner_join(scores_file, by = "Kids_First_Biospecimen_ID") 

## now model
for(i in 1:length(rna_library_list)){
  
  rna_library = rna_library_list[i]
  
  # find out the number of cancer group with this RNA library
  CLK1_group_n <- metadata_with_gsva %>%
    filter(data_type == rna_library) %>%
    pull(CLK1_group) %>% 
    unique() %>% length()
  
  # anova can only be run on factors with >=2 levels, so to avoid error, we give a if statement 
  if(CLK1_group_n>=2){
    CLK1_model_results <- gsva_anova_tukey(metadata_with_gsva, CLK1_group, rna_library, SIGNIFICANCE_THRESHOLD)
    
    # write out results
    readr::write_tsv(CLK1_model_results[["anova"]], cancer_group_anova_outpaths[[i]])
    readr::write_tsv(CLK1_model_results[["tukey"]],  cancer_group_tukey_outpaths[[i]])
    
    # print results for viewing 
    print(rna_library)
    print(head(CLK1_model_results))
  }
}

## Plot significant pathways in heatmap format
sign_pathways_hm <- CLK1_model_results[["tukey"]] %>% 
  filter(significant_tukey == TRUE) %>% 
  select(hallmark_name)

# First make significant pathways dataframe
sign_pathways_mat <- metadata_with_gsva %>% 
                     select(Kids_First_Biospecimen_ID, hallmark_name, gsea_score) %>% 
                     filter(hallmark_name %in% sign_pathways_hm$hallmark_name)

# Plot GSEA Heatmap
gsea_scores_hm <- ggplot(sign_pathways_mat, aes(y=hallmark_name, x=Kids_First_Biospecimen_ID, fill=gsea_score)) + 
  geom_raster() + 
  geom_tile() +  
  scale_fill_gradient(low="white", high="#DC3220") + 
  xlab("Sample") +
  ylab("Pathway") + 
  theme_Publication() + 
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = .7, vjust = .8), 
        axis.text.y = element_text(size = 10)) 
# View Heatmap
print(gsea_scores_hm)

# Path to save file as
file_tiff_heatmap_plot <- file.path(plots_dir,"heatmap_sign_gsea.tiff")

# Save plot as tiff
tiff(file_tiff_heatmap_plot, height = 1600, width = 3200, res = 300)

# Close the plotting device
dev.off()
