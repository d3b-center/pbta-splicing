################################################################################
# 08-model_and_plot-gsea.R
# written by Ammar Naqvi modified from OpenPBTA
#
# usage: Rscript 08-model_and_plot-gsea.R
################################################################################

## load libraries: 
suppressPackageStartupMessages({
  library("broom")
  library("dplyr")
  library("tidyverse")
  library("vroom")
})

## Magrittr pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1_specific")
util_dir <- file.path(root_dir, "analyses", "CLK1_specific","util")
figures_dir <- file.path(root_dir, "figures")

input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

# This script contains functions used to modeling GSVA scores
source(file.path(util_dir,"hallmark_models.R"))

# source function for theme for plots 
source(file.path(figures_dir, "theme_for_plots.R"))

# Significance testing universal threshold
SIGNIFICANCE_THRESHOLD <- 0.05

######### Define input files
## Metadata file (histologies/clinical data)
metadata_file <- file.path(data_dir, "histologies.tsv")

## GSEA scores
scores_file <- file.path(results_dir, "gsea_out.tsv")

######## Load input files
## histology file
metadata    <- readr::read_tsv(metadata_file, guess_max = 100000) %>% 
  
  ## filters to retrieve samples of interests 
  dplyr::filter(experimental_strategy == "RNA-Seq") %>% 
  dplyr::filter(!is.na(RNA_library)) %>%
  dplyr::filter(cohort == "PBTA") %>%
  dplyr::filter(RNA_library == "stranded") %>% 
  dplyr::filter(CNS_region == 'Midline') %>% 
  dplyr::filter(short_histology == 'HGAT') 

## merged rmats output files
rmats_file <-  file.path(data_dir,"rMATS_merged.single.SE.tsv.gz")
rmats_df <-  vroom(rmats_file, comment = "#",delim="\t") %>%
  
  # select specific samples and extract CLK1 exon 4  
  dplyr::filter(geneSymbol=="CLK1") %>% dplyr::filter(exonStart_0base=="200860124", exonEnd=="200860215") %>%  
  
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


scores_file <- readr::read_tsv(scores_file) 

######## Find out unique RNA library types
rna_library_list <- scores_file %>% pull(data_type) %>% unique()

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

## Merge histology metadata with each set of gsea scores
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
    readr::write_tsv(CLK1_model_results[["anova"]], CLK1_group_anova_outpaths[[i]])
    readr::write_tsv(CLK1_model_results[["tukey"]],  CLK1_group_tukey_outpaths[[i]])
    
    # print results for viewing 
    print(rna_library)
    print(head(CLK1_model_results))
  }
}

## plot significant pathways in heatmap format
sign_pathways_hm <- CLK1_model_results[["tukey"]] %>% filter(significant_tukey_bonf == TRUE) %>% 
                                                      dplyr::select(hallmark_name)

sign_pathways_mat <- metadata_with_gsva %>% 
                     dplyr::select(Kids_First_Biospecimen_ID, hallmark_name,gsea_score) %>% 
                     dplyr::filter(hallmark_name %in% sign_pathways_hm$hallmark_name)

sign_pathways_mat_grouping <- full_join(sign_pathways_mat, metadata) %>% select(Kids_First_Biospecimen_ID, hallmark_name,gsea_score, CLK1_group)

gsea_scores_hm <- ggplot(sign_pathways_mat_grouping, aes(y=hallmark_name, x=Kids_First_Biospecimen_ID, fill=gsea_score)) + 
  geom_raster() + geom_tile() +  scale_fill_gradient(low="white", high="#DC3220") + xlab("Sample") +
  ylab("Pathway") + theme_Publication() + theme(axis.text.x = element_text(size = 10,angle = 45,hjust = .7, vjust = .8),axis.text.y = element_text(size = 7)) 

file_tiff_heatmap_plot = file.path(plots_dir,"heatmap_sign_gsea-kegg.tiff")

# save plot tiff version
tiff(file_tiff_heatmap_plot, height = 2000, width = 2800, res = 300)
print(gsea_scores_hm)
dev.off()


