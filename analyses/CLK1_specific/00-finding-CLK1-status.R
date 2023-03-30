################################################################################
# 00-finding-CLK1-status.R
# Author: Ammar Naqvi
#
# usage: Rscript 00-finding-CLK1-status.R
################################################################################

## Load libraries: 
suppressPackageStartupMessages({
  # Library for data manipulation
  library("tidyverse")
  library("vroom")
})

## Set directories
# Input directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data", "v2")
analysis_dir <- file.path(root_dir, "analyses", "CLK1_specific")
input_dir   <- file.path(analysis_dir, "input")
# Output directories
results_dir <- file.path(analysis_dir, "results")


# Load input files
## Load metadata file path
metadata_file <- file.path(data_dir, "histologies.tsv")

## Load metadata file
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

# Compute quantiles to define high vs low Exon 4 PSI groups
quartiles_psi <- quantile(rmats_df$IncLevel1, probs=c(.25, .75), na.rm = FALSE)
IQR_psi <- IQR(rmats_df$IncLevel1)
lower_psi <- quartiles_psi[1] 
upper_psi <- quartiles_psi[2] 

# Determine Low CLK1 status
rmats_subset_samples_lowEx4  <- rmats_df %>% 
  dplyr::filter(rmats_df$IncLevel1 < lower_psi) %>% 
  rename(Kids_First_Biospecimen_ID = sample)
# Determine High CLK1 status
rmats_subset_samples_highEx4 <- rmats_df %>% 
  dplyr::filter(rmats_df$IncLevel1 > upper_psi) %>%  
  rename(Kids_First_Biospecimen_ID = sample)
# Combine High and Low CLK1 status
rmats_subset_samples_Ex4 <- rbind(rmats_subset_samples_lowEx4,rmats_subset_samples_highEx4)

# Create final metadata with CLK1 status
metadata <- inner_join(metadata,rmats_subset_samples_Ex4, by="Kids_First_Biospecimen_ID") %>%   
  mutate(CLK1_group = if_else(IncLevel1 < lower_psi, "Low", "High")) %>%
  select(Kids_First_Biospecimen_ID, CNS_region, short_histology, geneSymbol, IncLevel1, CLK1_group)

# Export the CLK1 status metadata
readr::write_tsv(metadata, file.path(results_dir, "CLK1_group_status.txt"))


