suppressPackageStartupMessages({
  library(tidyverse)
  library("vroom")
})

#Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

##set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "histology-specific-splicing")
results_dir <- file.path(analysis_dir, "results")

plots_dir <- file.path(analysis_dir, "plots")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

histology_file <- file.path(data_dir,"histologies.tsv")
histology_df <- vroom(histology_file)

summ_dmg_file <- "/Users/naqvia/Desktop/splicing-based_neoepitope_discovery/DMG_subset/summary_table.tsv"
summ_dmg_df <-  vroom(summ_dmg_file, comment = "#",delim="\t", col_names = TRUE, trim_ws = TRUE, show_col_types = FALSE)

H3_wt_df <- summ_dmg_df %>% filter(H3_WT=="Y") %>% mutate(Kids_First_Participant_ID= paste0("PT_",str_split(match_id, "_", simplify = TRUE)[,2]))
H3_1_k28m_df <- summ_dmg_df %>% filter(H3.1_K28M_sum=="Y") %>% mutate(Kids_First_Participant_ID= paste0("PT_",str_split(match_id, "_", simplify = TRUE)[,2]))
H3_3_k28m_df <- summ_dmg_df %>% filter(H3.3_K28M_sum=="Y") %>% mutate(Kids_First_Participant_ID= paste0("PT_",str_split(match_id, "_", simplify = TRUE)[,2]))

## rename cancer group for downstream scripts
histology_wt_df    <- histology_df %>% inner_join(H3_wt_df, by="Kids_First_Participant_ID") %>% mutate(cancer_group="H3 WT") 
histology_H31_k28m <- histology_df %>% inner_join(H3_1_k28m_df, by="Kids_First_Participant_ID") %>% mutate(cancer_group="H3.1 K28") 
histology_H33_k28m <- histology_df %>% inner_join(H3_3_k28m_df, by="Kids_First_Participant_ID")  %>% mutate(cancer_group="H3.3 K28") 

histology_modified <- rbind(histology_wt_df,histology_H31_k28m,histology_H33_k28m) 
write_tsv(histology_modified, file=file.path(results_dir,"histologies-dmg.tsv"))
