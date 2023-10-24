################################################################################
# 00-pre-process-for-plots.R
# script that takes the original histology file and maps/adds plot_group 
# columns to tsv file that is used in subsequent scripts.
#
# written by Ammar Naqvi
#
# usage: Rscript 00-pre-process-for-plots.R
################################################################################
# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
})

#Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "splicing_index")
input_dir   <- file.path(analysis_dir, "input")

# Load datasets and pre-process
hist_df <- read_tsv(file.path(data_dir,"histologies.tsv"), guess_max = 100000) %>% 
  # filter
  filter(experimental_strategy == "RNA-Seq",
         cohort == "PBTA",
         !is.na(pathology_diagnosis),
         composition != "Derived Cell Line") %>%
  # collapse reported gender to 3 groups
  mutate(reported_gender = case_when(reported_gender == "Not Reported" ~ "Unknown",
                                     TRUE ~ reported_gender),
         # update 7316-3066
         broad_histology = case_when(sample_id == "7316-3066" ~ "Tumor of cranial and paraspinal nerves", 
                                     broad_histology == "Other" ~ "Other tumor",
                                     cancer_group == "Glial-neuronal tumor" ~ "Neuronal and mixed neuronal-glial tumor",
                                     cancer_group %in% c("Cavernoma", "Malignant peripheral nerve sheath tumor") ~ "Benign tumor",
                                     TRUE ~ broad_histology),
         cancer_group = case_when(sample_id == "7316-3066" ~ "Neurofibroma/Plexiform",
                                  grepl("xanthogranuloma", pathology_free_text_diagnosis) & broad_histology == "Histiocytic tumor" & pathology_diagnosis == "Other" ~ "Juvenile xanthogranuloma",
                                  broad_histology == "Choroid plexus tumor" ~ "Choroid plexus tumor",
                                  cancer_group == "Glial-neuronal tumor" ~ "Glial-neuronal tumor NOS",
                                  cancer_group == "Low-grade glioma" ~ "Low-grade glioma/astrocytoma",
                                  cancer_group %in% c("High-grade glioma", "Astrocytoma", "Astroblastoma", "Glioblastoma", "Diffuse hemispheric glioma",
                                                      "Infant-type hemispheric glioma") ~ "High-grade glioma/astrocytoma",
                                  broad_histology == "Ependymal tumor" ~ "Ependymoma",
                                  broad_histology == "Other tumor" ~ "Other tumor",
                                  broad_histology == "Diffuse astrocytic and oligodendroglial tumor" & (is.na(cancer_group) | cancer_group == "Oligodendroglioma") ~ "High-grade glioma/astrocytoma",
                                  cancer_group %in% c("Non-germinomatous germ cell tumor", "Diffuse leptomeningeal glioneuronal tumor", "Malignant peripheral nerve sheath tumor") ~ NA_character_,
                                  broad_histology == "Meningioma" ~ "Meningioma",
                                  cancer_group == "Perineuroma" ~ "Neurofibroma/Plexiform",
                                  is.na(cancer_group) & broad_histology == "Tumor of cranial and paraspinal nerves" ~ "Neurofibroma/Plexiform",
                                  TRUE ~ cancer_group))

# add cancer/plot group mapping file 
map_file <- read_tsv(file.path(input_dir, "plot_mapping.tsv")) 

# add plot mapping to histlogy df
combined_hist_map <- hist_df %>%
  left_join(map_file, by = c("broad_histology", "cancer_group")) %>% 
  write_tsv(file.path(input_dir, "histologies-plot_group.tsv"))