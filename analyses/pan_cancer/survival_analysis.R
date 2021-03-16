################################################################################
# survival_analysis.R
# script that takes cluster members and does KM curves
# written by Ammar Naqvi
#
# usage: Rscript survival_analysis.R
################################################################################
library(dplyr)

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

file = "/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/results/CC_groups.remDups.txt"
cluster_mem  = read.delim(file, sep = "\t", header=TRUE)

# source function to compute survival
source(file.path("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/util", "survival_models.R"))
`%>%` <- dplyr::`%>%`

# set up directories
setwd("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/results")
input_dir <- "/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer" # you can set this to the directory which contains this code


data_dir  <- file.path("/Users/naqvia/Desktop/AS-DMG/data") # contains pbta-histologies file and all input .dat files
plots_dir <- file.path("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/plots_survival_remDups")  # output directory for plot

# create output directories
dir.create(plots_dir, recursive = T, showWarnings = F)

# read pbta-histology file
metadata <- readr::read_tsv(file.path(data_dir, "pbta-histologies.tsv"))
#metadata <- metadata %>%
#  dplyr::filter(broad_histology == "Diffuse astrocytic and oligodendroglial tumor",
#                composition == "Solid Tissue")


cluster_mem <- cluster_mem %>%
  dplyr::left_join(metadata, by = "Kids_First_Biospecimen_ID")


# Patients with multiple different molecular subtype might affect downstream analysis for survival so using random selection of 1-1 sample-patient matches
unique_match <- cluster_mem %>%
  dplyr::group_by(Kids_First_Participant_ID) %>%
  dplyr::summarize(Kids_First_Biospecimen_ID = sample(Kids_First_Biospecimen_ID, 1))

# select only independent primary WGS bs_ids to get overall survival
metadata_sub <- metadata %>%
  dplyr::filter(Kids_First_Biospecimen_ID %in% unique_match$Kids_First_Biospecimen_ID) %>%
  dplyr::left_join(cluster_mem[,c("sample_id","Cluster")], suffix = c("_DNA","_RNA")) %>%
  dplyr::filter(!is.na(Cluster))

# Kaplan-Meier for all HGG/DMG subtypes
kap_fit <- survival_analysis(metadata_sub,
                             ind_var = "Cluster",
                             test = "kap.meier",
                             metadata_sample_col = "Kids_First_Biospecimen_ID")

surv_pvalue(kap_fit$model, data = kap_fit$original_data)
spval<- surv_pvalue(kap_fit$model, data = kap_fit$original_data)

survminer::ggsurvplot(kap_fit$model,
                      pval = TRUE,
                      data = kap_fit$original_data,
                      risk.table = TRUE,
                      break.time.by = 500,
                      ggtheme = theme_Publication(),
                      risk.table.y.text.col = TRUE,
                      risk.table.y.text = FALSE)







