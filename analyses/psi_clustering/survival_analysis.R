################################################################################
# survival_analysis.R
# script that takes cluster members and does KM curves
# written by Ammar Naqvi
#
# usage: Rscript survival_analysis.R
################################################################################
library(dplyr)
library(ggthemes) 
# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

file = "/Users/naqvia/Desktop/AS-DMG/analyses/psi_clustering/input/CC_groups_remDups.v2.txt"
file = "/Users/naqvia/Desktop/pbta-splicing/analyses/psi_clustering/results/Clust_HGGs.txt"
cluster_mem  = read.delim(file, sep = "\t", header=TRUE)

# source function to compute survival
source(file.path("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/util", "survival_models.R"))
`%>%` <- dplyr::`%>%`

# read pbta-histology file
metadata <- readr::read_tsv("/Users/naqvia/Desktop/pbta-splicing/analyses/psi_clustering/input/pbta-histologies.tsv")

cluster_mem <- cluster_mem %>%
  dplyr::left_join(metadata, by = "Kids_First_Biospecimen_ID")

# Kaplan-Meier for all clusters
kap_fit <- survival_analysis(cluster_mem,
                             ind_var = "Cluster",
                             test = "kap.meier",
                             metadata_sample_col = "Kids_First_Biospecimen_ID")

survminer::ggsurvplot(kap_fit$model,
                      pval = TRUE,
                      data = kap_fit$original_data,
                      risk.table = TRUE,
                      break.time.by = 500,
                      #ggtheme = theme_Publication(),
                      risk.table.y.text.col = TRUE,
                      risk.table.y.text = FALSE)

## for SI (HGGs-- low vs high)
file = "/Users/naqvia/Desktop/pbta-splicing/analyses/splicing_index/results/high_psi.hgg.tmp"
si_mem  = read.delim(file, sep = "\t", header=TRUE)

# source function to compute survival
source(file.path("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/util", "survival_models.R"))
`%>%` <- dplyr::`%>%`

# read pbta-histology file
metadata <- readr::read_tsv("/Users/naqvia/Desktop/pbta-splicing/analyses/psi_clustering/input/pbta-histologies.tsv")

si_mem <- si_mem %>%
  dplyr::left_join(metadata, by = "Kids_First_Biospecimen_ID")

# Kaplan-Meier for all clusters
kap_fit <- survival_analysis(si_mem,
                             ind_var = "SI",
                             test = "kap.meier",
                             metadata_sample_col = "Kids_First_Biospecimen_ID")

survminer::ggsurvplot(kap_fit$model,
                      pval = TRUE,
                      data = kap_fit$original_data,
                      risk.table = TRUE,
                      break.time.by = 500,
                      ggtheme = theme_Publication(),
                      risk.table.y.text.col = TRUE,
                      risk.table.y.text = FALSE)

## separate out based on hist
file= "/Users/naqvia/Desktop/pbta-splicing/analyses/psi_clustering/input/CC_groups_remDups.v2.txt"
cluster_mem  = read.delim(file, sep = "\t", header=TRUE)

hist_file = "/Users/naqvia/Desktop/pbta-splicing/analyses/psi_clustering/input/v19_plus_20210311_pnoc.tsv"
hist = read.delim(hist_file, sep = "\t", header=TRUE)

right_join_w_CC = right_join(hist, cluster_mem, by = "Kids_First_Biospecimen_ID")

lgat_cc = filter(right_join_w_CC, short_histology=="LGAT")
med_cc = filter(right_join_w_CC, short_histology=="Medulloblastoma")
hgat_cc = filter(right_join_w_CC, short_histology=="HGAT")
gang_cc = filter(right_join_w_CC, short_histology=="Ganglioglioma")
atrt_cc = filter(right_join_w_CC, short_histology=="ATRT")

# Kaplan-Meier for ATRT clusters
kap_fit <- survival_analysis(atrt_cc,
                             ind_var = "CC",
                             test = "kap.meier",
                             metadata_sample_col = "Kids_First_Biospecimen_ID")

survminer::ggsurvplot(kap_fit$model,
                      pval = TRUE,
                      data = kap_fit$original_data,
                      risk.table = TRUE,
                      break.time.by = 500,
                      #ggtheme = theme_Publication(),
                      risk.table.y.text.col = TRUE,
                      risk.table.y.text = FALSE)

# Kaplan-Meier for gangl clusters
kap_fit <- survival_analysis(gang_cc,
                             ind_var = "CC",
                             test = "kap.meier",
                             metadata_sample_col = "Kids_First_Biospecimen_ID")

survminer::ggsurvplot(kap_fit$model,
                      pval = TRUE,
                      data = kap_fit$original_data,
                      risk.table = TRUE,
                      break.time.by = 500,
                      #ggtheme = theme_Publication(),
                      risk.table.y.text.col = TRUE,
                      risk.table.y.text = FALSE)

# Kaplan-Meier for hgat clusters
kap_fit <- survival_analysis(hgat_cc,
                             ind_var = "CC",
                             test = "kap.meier",
                             metadata_sample_col = "Kids_First_Biospecimen_ID")

survminer::ggsurvplot(kap_fit$model,
                      pval = TRUE,
                      data = kap_fit$original_data,
                      risk.table = TRUE,
                      break.time.by = 500,
                      #ggtheme = theme_Publication(),
                      risk.table.y.text.col = TRUE,
                      risk.table.y.text = FALSE)

# Kaplan-Meier for lgat clusters
kap_fit <- survival_analysis(lgat_cc,
                             ind_var = "CC",
                             test = "kap.meier",
                             metadata_sample_col = "Kids_First_Biospecimen_ID")

survminer::ggsurvplot(kap_fit$model,
                      pval = TRUE,
                      data = kap_fit$original_data,
                      risk.table = TRUE,
                      break.time.by = 500,
                      #ggtheme = theme_Publication(),
                      risk.table.y.text.col = TRUE,
                      risk.table.y.text = FALSE)

# Kaplan-Meier for medullo clusters 
kap_fit <- survival_analysis(med_cc,
                             ind_var = "CC",
                             test = "kap.meier",
                             metadata_sample_col = "Kids_First_Biospecimen_ID")

survminer::ggsurvplot(kap_fit$model,
                      pval = TRUE,
                      data = kap_fit$original_data,
                      risk.table = TRUE,
                      break.time.by = 500,
                      #ggtheme = theme_Publication(),
                      risk.table.y.text.col = TRUE,
                      risk.table.y.text = FALSE)
