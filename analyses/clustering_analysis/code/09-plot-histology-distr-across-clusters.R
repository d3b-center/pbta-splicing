################################################################################
# 09-plot-histology-distr-across-clusters.R
# Plot the distribution of histologies across clusters
#
# Author: Ammar Naqvi, Jo Lynne Rokita
################################################################################

## libraries 
suppressPackageStartupMessages({
  library("tidyverse")
  library("vroom")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`


## directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "clustering_analysis")

plots_dir <- file.path(analysis_dir, "plots")
input_dir <- file.path(analysis_dir, "input")
output_dir <- file.path(analysis_dir, "output")

plots_dir <- file.path(analysis_dir, "plots")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

## theme for all plots
# source function for theme for plots survival
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## filepaths 
optimal_cluster_tsv <- file.path(output_dir, "optimal_clustering", "ccp_optimal_clusters.tsv")
cluster_membership_tsv <-  file.path(output_dir, "cluster_members_by_cancer_group_subtype.tsv")

# read in palette and cluster file
cluster_df <- read_tsv(optimal_cluster_tsv) %>%
  dplyr::rename(Kids_First_Biospecimen_ID = sample)

histologies_df <- read_tsv(file.path(root_dir,"analyses", "cohort_summary", "results", "histologies-plot-group.tsv"), guess_max = 100000) %>%
  dplyr::select(Kids_First_Biospecimen_ID, broad_histology, cancer_group, plot_group, plot_group_hex, molecular_subtype) %>%
  dplyr::filter(Kids_First_Biospecimen_ID %in% cluster_df$Kids_First_Biospecimen_ID) %>%
  unique() %>%
  dplyr::right_join(cluster_df)

color_df <- histologies_df %>%
  dplyr::select(plot_group_hex, plot_group) %>%
  dplyr::filter(!is.na(plot_group)) %>%
  unique()
cols <- as.character(color_df$plot_group_hex)
names(cols) <- as.character(color_df$plot_group)

# create plot
pdf(file.path(plots_dir, "cluster_membership.pdf"), height = 4, width = 7)
ggplot(histologies_df, aes(fill=plot_group, x= factor(cluster_assigned))) +
  geom_bar(stat="count", position="stack") + 
  xlab("Cluster") + ylab("Frequency") +
  scale_fill_manual("Histology", values = cols) + 
  theme_Publication()
dev.off()

## stratify by molecular subtype
pdf(file.path(plots_dir, "cluster_membership-MB.pdf"), height = 4, width = 7)
histologies_MB_plot <-  histologies_df %>% filter(plot_group=='Medulloblastoma') %>% 
  ggplot(aes(fill=molecular_subtype, x= factor(cluster_assigned))) +
  geom_bar(stat="count", position="stack") + 
  xlab("Cluster") + ylab("Frequency") +
  #scale_fill_manual("Histology", values = cols) + 
  theme_Publication()
histologies_MB_plot
dev.off()

pdf(file.path(plots_dir, "cluster_membership-ATRT.pdf"), height = 4, width = 7)
histologies_ATRT_plot <-  histologies_df %>% filter(plot_group=='Atypical Teratoid Rhabdoid Tumor') %>% 
  ggplot(aes(fill=molecular_subtype, x= factor(cluster_assigned))) +
  geom_bar(stat="count", position="stack") + 
  xlab("Cluster") + ylab("Frequency") +
  #scale_fill_manual("Histology", values = cols) + 
  theme_Publication()
histologies_ATRT_plot
dev.off()

pdf(file.path(plots_dir, "cluster_membership-EPN.pdf"), height = 4, width = 7)
histologies_EPN_plot <-  histologies_df %>% filter(plot_group=='Ependymoma') %>% 
  ggplot(aes(fill=molecular_subtype, x= factor(cluster_assigned))) +
  geom_bar(stat="count", position="stack") + 
  xlab("Cluster") + ylab("Frequency") +
  #scale_fill_manual("Histology", values = cols) + 
  theme_Publication()
histologies_EPN_plot
dev.off()

pdf(file.path(plots_dir, "cluster_membership-DMG.pdf"), height = 4, width = 7)
histologies_DMG_plot <-  histologies_df %>% filter(plot_group=='DIPG or DMG') %>% 
  ggplot(aes(fill=molecular_subtype, x= factor(cluster_assigned))) +
  geom_bar(stat="count", position="stack") + 
  xlab("Cluster") + ylab("Frequency") +
  #scale_fill_manual("Histology", values = cols) + 
  theme_Publication()
histologies_DMG_plot
dev.off()

pdf(file.path(plots_dir, "cluster_membership-OHGG.pdf"), height = 4, width = 7)
histologies_HGG_plot <-  histologies_df %>% filter(plot_group=='Other high-grade glioma') %>% 
  ggplot(aes(fill=molecular_subtype, x= factor(cluster_assigned))) +
  geom_bar(stat="count", position="stack") + 
  xlab("Cluster") + ylab("Frequency") +
  #scale_fill_manual("Histology", values = cols) + 
  theme_Publication()
histologies_HGG_plot
dev.off()

pdf(file.path(plots_dir, "cluster_membership-LGG.pdf"), height = 4, width = 14) ## too many subtypes! 
histologies_LGG_plot <-  histologies_df %>% filter(plot_group=='Low-grade glioma') %>% 
  ggplot(aes(fill=molecular_subtype, x= factor(cluster_assigned))) +
  geom_bar(stat="count", position="stack") + 
  xlab("Cluster") + ylab("Frequency") +
  #scale_fill_manual("Histology", values = cols) + 
  theme_Publication()
histologies_LGG_plot
dev.off()

# write out bs id with cancer group, subtype, and cluster
histologies_df %>%
  dplyr::select(Kids_First_Biospecimen_ID, cancer_group, plot_group, molecular_subtype, cluster_assigned) %>%
  write_tsv(cluster_membership_tsv)
