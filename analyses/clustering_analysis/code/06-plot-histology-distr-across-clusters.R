################################################################################
# 06-plot-histology-distr-across-clusters.R
# Plot the distribution of histologies across clusters
#
# Author: Ammar Naqvi 
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

## input 
optimal_cluster_tsv = "ccp_optimal_clusters.tsv"
histology_file = "histologies.tsv"

color_df <- vroom(file.path(root_dir,"palettes/short_histology_color_palette.tsv"), delim="\t", col_names = TRUE, trim_ws = TRUE, show_col_types = FALSE)
cols <- as.character(color_df$hex_code)
names(cols) <- as.character(color_df$short_histology)

histology_w_clusters_df$short_histology <- gsub("^HGAT", "HGG", histology_w_clusters_df$short_histology)
histology_w_clusters_df$short_histology <- gsub("^LGAT", "LGG", histology_w_clusters_df$short_histology)

names(cols) <- gsub("^HGAT", "HGG",names(cols))
names(cols) <- gsub("^LGAT", "LGG",names(cols))

pdf(file.path(plots_dir, "cluster_membership.pdf"), height = 4, width = 8)
ggplot(histology_w_clusters_df, aes(fill=short_histology, x= factor(cluster_assigned))) +
  geom_bar(stat="count", position="stack") + 
  xlab("Cluster") + ylab("Frequency") +
  scale_fill_manual("Histology", values = cols,labels=plot_labels[['short_histology']]) + 
  theme_Publication()
dev.off()
