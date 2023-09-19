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

clusters_df <- vroom(file.path(output_dir,"optimal_clustering",optimal_cluster_tsv),show_col_types = FALSE) %>% 
  rename("Kids_First_Biospecimen_ID" = sample)

histology_w_clusters_df <- vroom(file.path(data_dir,histology_file),show_col_types = FALSE) %>% 
  full_join(clusters_df, by="Kids_First_Biospecimen_ID" ) %>% 
  filter(!is.na(as.numeric(cluster_assigned)))

color_df <- vroom(file.path(root_dir,"palettes/histology_label_color_table.tsv"), delim="\t", col_names = TRUE, trim_ws = TRUE, show_col_types = FALSE)
cols <- as.character(color_df$cancer_group_hex_codes)
names(cols) <- as.character(color_df$short_histology)

histology_w_clusters_df$short_histology <- gsub("^HGAT", "HGG", histology_w_clusters_df$short_histology)
histology_w_clusters_df$short_histology <- gsub("^LGAT", "LGG", histology_w_clusters_df$short_histology)

names(cols) <- gsub("^HGAT", "HGG",names(cols))
names(cols) <- gsub("^LGAT", "LGG",names(cols))



tiff(file.path(plots_dir, "cluster_membership.tiff"), height = 1200, width = 2200, units = "px", res = 300)
ggplot(histology_w_clusters_df, aes(fill=short_histology, x=cluster_assigned)) +
  geom_bar(stat="count", position="stack") + 
  xlab("Cluster") + ylab("Frequency") +
  scale_fill_manual("Histology", values = cols) + theme_Publication()
dev.off()
