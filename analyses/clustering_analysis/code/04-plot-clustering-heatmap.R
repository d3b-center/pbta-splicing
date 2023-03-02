# Author: Komal S. Rathi
# Function: Get ccp clustering plot for a specific combination of distance + algorithm + % variable genes 

# load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
  library(pheatmap)
  library(RColorBrewer)
})

# source function to perform dip test
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "clustering_analysis")
output_dir <- file.path(analysis_dir, "output", "ccp_output")

# parse parameters     
option_list <- list(
  make_option(c("--ccp_output"), type = "character",
              help = "CCP output from running ConsensusClusterPlus (.rds)"),
  make_option(c("--input_clin"), type = "character",
              help = "histologies (.tsv)"),
  make_option(c("--n_cluster"), type = "character",
              help = "number of clusters to investigate"),
  make_option(c("--prefix"), type = "character",
              help = "prefix of output files")
)
opt <- parse_args(OptionParser(option_list = option_list))
ccp_output <- opt$ccp_output %>% readRDS()
input_clin <- opt$input_clin %>% readr::read_tsv()
n_cluster <- as.numeric(opt$n_cluster)
prefix <- opt$prefix

# get CCP matrix
CC_consensus_mat <- ccp_output[[n_cluster]]$consensusMatrix %>%
  as.data.frame()
colnames(CC_consensus_mat) <- ccp_output[[n_cluster]]$consensusClass %>% names()

# get CCP tree
CC_tree <- ccp_output[[n_cluster]]$consensusTree

# get cluster assignments for samples
CC_class <- ccp_output[[n_cluster]]$consensusClass
CC_class <- data.frame(cluster_class = CC_class)

# color palette for short histology
palettes_dir <- "../../palettes/"
palette_file <- file.path(palettes_dir, "short_histology_color_palette.tsv") %>% read_tsv()

# add short histology and color code to clusters
CC_annot <- CC_class %>%
  rownames_to_column("Kids_First_Biospecimen_ID") %>%
  inner_join(input_clin, by = "Kids_First_Biospecimen_ID") %>%
  inner_join(palette_file, by = "short_histology") %>%
  dplyr::select(Kids_First_Biospecimen_ID, cluster_class, plot_group_display, hex_code) %>%
  dplyr::mutate(cluster_class = as.character(cluster_class)) %>%
  column_to_rownames("Kids_First_Biospecimen_ID")

# rename annotation columns
CC_annot <- CC_annot %>%
  dplyr::rename("Short Histology" = "plot_group_display",
                "Cluster Class" = "cluster_class")

# create annotation for cluster class
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
l <- gg_color_hue(n_cluster)
names(l) <- as.character(1:n_cluster)
mycolors <- list()
mycolors[['Cluster Class']] <- l 

# create annotation for short histology
short_histology_palettes <- CC_annot %>%
  dplyr::select(`Short Histology`, hex_code) %>%
  unique()
mycolors[['Short Histology']] <- short_histology_palettes$hex_code
names(mycolors[['Short Histology']]) <- short_histology_palettes$`Short Histology`

# remove colors from annotation table
CC_annot$hex_code <- NULL

# create heatmap
pheatmap(CC_consensus_mat,
         color = colorRampPalette(brewer.pal(n = 9, name = "Blues"))(9),
         fontsize = 10,
         main = "Consensus Matrix",
         show_colnames = F, 
         show_rownames = F,
         annotation_col = CC_annot, 
         annotation_colors = mycolors, 
         cluster_cols = CC_tree, 
         filename = file.path(output_dir, paste0(prefix, '_ccp_heatmap.tiff')), 
         width = 10, height = 8)
