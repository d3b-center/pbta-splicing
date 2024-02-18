# Author: Komal S. Rathi, Jo Lynne Rokita
# Function: Get ccp clustering plot for a specific combination of distance + algorithm + % variable genes 

# load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
#  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
})

# source function to perform dip test
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "clustering_analysis")
output_dir <- file.path(analysis_dir, "plots")

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
  as.data.frame() %>% 
  filter(if_any(everything(), ~ !is.na(.)))

colnames(CC_consensus_mat) <- ccp_output[[n_cluster]]$consensusClass %>% names()
rownames(CC_consensus_mat) <- ccp_output[[n_cluster]]$consensusClass %>% names()

# convert to matrix
CC_consensus_mat <- as.matrix(CC_consensus_mat)

# get CCP tree
CC_tree <- ccp_output[[n_cluster]]$consensusTree 

# get cluster assignments for samples
CC_class <- ccp_output[[n_cluster]]$consensusClass
CC_class <- data.frame(cluster_class = CC_class) 

stopifnot(identical(colnames(CC_consensus_mat), rownames(CC_class)))


# add short histology and color code to clusters
CC_annot <- CC_class %>%
  rownames_to_column("Kids_First_Biospecimen_ID") %>%
  inner_join(input_clin, by = "Kids_First_Biospecimen_ID") %>%
  dplyr::select(Kids_First_Biospecimen_ID, cluster_class, plot_group, plot_group_hex) %>%
  # make sure cluster legend maintains its order
  dplyr::mutate(cluster_class = factor(cluster_class, levels = as.character(1:n_cluster))) %>%
  column_to_rownames("Kids_First_Biospecimen_ID") 
  

# rename annotation columns
CC_annot <- CC_annot %>%
  dplyr::rename("Histology" = "plot_group",
                "Cluster" = "cluster_class")

# create annotation for cluster class
# colors for marking different clusters
thisPal <- c("#bd18ea","#FB9A99","#FDBF6F","#B2DF8A","#B15928",
             "#E31A1C","#33A02C","#FF7F00","#A6CEE3","#6A3D9A",
             "#CAB2D6","#1F78B4","#FFFF99","#2ef4ca","#f4cced",
             "#f4cc03", #lightorange
             "#05188a", #navy,
             "#e5a25a", #light brown
             "#06f106", #bright green
             "#85848f", #med gray
             "#000000", #black
             "#076f25", #dark green
             "#93cd7f",#lime green
             "#4d0776" #dark purple
)

l <- thisPal[1:n_cluster]
names(l) <- as.character(1:n_cluster)
mycolors <- list()
mycolors[['Cluster']] <- l 

# create annotation for short histology
short_histology_palettes <- CC_annot %>%
  dplyr::select(Histology,plot_group_hex) %>%
  unique()

mycolors[['Histology']] <- short_histology_palettes$plot_group_hex
names(mycolors[['Histology']]) <- short_histology_palettes$Histology

# remove colors from annotation table
CC_annot$plot_group_hex <- NULL

#CC_annot <- CC_annot #%>% 
#  filter(if_any(everything(), ~ !is.na(.)))

# Reorder CC_annot to match the column order of CC_consensus_mat
CC_annot <- CC_annot[colnames(CC_consensus_mat), ]

# Assuming output_dir and prefix variables are defined
# Define the dimensions in inches for 300 DPI resolution
# For example, for a 10x10 inch image at 300 DPI:
width_in_inches <- 10
height_in_inches <- 7
res_dpi <- 300

# Convert dimensions to pixels
width_in_pixels <- width_in_inches * res_dpi
height_in_pixels <- height_in_inches * res_dpi

# Open a TIFF device
tiff(file.path(output_dir, paste0(prefix, '_ccp_heatmap.tiff')), 
     res = res_dpi, 
     width = width_in_pixels, 
     height = height_in_pixels, 
     units = "px")
# print heatmap
print(pheatmap(CC_consensus_mat,
         color = colorRampPalette(brewer.pal(n = 9, name = "Blues"))(9),
         fontsize = 12,
         main = "",
         name = "Similarity",
         show_colnames = F, 
         show_rownames = F, 
         annotation_col = CC_annot,
         annotation_colors = mycolors, 
         cluster_cols = CC_tree, 
         cluster_rows = CC_tree, 
         treeheight_row = 0))
dev.off()

