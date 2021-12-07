# Author: Komal S. Rathi
# Function: Brain cell type proportion analysis using BRETIGEA 

# load libraries
suppressPackageStartupMessages({
  library(ggplot2)
})

# root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analyses_dir <- file.path(root_dir, 'analyses', 'brain_cells_prop')
output_dir <- file.path(analyses_dir, 'output')

# source function to calculate and plot brain cell type proportions
source(file.path(analyses_dir, 'util', 'brain_cells_prop.R'))

p <- brain_cells_prop(expr = file.path(root_dir, 'data', 'tpm_matrix.rds'), 
                   meta = file.path(root_dir, 'data', 'metadata.rds'), 
                   method = "t.test")
ggsave(filename = file.path(output_dir, 'bretigea_plot.pdf'), plot = p, device = "pdf", width = 10, height = 7)


