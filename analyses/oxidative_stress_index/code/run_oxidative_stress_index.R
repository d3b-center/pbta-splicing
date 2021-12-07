# Author: Komal S. Rathi
# GSVA using Oxidative Stress genes

# load libraries
suppressPackageStartupMessages({
  library(ggplot2)
})

# root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analyses_dir <- file.path(root_dir, 'analyses', 'oxidative_stress_index')
output_dir <- file.path(analyses_dir, 'output')

# source function to calculate and plot oxidative stress index
source(file.path(analyses_dir, 'util', 'oxidative_stress_index.R'))

p <- oxidative_stress_index(expr = file.path(root_dir, 'data', 'tpm_matrix.rds'), 
                            meta = file.path(root_dir, 'data', 'metadata.rds'), 
                            method = "t.test") 
ggsave(filename = file.path(output_dir, 'oxidative_stress_index_plot.pdf'), device = "pdf", plot = p,width = 10, height = 7)
