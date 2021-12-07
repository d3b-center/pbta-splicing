# Author: Komal S. Rathi
# GSVA using DNA Repair gene sets

# load libraries
suppressPackageStartupMessages({
  library(ggplot2)
})

# root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analyses_dir <- file.path(root_dir, 'analyses', 'dna_repair_index')
output_dir <- file.path(analyses_dir, 'output')

# source function to calculate and plot DNA repair index
source(file.path(analyses_dir, 'util', 'dna_repair_index.R'))

p <- dna_repair_index(expr = file.path(root_dir, 'data', 'tpm_matrix.rds'), 
                      meta = file.path(root_dir, 'data', 'metadata.rds'), 
                      method = "t.test") 
ggsave(filename = file.path(output_dir, 'dna_repair_index_plot.pdf'), device = "pdf", plot = p,width = 10, height = 7)

