# Author: Komal S. Rathi
# GSVA using proliferative genes
# cell cycle genes obtained from table S2 of PMID: 30041684

# load libraries
suppressPackageStartupMessages({
  library(ggplot2)
})

# root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analyses_dir <- file.path(root_dir, 'analyses', 'proliferative_index')
output_dir <- file.path(analyses_dir, 'output')

# source function to calculate and plot proliferative index
source(file.path(analyses_dir, 'util', 'proliferative_index.R'))

# read cell cycle genes
cell_cycle_sig <- file.path(analyses_dir, 'input', 'cell-cycle-genes.txt')
cell_cycle_sig <- read.delim(cell_cycle_sig, header = F, stringsAsFactors = F)
cell_cycle_sig <- list(cell_cycle = cell_cycle_sig$V1)

p <- proliferative_index(expr = file.path(root_dir, 'data', 'tpm_matrix.rds'), 
                         meta = file.path(root_dir, 'data', 'metadata.rds'), 
                         sig = cell_cycle_sig, 
                         method = "t.test") 
ggsave(filename = file.path(output_dir, 'proliferative_index_plot.pdf'), device = "pdf", plot = p, width = 10, height = 7)

