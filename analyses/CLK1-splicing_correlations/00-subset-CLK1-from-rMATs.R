################################################################################
# 00-subset-CLK1-from-rMATs.R
# written by Ammar Naqvi
#
# usage: Rscript 00-subset-CLK1-from-rMATs.R
################################################################################

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})


## Set directories
# Input directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing_correlations")
results_dir <- file.path(analysis_dir, "results")

# Specify file paths
rmats_file <- file.path(data_dir, "splice-events-rmats.tsv.gz")

## Load rmats file
rmats_df <-  fread(rmats_file) 

rmats_clk1 <- rmats_df %>%
  # Select CLK1 gene
  filter(geneSymbol=="CLK1") %>%
  write_tsv(file.path(results_dir, "clk1-splice-events-rmats.tsv"))

rmats_nf1 <- rmats_df %>%
  # Select NF1 gene
  filter(geneSymbol=="NF1") %>%
  write_tsv(file.path(results_dir, "nf1-splice-events-rmats.tsv"))
