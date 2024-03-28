################################################################################
# 00-subset-CLK1-from-rMATs.R
# written by Ammar Naqvi
#
# usage: Rscript 00-subset-CLK1-from-rMATs.R
################################################################################

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(vroom)
  library(data.table)
})


## Set directories
# Input directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing_correlations")
input_dir <- file.path(analysis_dir, "input")

# Specify file paths
rmats_file <- file.path(data_dir, "splice-events-rmats.tsv.gz")

## output file
rmats_file <- file.path(data_dir, "splice-events-rmats.tsv.gz")

## Load rmats file
rmats_df <-  fread(rmats_file) %>%
  # Select CLK1 gene
  filter(geneSymbol=="CLK1") 

write_tsv(rmats_df,file.path(input_dir, "CLK1-rmats.tsv"))