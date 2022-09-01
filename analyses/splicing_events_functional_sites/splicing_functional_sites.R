################################################################################
# splicing_functional_sites.R
# script that plots a table of splicing events overlapping uniprot sites
# written by Ammar Naqvi
#
# usage: Rscript splicing_functional_sites.R <file>
################################################################################

suppressPackageStartupMessages({
  library("ggstatsplot")
  library("dplyr")
  library("ggplot2")
  library("dplyr")
  library("tidyverse")

})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "splicing_events_functional_sites")

input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

## output files for final plots
file_dpsi_skip_plot <- file.path(analysis_dir, "plots", "dPSI_across_functional_sites_pos.pdf")
file_dpsi_incl_plot <- file.path(analysis_dir, "plots", "dPSI_across_functional_sites_neg.pdf")

## retrieve psi values from tables
dpsi_unip_pos <- read.table(file_pos, header=TRUE,sep = "\t") ## read table of recurrent functional splicing (skipping)
dpsi_unip_neg <- read.table(file_neg, header=TRUE,sep = "\t") ## read table of recurrent functional splicing (inclusion)

## ggstatplot across functional sites
set.seed(123)
plot1 <- ggstatsplot::ggbetweenstats(
  data = dpsi_unip_neg, 
  x = Uniprot, 
  y = dPSI,
  k = 3,
  nboot = 15,
  outlier.label = SpliceID, # label to attach to outlier values
  outlier.label.args = list(color = "red"), # outlier point label color
  notch = TRUE,
  mean.ci = TRUE,
  outlier.tagging = TRUE,
  #title = "Tumor supressors",
  type = "robust",
  xlab = "Unipro-defined Site",
  pairwise.comparisons = TRUE,
  messages = FALSE
) + theme_Publication()

set.seed(123)
plot2 <- ggstatsplot::ggbetweenstats(
  data = dpsi_unip_pos, 
  x = Uniprot, 
  y = dPSI,
  k = 3,
  nboot = 15,
  outlier.label = SpliceID, # label to attach to outlier values
  outlier.label.args = list(color = "red"), # outlier point label color
  notch = TRUE,
  mean.ci = TRUE,
  #outlier.tagging = TRUE,
  type = "robust",
  xlab = "Unipro-defined Site",
  pairwise.comparisons = TRUE,
  messages = FALSE
) + theme_Publication()

## save plots 
ggsave(file_dpsi_skip_plot, width = 15, height = 5)
ggsave(file_dpsi_incl_plot, width = 15, height = 5)
