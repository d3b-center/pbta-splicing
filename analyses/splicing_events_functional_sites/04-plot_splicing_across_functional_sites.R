################################################################################
# 04-plot_splicing_across_functional_sites.R
# script that plots a table of splicing events overlapping uniprot sites
# written by Ammar Naqvi
#
# usage: Rscript 04-plot_splicing_across_functional_sites.R 
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

plots_dir <- file.path(analysis_dir, "plots")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

##theme for all plots
# source function for theme for plots survival
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## output files for final plots
file_dpsi_skip_plot <- file.path(analysis_dir, "plots", "dPSI_across_functional_sites_pos.pdf")
file_dpsi_incl_plot <- file.path(analysis_dir, "plots", "dPSI_across_functional_sites_neg.pdf")

## retrieve psi values from tables
file_psi_pos_func <-  file.path(results_dir,"splicing_events.total.pos.intersectUnip.ggplot.txt")
file_psi_neg_func <-  file.path(results_dir,"splicing_events.total.neg.intersectUnip.ggplot.txt")

dpsi_unip_pos <- read.table(file_psi_pos_func, header=TRUE,sep = "\t") ## read table of recurrent functional splicing (skipping)
dpsi_unip_neg <- read.table(file_psi_neg_func, header=TRUE,sep = "\t") ## read table of recurrent functional splicing (inclusion)

## ggstatplot across functional sites
set.seed(123)
plot_incl <- ggstatsplot::ggbetweenstats(
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
  type = "robust",
  xlab = "Unipro-defined Site",
  pairwise.comparisons = FALSE,
  messages = FALSE
) + theme_Publication()

# Save plot as PDF
pdf(file_dpsi_incl_plot, 
    width = 15, height = 5)
plot_incl
dev.off()


set.seed(123)
plot_skip <- ggstatsplot::ggbetweenstats(
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
  pairwise.comparisons = FALSE,
  messages = FALSE
) + theme_Publication()

# Save plot as PDF
pdf(file_dpsi_skip_plot, 
    width = 15, height = 5)
plot_skip
dev.off()

