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
file_dpsi_kinase_plot <- file.path(analysis_dir, "plots", "dPSI_across_functional_sites_kinase.pdf")

## retrieve psi values from tables
file_psi_pos_func <-  file.path(results_dir,"splicing_events.total.pos.intersectUnip.ggplot.txt")
file_psi_neg_func <-  file.path(results_dir,"splicing_events.total.neg.intersectUnip.ggplot.txt")

## read table of recurrent functional splicing (skipping)
dpsi_unip_pos <- read.table(file_psi_pos_func, header=TRUE,sep = "\t") %>% mutate(gene=str_match(SpliceID, "(\\w+[\\.\\d]*)\\_")[, 2]) 

## read table of recurrent functional splicing (inclusion) 
dpsi_unip_neg <- read.table(file_psi_neg_func, header=TRUE,sep = "\t") %>% mutate(gene=str_match(SpliceID, "(\\w+[\\.\\d]*)\\_")[, 2]) 

## ggstatplot across functional sites
set.seed(123)
plot_incl <- ggstatsplot::ggbetweenstats(
  data = dpsi_unip_neg, 
  x = Uniprot, 
  y = dPSI,
  k = 3,
  nboot = 15,
  outlier.label = gene, # label to attach to outlier values
  outlier.label.args = list(color = "red", size=1.8), # outlier point label color
  notch = TRUE,
  mean.ci = TRUE,
  outlier.tagging = TRUE,
  type = "robust",
  xlab = "Unipro-defined Site",
  pairwise.comparisons = FALSE,
  options(ggrepel.max.overlaps = 10),
  messages = FALSE
) + theme_Publication() + labs(y=expression(Delta*PSI)) 

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
  outlier.label = gene, # label to attach to outlier values
  outlier.label.args = list(color = "red", size = 1.8), # outlier point label color
  notch = TRUE,
  mean.ci = TRUE,
  outlier.tagging = TRUE,
  type = "robust",
  xlab = "Unipro-defined Site",
  pairwise.comparisons = FALSE,
  options(ggrepel.max.overlaps = 10),
  messages = FALSE) + 
  theme_Publication() + labs(y=expression(Delta*PSI)) 

# Save plot as PDF
pdf(file_dpsi_skip_plot, 
    width = 15, height = 5)
plot_skip
dev.off()

# kinase gene list
known_kinase_file <- file.path(input_dir,'kinase_known.txt')
known_kinase_df <- read.table(known_kinase_file,header=FALSE) %>% 
  rename('gene' = V1)

dpsi_unip_pos_kinase <- dplyr::inner_join(dpsi_unip_pos, known_kinase_df, by='gene') %>% mutate(Preference="Skipping")
dpsi_unip_neg_kinase <- dplyr::inner_join(dpsi_unip_neg, known_kinase_df, by='gene') %>% mutate(Preference="Inclusion")
dpsi_unip_both_kinase <- rbind(dpsi_unip_pos_kinase,dpsi_unip_neg_kinase) 

## make sina plot
set.seed(45)
kinase_dpsi_plot <- ggplot(dpsi_unip_both_kinase,aes(Preference,dPSI)) +  
  ylab(expression(bold("dPSI"))) +
  geom_violin() +
  ggforce::geom_sina(aes(color = Preference), size = 2,method="density") +
  geom_label_repel(box.padding = 0.5, min.segment.length = 0.5,max.overlaps =Inf, aes(label = gene), data=dpsi_unip_both_kinase %>% subset(gene=='CLK1'), size=2) +
  scale_color_manual(name = "Preference", values = c(Skipping = "#0C7BDC", Inclusion = "#FFC20A")) + 
  theme_Publication() + labs(y=expression(Delta*PSI)) + 
  stat_compare_means(position = "identity", label.x = 1) + 
  theme(legend.position="none")

pdf(file_dpsi_kinase_plot, 
    width = 5, height = 5)
kinase_dpsi_plot
dev.off()