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
#file_dpsi_skip_plot <- file.path(analysis_dir, "plots", "dPSI_across_functional_sites_pos.pdf")
#file_dpsi_incl_plot <- file.path(analysis_dir, "plots", "dPSI_across_functional_sites_neg.pdf")
file_dpsi_plot <- file.path(analysis_dir, "plots", "dPSI_across_functional_sites.pdf")


## retrieve psi values from tables
file_psi_pos_func <-  file.path(results_dir,"splicing_events.total.pos.intersectUnip.ggplot.txt")
file_psi_neg_func <-  file.path(results_dir,"splicing_events.total.neg.intersectUnip.ggplot.txt")

## get gene_lsits
gene_list_df = vroom(file.path(input_dir,"gene_lists.tsv")) %>% 
  ## reformat and rename to plot
  pivot_longer(cols=c('Kinase','Epigenetic','SWISNF','Cancer','SF'), 
               names_to = "Annotation",
               values_to = "gene") %>% 
  filter(Annotation != 'SWISNF')

## read table of recurrent functional splicing (skipping)
dpsi_unip_pos <- read.table(file_psi_pos_func, header=TRUE,sep = "\t") %>% 
  mutate(gene=str_match(SpliceID, "(\\w+[\\.\\d]*)\\_")[, 2]) %>%
  mutate('Type'="skipping")

## read table of recurrent functional splicing (inclusion) 
dpsi_unip_neg <- read.table(file_psi_neg_func, header=TRUE,sep = "\t") %>% 
  mutate(gene=str_match(SpliceID, "(\\w+[\\.\\d]*)\\_")[, 2]) %>% 
  mutate('Type'="inclusion")

dpsi_unip_total <- rbind(dpsi_unip_pos,dpsi_unip_neg) %>% 
  full_join(gene_list_df,by='gene',relationship = "many-to-many") %>% 
  replace(is.na(.), "Other") %>% 
  mutate(dPSI=abs(as.numeric(dPSI))) %>% 
  mutate(outlier = ifelse(find_outlier(dPSI), dPSI, NA))


dpsi_unip_total$Annotation <- factor(dpsi_unip_total$Annotation, levels = c('Kinase','Epigenetic','SF','Cancer','Other')) 

## make plot
set.seed(45)
sbi_vs_incl_plot <- ggplot(dpsi_unip_total,aes(Uniprot,dPSI)) +  
  xlab(expression(bold("Uniprot Annotation"))) + ylab(expression(bold("dPSI"))) +
  facet_grid(rows = vars(Annotation)) + 
  ggforce::geom_sina(aes(color = Type), size = 1,method="density") +
  scale_color_manual(name = "Type", values = c(skipping = "#FFC20A", inclusion = "#0C7BDC")) + 
  theme_Publication()

# Save plot as PDF
pdf(file_dpsi_plot, 
    width = 15, height = 5)
sbi_vs_incl_plot
dev.off()