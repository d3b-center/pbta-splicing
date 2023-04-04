################################################################################
# 03-SRSF11_plots.R
# written by Ammar Naqvi
#
# This script uses SRSF11 specific input files to compute and plot correlations 
# between PSI vs expression and also makes barplot of splicing changes across 
# samples. 
#
# usage: Rscript 03-SRSF11_plots.R
################################################################################

suppressPackageStartupMessages({
  library("reshape2")
  library("tidyverse")
  library("ggpubr")
  library("ggplot2")
  library("vroom")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "splicing-factor_dysregulation")

input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

##theme for all plots
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## output files for final plots
file_SRSF11_plot <- file.path(analysis_dir, "plots", "SRSF11_hgg_stacked.tiff")
file_SRSF11_corr_plot <- file.path(analysis_dir, "plots", "corr_rna_vs_psi_SRSF11.tiff")

## get SRSF11 psi values in tumors and ctrls
psi_tumors_file <- "/se.srsf11.incl.txt"
psi_ctrl_file <- "/se.srsf11.ctrl.incl.txt"

psi_tab  <-  read.delim(paste0(input_dir, psi_tumors_file), sep = "\t", header=FALSE) 
colnames(psi_tab) = c("Kids_First_Biospecimen_ID","psi","type")

psi_ctrl_tab  <-  read.delim(paste0(input_dir, psi_ctrl_file), sep = "\t", header=FALSE) 
colnames(psi_ctrl_tab) = c("Kids_First_Biospecimen_ID","psi","type")  
psi_ctrl_tab <- psi_ctrl_tab %>% filter(Kids_First_Biospecimen_ID=="control-BS")

## intersect just midline HGGs
clin_file = file.path(data_dir,"histologies.tsv")
clin_df  <-  vroom(clin_file, comment = "#",delim="\t")

psi_with_midline_filter <- clin_df %>% 
  # Only include RNA-Seq samples 
  filter(experimental_strategy == "RNA-Seq") %>% 
  # Only include samples from PBTA cohort
  filter(cohort == "PBTA") %>%
  # Include samples that are "stranded"
  filter(RNA_library == "stranded") %>% 
  # Include those samples originated from the Midline
  filter(CNS_region == 'Midline') %>% 
  # Include only HGAT sample histology
  filter(short_histology == 'HGAT') %>% inner_join(psi_tab, by="Kids_First_Biospecimen_ID") %>% 
  dplyr::select("Kids_First_Biospecimen_ID","psi","type") %>% distinct()

psi_with_midline_filter_and_ctrl <- rbind(psi_ctrl_tab,psi_with_midline_filter)

## stacked barplot 
stacked_barplot_SRSF11 <- psi_with_midline_filter_and_ctrl %>% 
  arrange(desc(psi)) %>%
  mutate(sample=fct_reorder(Kids_First_Biospecimen_ID,psi)) %>% 
  ggplot(aes(x = Kids_First_Biospecimen_ID,y = psi, fill= type )) +
  geom_bar(position="stack", stat="identity")    + 
  scale_fill_manual(values=c("#FFC20A","#0C7BDC"))      + 
  theme_Publication() +  ylab("Isoform %") + xlab("Sample") + 
  theme(axis.text.x=element_blank()) 


# Save plot as tiff
tiff(file_SRSF11_plot, 
    res = 600, width = 16, height = 8, units = "in")
stacked_barplot_SRSF11
dev.off()

