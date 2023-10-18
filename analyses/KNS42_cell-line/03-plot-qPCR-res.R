################################################################################
# 03-plot-qPCR-res.R
# script that generates line graph and compares Treatment absorbance levels
# written by Ammar Naqvi
#
# usage: Rscript 03-plot-qPCR-res.R
################################################################################

## libraries used 
suppressPackageStartupMessages({
  library("vroom")
  library("tidyverse")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "KNS42_cell-line")
input_dir <- file.path(analysis_dir, "input")
plots_dir <- file.path(analysis_dir, "plots")

## theme for all plots
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir,"theme_for_plots.R"))

## output for plot
file_plot = file.path(plots_dir,"qPCR-morp.pdf")

qpcr_res_file = file.path(input_dir,"qpcr-res.tsv")
qpcr_fc_df <- vroom(qpcr_res_file, delim="\t") %>%  pivot_longer(cols=c('Control','1uM','5uM','10uM'), 
                                                                 names_to = "Concentration",
                                                                 values_to = "Fold change") 

# specify order
qpcr_fc_df$Concentration <- factor(qpcr_fc_df$Concentration, levels = c('Control','1uM','5uM','10uM'))


# Save plot as PDF
pdf(file_plot, width =8, height = 5)

# generat grouped barpot
ggplot(qpcr_fc_df, aes(y=`Fold change`, x=Concentration, fill=`Splice junction`)) + 
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values=c("#FFC20A","#0C7BDC")) + 
  theme_Publication()

dev.off()