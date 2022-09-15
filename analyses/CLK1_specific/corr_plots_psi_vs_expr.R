################################################################################
# corr_plots_psi_vs_expr.R
# written by Ammar Naqvi
#
# usage: Rscript corr_plots_psi_vs_expr.R <file>
################################################################################

suppressPackageStartupMessages({
  library("ggplot2")
  library("dplyr")
  library("tidyverse")
  
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1_specific")

input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

## output files for final plots
file_corr_plot <- file.path(analysis_dir, "plots", "corr_CLK1_spl_vs_expr.pdf")

input_file = file.path(results_dir, "corr_diffSplice_diffExpr.corr.v2.txt")
tab=read.delim(file,header=TRUE,sep = "\t", row.names = NULL)

## make corr plot
p <- ggplot(tab, aes(x = row.names,y=Gene.1)) + 
  geom_boxplot(fill="lightgreen") + coord_flip() + geom_dotplot(binaxis='y', stackdir='centerwhole', dotsize=.3) +
  theme_Publication() + 
  ggtitle("Splicing vs Expression") +
  xlab("") + ylab("Corr") 
p 

## save plot
ggsave(file_corr_plot, width = 15, height = 5)
