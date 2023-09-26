################################################################################
# 01-plot_ont_vs_short_CLK1-Ex4_psi.R
# 
# Plotting script that takes in Exon 4 PSI from short vs ONT reads
#
# written by Ammar Naqvi
#
# usage: Rscript 01-plot_ont_vs_short_CLK1-Ex4_psi.R
################################################################################

suppressPackageStartupMessages({
  library("ggplot2")
  library("dplyr")
  library("tidyverse")
  library("viridis")
  library("RColorBrewer")
  library("vroom")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

##directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "ONT_visualization")

results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")
input_dir <- file.path(analysis_dir, "input")

if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

# theme for all plots
# source functions
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## ont vs short tsv
ont_vs_rmats_comparison_file <- "ont_vs_rmats.ggplot.tsv"
ont_rmats_df  <-  vroom(file.path(input_dir,ont_vs_rmats_comparison_file), comment = "#",delim="\t", col_names = TRUE, trim_ws = TRUE, show_col_types = FALSE)

# stacked barplot
tiff(file.path(plots_dir, "isoform_stackedbarplot.tiff"), height = 1800, width = 1600, units = "px", res = 300)
ggplot(ont_rmats_df, aes(fill=isoform, y=PSI, x=type)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("#FFC20A","#0C7BDC")) +
  ggtitle("long vs short RNA-seq") +
  facet_wrap(~cell_lines) +
  xlab("") + theme_Publication()
dev.off()


