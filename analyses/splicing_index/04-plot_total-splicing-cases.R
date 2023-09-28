# 04-plot_total-splicing-cases.R
# written by Ammar Naqvi
#
# This script plots total splicing cases across samples
#
# usage: Rscript 04-plot_total-splicing-cases.R
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
analysis_dir <- file.path(root_dir, "analyses", "splicing_index")
results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")

if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

# theme for all plots
# source functions
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## define output files
piechart_plot_of_splicing_case <- file.path(plots_dir,"piechart_splice-types.tiff")

## get and setup input
splice_case_SE_file  <- file.path(results_dir,"splice_events.diff.SE.txt")
splice_case_RI_file  <- file.path(results_dir,"splice_events.diff.RI.txt")
splice_case_A5SS_file  <- file.path(results_dir,"splice_events.diff.A5SS.txt")
splice_case_A3SS_file  <- file.path(results_dir,"splice_events.diff.A3SS.txt")

splice_case_SE_df  <-  vroom(splice_case_SE_file ,delim="\t") 
splice_case_RI_df  <-  vroom(splice_case_RI_file ,delim="\t") 
splice_case_A5SS_df  <-  vroom(splice_case_A5SS_file ,delim="\t") 
splice_case_A3SS_df  <-  vroom(splice_case_A3SS_file ,delim="\t") 

splice_case_total <- rbind(splice_case_SE_df,splice_case_RI_df,splice_case_A5SS_df,splice_case_A3SS_df)
splice_case_counts_df <- splice_case_total %>% dplyr::count(Case, Type) 

# Specify colors
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#AA4499")
                                        
case_colors = list(
  Case = c(A3SS=safe_colorblind_palette[1], A5SS=safe_colorblind_palette[3],RI=safe_colorblind_palette[4],SE=safe_colorblind_palette[5]) )


piechart_plot<- ggplot(data = splice_case_counts_df, aes(x = "", y = n, fill = Case )) + 
  geom_bar(stat = "identity", position = position_fill()) +
  geom_text(aes(label = n), position = position_fill(vjust = 0.5)) +
  coord_polar(theta = "y",start=0) +
  facet_wrap(~ Type, ncol=1)  +
  scale_fill_manual(name = "Splicing Case",values = case_colors[['Case']]) +
  xlab("") + ylab("") +
  theme_Publication() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))

# save plot tiff version
tiff(piechart_plot_of_splicing_case, height =2000, width = 1500, res = 300)
print(piechart_plot)
dev.off()
