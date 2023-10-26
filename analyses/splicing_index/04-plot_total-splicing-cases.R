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
plot_path <- file.path(plots_dir,"splice-types.pdf")

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
splice_case_counts_df <- splice_case_total %>% dplyr::count(Case, Type) %>% arrange(Type)

lolliplot_plot <- ggplot(splice_case_counts_df, aes(x=Case, y=n)) +
  geom_segment( aes(x=Case, xend=Case, y=0, yend=n))+
  geom_point( color="black", size=3) +
  #scale_fill_manual(name = "Splicing Case",values = case_colors[['Case']]) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  theme_Publication() + 
  xlab("Type of splice event") +
  ylab("Number of splice variants") + 
  facet_wrap(Type ~ .) +
  theme(
    axis.title=element_text(size=14,face="bold"),
    axis.text.x = element_text(angle = 75, hjust = 1),
  )
  

# save plot tiff version
pdf(plot_path, height = 4, width = 5.5)
print(lolliplot_plot)
dev.off()
