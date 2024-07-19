# 09-plot_total-splicing-cases.R
# written by Ammar Naqvi
#
# This script plots total splicing cases comparing treated vs ctrl
#
# usage: Rscript 08-plot_total-splicing-cases.R
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
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing-impact-morpholino")
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
## rmats file
rmats_merged_file  <- file.path(data_dir,"morpholno.merged.rmats.tsv")

## extract strong splicing changes 
splicing_df  <-  vroom(rmats_merged_file, comment = "#", delim="\t") %>% 
  filter(FDR < 0.05 & PValue < 0.05) 

## extract strong differential splicing cases (dPSI >= |.10|)
splicing_df_ES <- splicing_df %>% filter(IncLevelDifference  >= .10) %>% mutate(Preference="Skipping")
splicing_df_EI <- splicing_df %>% filter(IncLevelDifference <= -.10) %>% mutate(Preference="Inclusion",
                                                                                IncLevelDifference = abs(IncLevelDifference) )

splice_case_total <- rbind(splicing_df_ES,splicing_df_EI)

splice_case_counts_df <- splice_case_total %>% dplyr::count(splicing_case, Preference) %>% arrange(Preference)

lolliplot_plot <- ggplot(splice_case_counts_df, aes(x=splicing_case, y=n)) +
  geom_segment( aes(x=splicing_case, xend=splicing_case, y=0, yend=n))+
  geom_point(color="black", size=4) +
  coord_flip() +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme_Publication() + 
  xlab("Type of splice event") +
  ylab("Number") + 
  facet_wrap(~Preference) +
  theme(
    axis.title=element_text(size=14,face="bold"),
    axis.text.x = element_text(angle = 75, hjust = 1),
  )


# save plot pdf version
pdf(plot_path, height = 3, width = 4.5)
print(lolliplot_plot)
dev.off()

