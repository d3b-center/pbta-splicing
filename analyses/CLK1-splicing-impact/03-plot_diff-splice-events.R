################################################################################
# 03-plot_diff-splice-events.R
# written byAmmar Naqvi and Jo Lynne Rokita
#
# usage: Rscript 03-plot_diff-splice-events.R
################################################################################

suppressPackageStartupMessages({
  library("ggplot2")
  library("tidyverse")
  library("ggrepel")
  library("vroom")
  library("PMCMRplus")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories and file paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing-impact")

results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

plots_dir <- file.path(analysis_dir, "plots")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

# set file path
plot_file = file.path(plots_dir,"dPSI_volcano_morpholino.pdf")

##theme for all plots
# source function for theme for plots survival
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## define output files
piechart_plot_of_splicing_case <- file.path(plots_dir,"piechart_splice-types.tiff")
file_dpsi_plot <- file.path(plots_dir,"dPSI_distr.pdf")

## get and setup input
## rmats file
rmats_merged_file  <- file.path(analysis_dir,"input","morpholno.merged.rmats.tsv")

## extract strong splicing changes 
splicing_df  <-  vroom(rmats_merged_file, comment = "#", delim="\t") %>% 
                 filter(FDR < 0.05 & PValue < 0.05) 

## extract strong differential splicing cases (dPSI >= |.10|)
splicing_df_ES <- splicing_df %>% filter(IncLevelDifference  >= .10) %>% mutate(Preference="Skipping")
splicing_df_EI <- splicing_df %>% filter(IncLevelDifference <= -.10) %>% mutate(Preference="Inclusion",
                                                                                IncLevelDifference =IncLevelDifference*-1 )

psi_comb <- rbind(splicing_df_EI,splicing_df_ES)

## ggstatplot across functional sites
set.seed(123)
plot_dsp <-  ggplot(psi_comb,aes(splicing_case, IncLevelDifference*100) ) +  
  ylab(expression(Delta*"PSI"))+
  ggforce::geom_sina(aes(color = Preference, alpha = 0.4), pch = 16, size = 4, method="density") +
  geom_boxplot(outlier.shape = NA, color = "black", size = 0.5, coef = 0, aes(alpha = 0.4)) +
  facet_wrap("Preference") +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("A3SS", "A5SS"),
                                                                c("A3SS", "MXE"),
                                                                c("A3SS", "RI"),
                                                                c("A3SS", "SE"),
                                                                c("A5SS", "MXE"),
                                                                c("A5SS", "RI"),
                                                                c("A5SS", "SE"),
                                                                c("MXE", "RI"),
                                                                c("MXE", "SE"),
                                                                c("SE", "A5SS"),
                                                                c("SE", "A3SS"),
                                                                c("SE", "RI"))) + 
  scale_color_manual(name = "Preference", values = c(Inclusion = "#FFC20A", Skipping = "#0C7BDC"))  + 
  theme_Publication() + 
  labs(x= "Splicing Case") + 
  theme(legend.position="none") +
  ylim(c(0,170))


# Save plot as PDF
pdf(file_dpsi_plot, 
    width = 15, height = 7)
plot_dsp
dev.off()
