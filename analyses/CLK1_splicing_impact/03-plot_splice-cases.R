################################################################################
# 04-plot_splice-cases.R
# written byAmmar Naqvi and Jo Lynne Rokita
#
# usage: Rscript 04-plot_splice-cases.R
################################################################################

suppressPackageStartupMessages({
  library("ggplot2")
  library("dplyr")
  library("tidyverse")
  library("ggrepel")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories and file paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1_splicing_impact")

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
file_dpsi_es_plot <- file.path(plots_dir,"dPSI_distr_es.pdf")
file_dpsi_ei_plot <- file.path(plots_dir,"dPSI_distr_ei.pdf")


## get and setup input
rmats_merged_file  <- file.path(input_dir,"morpholno.merged.rmats.tsv")
splicing_df  <-  vroom(rmats_merged_file, comment = "#", delim="\t") %>% 
                 filter( (IncLevelDifference >= .20) | (IncLevelDifference <= -.20))  %>%
                 mutate(Type = case_when(IncLevelDifference >= .20 ~ "CLK1-Ex4 Mediated Skipping", 
                                         IncLevelDifference <= -.20 ~ "CLK1-Ex4 Mediated Inclusion"))

splice_case_counts_df <- splicing_df %>% dplyr::count(splicing_case, Type) 

# Specify colors
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#AA4499","#D55E00")


case_colors = list(
  splicing_case = c(A3SS=safe_colorblind_palette[1], A5SS=safe_colorblind_palette[3],RI=safe_colorblind_palette[4],SE=safe_colorblind_palette[5], MXE=safe_colorblind_palette[6]) ) 

piechart_plot<- ggplot(data = splice_case_counts_df, aes(x = "", y = n, fill = splicing_case )) + 
  geom_bar(stat = "identity", position = position_fill()) +
  geom_text(aes(label = n), position = position_fill(vjust = 0.5)) +
  coord_polar(theta = "y",start=0) +
  facet_wrap(~ Type, ncol=1)  +
  scale_fill_manual(name = "Splicing Case",values = case_colors[['splicing_case']]) +
  xlab("") + ylab("") +
  theme_Publication() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=10))

# save plot tiff version
tiff(piechart_plot_of_splicing_case, height =2000, width = 1500, res = 300)
print(piechart_plot)
dev.off()


splicing_df_ES <- splicing_df %>% filter(IncLevelDifference >= .20) 
splicing_df_EI <- splicing_df %>% filter( (IncLevelDifference <= -.20) 
                                                        
## ggstatplot across functional sites
set.seed(123)
plot_es <- ggstatsplot::ggbetweenstats(
  data = splicing_df_ES, 
  x = splicing_case, 
  y = IncLevelDifference,
  k = 3,
  nboot = 15,
  outlier.label = geneSymbol, # label to attach to outlier values
  outlier.label.args = list(color = "red", size=1.8), # outlier point label color
  notch = TRUE,
  mean.ci = TRUE,
  outlier.tagging = TRUE,
  type = "robust",
  pairwise.comparisons = FALSE,
  messages = FALSE
) + theme_Publication() + labs(y=expression(Delta*PSI)) 

set.seed(123)
plot_ei <- ggstatsplot::ggbetweenstats(
  data = splicing_df_EI, 
  x = splicing_case, 
  y = IncLevelDifference,
  k = 3,
  nboot = 15,
  outlier.label = geneSymbol, # label to attach to outlier values
  outlier.label.args = list(color = "red", size=1.8), # outlier point label color
  notch = TRUE,
  mean.ci = TRUE,
  outlier.tagging = TRUE,
  type = "robust",
  pairwise.comparisons = FALSE,
  messages = FALSE
) + theme_Publication() + labs(y=expression(Delta*PSI))

# Save plot as PDF
pdf(file_dpsi_es_plot, 
    width = 15, height = 5)
plot_es
dev.off()

# Save plot as PDF
pdf(file_dpsi_ei_plot, 
    width = 15, height = 5)
plot_ei
dev.off()
