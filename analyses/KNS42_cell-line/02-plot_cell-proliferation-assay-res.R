################################################################################
# 02-plot_cell-proliferation-assay-res.R
# script that generates line graph and compares Treatment absorbance levels
# written by Ammar Naqvi
#
# usage: Rscript 02-plot_cell-proliferation-assay-res.R 
################################################################################

## libraries used 
suppressPackageStartupMessages({
library("vroom")
library("tidyverse")
library("ggpubr")
library("rstatix")
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
figures_dir = "~/d3b_coding/pbta-splicing/figures/theme_for_plots.R"
source(figures_dir)

## output for plot
file_line_plot = file.path(plots_dir,"cell_prolif-line.pdf")

cell_prolif_res_file = file.path(input_dir,"cell_prolif_res.tsv")
cell_prolif_df <- vroom(file, delim=) %>% 
  pivot_longer(
    c(-Time), 
    names_pattern = "(.+)_([0-9])", 
    names_to = c("Treatment", "Rep"),
    values_to = "Absorbance"
  ) 

## stat test with ref
stat.test_ref <- cell_prolif_df %>%
  group_by(Time) %>%
  t_test(Absorbance ~ Treatment, ref.group = "CLK1") 


# create a line plot with error bars (mean +/- sd)
lp <- ggline(
  cell_prolif_df, x = "Time", y = "Absorbance", add = "mean_sd", 
  color = "Treatment", palette = c("red", "blue", "darkgreen")
)

# add p-values onto the line plots
stat.test_ref <- stat.test_ref %>%
  add_xy_position(fun = "mean_sd", x = "Time") 

lp + stat_pvalue_manual(
  stat.test_ref,  label = "p.adj.signif", 
  tip.length = 0, linetype  = "blank"
)

# move down the significance levels using vjust
lp + stat_pvalue_manual(
  stat.test,  label = "p.adj.signif", 
  linetype  = "blank", hjust = 0, vjust = 9
) + theme_Publication()



## theme for all plots
# source function for theme for plots survival
figures_dir <- file.path(root_dir, "figures")

figures_dir = "~/d3b_coding/pbta-splicing/figures/theme_for_plots.R"
source(figures_dir)


file = "/Users/naqvia/Documents/cell_prolif_res.tsv"
cell_prolif_df <- vroom(file, delim=) %>% 
  pivot_longer(
    c(-Time), 
    names_pattern = "(.+)_([0-9])", 
    names_to = c("Treatment", "Rep"),
    values_to = "Absorbance"
  ) 

## stat test with ref
stat.test_ref <- cell_prolif_df %>%
  group_by(Time) %>%
  t_test(Absorbance ~ Treatment, ref.group = "CLK1") 

# add p-values onto the line plots
stat.test_ref <- stat.test_ref %>%
  add_xy_position(fun = "mean_sd", x = "Time") 

# create a line plot with error bars (mean +/- sd)
lp <- ggline(
  cell_prolif_df, x = "Time", y = "Absorbance", add = "mean_sd", 
  color = "Treatment", palette = c("red", "blue", "darkgreen")
)



lp + stat_pvalue_manual(
  stat.test_ref,  label = "p.adj.signif", 
  tip.length = 0, linetype  = "blank"
)

# Save plot as PDF
pdf(file_line_plot, width = 10, height = 5)

# move down the significance levels using vjust
lp + stat_pvalue_manual(
  stat.test,  label = "p.adj.signif", 
  linetype  = "blank", hjust = 0, vjust = 7
) + theme_Publication()

dev.off()

