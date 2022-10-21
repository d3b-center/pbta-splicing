################################################################################
# 02-RSF11_plots.R
# written by Ammar Naqvi
#
# This script uses SRSF11 specific input files to compute and plot correlations 
# between PSI vs expression and also makes barplot of splicing changes across 
# samples. 
#
# usage: Rscript SRSF11_plot.R
################################################################################


suppressPackageStartupMessages({
  library("sva")
  library("EnhancedVolcano")
  library("DESeq2")
  library("reshape2")
  library("tidyverse")
  library("ggpubr")
})


# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "mRNA_diff_expr")

input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

##theme for all plots
# source function for theme for plots survival
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## output files for final plots
file_SRSF11_plot <- file.path(analysis_dir, "plots", "SRSF11_hgat_stacked.png")
file_SRSF11_corr_plot <- file.path(analysis_dir, "plots", "corr_rna_vs_psi_SRSF11.png")

## get SRSF11 psi table
file <- "/SRSF11_psi.txt"
psi_tab  <-  read.delim(paste0(input_dir, file), sep = "\t", header=TRUE)

## transform and melt data for ggplot
melt_psi_tab <- melt(psi_tab, id=c("sample","id","histology"), variable.name =c("Type"))

## stacked barplot 
melt_psi_tab %>% 
  #arrange(desc(value)) %>%
  mutate(sample=fct_reorder(sample,value)) %>% 
  ggplot(aes(x = sample,y = value, fill= Type ))+ 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values=c("red",
                                                                           "blue")) 
# Save plot as PNG
png(file_SRSF11_plot, 
    res = 800, width = 16, height = 8, units = "in")
last_plot()
dev.off()

## scatter plot of psi vs pct expression
file <- "/SRSF11_splicing_vs_expr.txt"
tab  <-  read.delim(paste0(input_dir, file), sep = "\t", header=TRUE)
ggscatter(tab, x="Expr", y="dPSI", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          add.params = list(color = "blue",
                            fill = "lightgray"),
          ticks = TRUE,
          xticks.by = .1, yticks.by = .1,
          xlab = "Expr", ylab = "dPSI") + theme_Publication()


# Save plot as PNG
png(file_SRSF11_corr_plot, 
    res = 800, width = 4, height = 4, units = "in")
last_plot()
dev.off()


