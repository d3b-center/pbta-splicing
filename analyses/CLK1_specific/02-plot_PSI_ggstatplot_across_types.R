################################################################################
# PSI_ggstatplot_across_types.R
# script that plots a table of splicing events overlapping uniprot sites
# written by Ammar Naqvi
#
# usage: Rscript splicing_functional_sites.R <file>
################################################################################

suppressPackageStartupMessages(library(ggthemes))
suppressPackageStartupMessages(library(extrafont))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggstatsplot))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
root_dir <- "/Users/naqvia/Desktop/pbta-splicing_git/pbta-splicing"
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1_specific")

input_dir <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")

## get data files 
file_pos   = file.path(results_dir,"rMATS_output.anno.pos.txt")
file_neg   = file.path(results_dir,"rMATS_output.anno.neg.txt")
file_combo = file.path(results_dir,"rMATS_output.anno.txt") 

dpsi_pos <- read.table(file_pos, header=TRUE,sep = "\t")
dpsi_neg <- read.table(file_neg, header=TRUE,sep = "\t")
dpsi_combo <- read.table(file_combo, header=TRUE,sep = "\t")

set.seed(123)
plot1 <- ggstatsplot::ggbetweenstats(
  data = dpsi_pos, 
  x = Type, 
  y = IncLevelDifference,
  k = 3,
  nboot = 15,
  outlier.label = geneSymbol, # label to attach to outlier values
  outlier.label.args = list(color = "red"), # outlier point label color
  notch = TRUE,
  mean.ci = TRUE,
  outlier.tagging = TRUE,
  type = "robust",
  xlab = "Unipro-defined Site",
  ylab = "dPSI",
  pairwise.comparisons = FALSE,
  messages = FALSE
) + theme_Publication()

plot_file = file.path(plots_dir,"dPSI_across_gene-types_neg.pdf") 
ggsave(plot_file, width = 15, height = 5)

set.seed(123)
plot2 <- ggstatsplot::ggbetweenstats(
  data = dpsi_neg, 
  x = Type, 
  y = IncLevelDifference,
  k = 3,
  nboot = 15,
  outlier.label = geneSymbol, # label to attach to outlier values
  outlier.label.args = list(color = "red"), # outlier point label color
  notch = TRUE,
  mean.ci = TRUE,
  outlier.tagging = TRUE,
  type = "robust",
  xlab = "Type",
  ylab = "PSI",
  pairwise.comparisons = FALSE,
  messages = FALSE
) + theme_Publication() 

plot_file = file.path(plots_dir,"dPSI_across_gene-types_pos.pdf") 
ggsave(plot_file, width = 15, height = 5)

dodge <- position_dodge(width = 1)
plot_comb <- ggplot(dpsi_combo, aes(x = Type, y = IncLevelDifference, fill=Type)) + 
  geom_violin(trim=FALSE,  position = dodge, color="black") +
  theme_Publication() + 
  geom_boxplot(width=.1, outlier.colour="NA", position="dodge")