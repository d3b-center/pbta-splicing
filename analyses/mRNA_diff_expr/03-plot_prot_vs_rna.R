################################################################################
# corr_plots.R
# written by Ammar Naqvi
#
# This script computes and plots correlations between RNA expression and 
# proteommic z-scores obtained from CPTACT portal for dysregulated splicing 
# factors RBM15, NOVA2 and MSI1-- that was previously identified in 
# 01-volcano_plot_mRNA.R script.
#
# usage: Rscript SRSF11_plot.R
################################################################################


suppressPackageStartupMessages({
  library(corrplot)
  library(plyr)
  library(gridExtra)
  library(grid)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
})

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "mRNA_diff_expr")

input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

file_RBM15 = "RBM15_prot_rna.txt"
file_MSI1  = "MSI1_prot_rna.txt"
file_NOVA2 = "NOVA2_prot_rna.txt"
file_SAMD4A = "SAMD4A_prot_rna.txt"

## get SF protein vs rna values
cptac_prot_RBM15 <-  read.delim(paste0(input_dir, "/",file_RBM15), sep = "\t", header=TRUE)
cptac_prot_MSI1  <-  read.delim(paste0(input_dir, "/",file_MSI1),  sep = "\t", header=TRUE)
cptac_prot_NOVA2 <-  read.delim(paste0(input_dir, "/",file_NOVA2), sep = "\t", header=TRUE)
cptac_prot_SAMD4A <-  read.delim(paste0(input_dir, "/",file_SAMD4A), sep = "\t", header=TRUE)

## make scatter plots of each SF 
plot_RBM15 <- ggscatter(cptac_prot_RBM15, x="rna", y="proteo", 
                        add = "reg.line", conf.int = TRUE, 
                        cor.coef = TRUE, cor.method = "pearson",
                        add.params = list(color = "red",
                                          fill = "pink"),
                        ticks = TRUE,
                        title = "RBM15",
                        xlab = "RNA", ylab = "Proteo")

plot_MSI1<- ggscatter(cptac_prot_MSI1, x="rna", y="proteo", 
                        add = "reg.line", conf.int = TRUE, 
                        cor.coef = TRUE, cor.method = "pearson",
                        add.params = list(color = "red",
                                          fill = "pink"),
                        ticks = TRUE,
                        title = "MSI1",
                        xlab = "RNA", ylab = "Proteo")

plot_NOVA2 <- ggscatter(cptac_prot_NOVA2, x="rna", y="proteo", 
                        add = "reg.line", conf.int = TRUE, 
                        cor.coef = TRUE, cor.method = "pearson",
                        add.params = list(color = "red",
                                          fill = "pink"),
                        ticks = TRUE,
                        title = "NOVA2",
                        xlab = "RNA", ylab = "Proteo")

plot_SAMD4A <- ggscatter(cptac_prot_SAMD4A, x="rna", y="proteo", 
                        add = "reg.line", conf.int = TRUE, 
                        cor.coef = TRUE, cor.method = "pearson",
                        add.params = list(color = "red",
                                          fill = "pink"),
                        ticks = TRUE,
                        title = "SAMD4A",
                        xlab = "RNA", ylab = "Proteo")


## arrange all plots in one grid
grid.arrange(plot_RBM15,plot_SAMD4A,plot_MSI1, plot_NOVA2,ncol=1, widths=c(.3)) 
grid.rect(width = .98, height = .98, gp = gpar(lwd = 2, col = "black", fill = NA))

# Save plot as PNG
file_corr_prot_plot <- file.path(analysis_dir, "plots", "corr_SF_prot_rna.png")

corr_plots <- arrangeGrob(plot_RBM15,plot_SAMD4A,plot_MSI1, plot_NOVA2, ncol=1) #generates g
ggsave(file=file_corr_prot_plot, corr_plots) #saves corr_phos_plot


