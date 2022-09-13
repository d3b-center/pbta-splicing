################################################################################
# corr_plots.R
# written by Ammar Naqvi
#
# usage: Rscript SRSF11_plot.R
################################################################################


suppressPackageStartupMessages({
  library(corrplot)
  library(limma)
  library(Hmisc)
  library(corrplot)
  library(plyr)
  library(gridExtra)
  library(grid)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
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

## output files for final plots
file_SRSF11_corr_prot_plot <- file.path(analysis_dir, "plots", "corr_SRSF11_prot.png")
file_RBM5_corr_prot_plot <- file.path(analysis_dir, "plots", "corr_RBM5_prot.png")

file_SRSF11_corr_phos_plot <- file.path(analysis_dir, "plots", "corr_SRSF11_phos.png")
file_RBM5_corr_phos_plot <- file.path(analysis_dir, "plots", "corr_RBM5_phos.png")

## get SRSF11 psi table
file <- "/SRSF11.rna_vs_prot.cptac.txt"
tab_rna_vs_prot  <-  read.delim(paste0(input_dir, file), sep = "\t", header=TRUE)

plot_prot_S11 <-ggscatter(tab_rna_vs_prot, x="RNA", y="Protein", 
                          add = "reg.line", conf.int = TRUE, 
                          cor.coef = TRUE, cor.method = "pearson",
                          add.params = list(color = "red",
                                            fill = "pink"),
                          ticks = TRUE,
                          xlab = "RNA", ylab = "Protein") + theme_Publication()


ggsave(
  file_SRSF11_corr_prot_plot,
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 5000,
  height = 5000,
  units = "px",
  dpi = 800,
  limitsize = TRUE,
  bg = NULL
)


## get SRSF11 phosp table
file <- "/SRSF11rna_vs_phos.txt"
tab_rna_vs_phosp  <-  read.delim(paste0(input_dir, file), sep = "\t", header=TRUE)


plot_phos1 <- ggscatter(tab_rna_vs_phosp, x="RNA", y="S433", 
                        add = "reg.line", conf.int = TRUE, 
                        cor.coef = TRUE, cor.method = "pearson",
                        add.params = list(color = "red",
                                          fill = "pink"),
                        ticks = TRUE,
                        #xticks.by = .1, yticks.by = .1,
                        xlab = "RNA", ylab = "S433") + theme_Publication()

plot_phos2 <-ggscatter(tab_rna_vs_phosp, x="RNA", y="S482", 
                       add = "reg.line", conf.int = TRUE, 
                       cor.coef = TRUE, cor.method = "pearson",
                       add.params = list(color = "red",
                                         fill = "pink"),
                       ticks = TRUE,
                       #xticks.by = .1, yticks.by = .1,
                       xlab = "RNA", ylab = "S482") + theme_Publication()

plot_phos3 <- ggscatter(tab_rna_vs_phosp, x="RNA", y="S448", 
                        add = "reg.line", conf.int = TRUE, 
                        cor.coef = TRUE, cor.method = "pearson",
                        add.params = list(color = "red",
                                          fill = "pink"),
                        ticks = TRUE,
                        #xticks.by = .1, yticks.by = .1,
                        xlab = "RNA", ylab = "S448") + theme_Publication()

grid.arrange(plot_prot_S11, plot_phos1, plot_phos2,plot_phos3, ncol=4, widths=c(2.3,2.3, 2.3, 2.3))
grid.rect(width = .98, height = .98, gp = gpar(lwd = 2, col = "black", fill = NA))


file = "/Users/naqvia/Desktop/pbta-splicing_git/pbta-splicing/analyses/splicing_expr_corr/RBM5rna_vs_phos.txt"
tab=read.table(file,header=TRUE,sep = "\t")

tab_mod <- tab %>%
  na.omit()

scatter_proteo <- ggscatter(tab_mod, x="rna", y="proteo", 
                            add = "reg.line", conf.int = TRUE, 
                            cor.coef = TRUE, cor.method = "pearson",
                            add.params = list(color = "red",
                                              fill = "pink"),
                            ticks = TRUE,
                            #xticks.by = .1, yticks.by = .1,
                            xlab = "RNA", ylab = "Proteo") + theme_Publication()

ggsave(
  file_RBM5_corr_prot_plot,
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 5000,
  height = 5000,
  units = "px",
  dpi = 800,
  limitsize = TRUE,
  bg = NULL
)

scatter_phosp1 <- ggscatter(tab_mod, x="rna", y="S624", 
                            add = "reg.line", conf.int = TRUE, 
                            cor.coef = TRUE, cor.method = "pearson",
                            add.params = list(color = "red",
                                              fill = "pink"),
                            ticks = TRUE,
                            #xticks.by = .1, yticks.by = .1,
                            xlab = "RNA", ylab = "S624") + theme_Publication()

scatter_phosp2 <- ggscatter(tab_mod, x="rna", y="S78", 
                            add = "reg.line", conf.int = TRUE, 
                            cor.coef = TRUE, cor.method = "pearson",
                            add.params = list(color = "red",
                                              fill = "pink"),
                            ticks = TRUE,
                            #xticks.by = .1, yticks.by = .1,
                            xlab = "RNA", ylab = "S78") + theme_Publication()

grid.arrange(scatter_proteo, scatter_phosp1, scatter_phosp2, ncol=3, widths=c(2.3, 2.3, 2.3))
grid.rect(width = .98, height = .98, gp = gpar(lwd = 2, col = "black", fill = NA))
