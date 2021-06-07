# Author: Komal S. Rathi
# Function: Brain cell type proportion analysis using BRETIGEA 

# load libraries
library(BRETIGEA, quietly = TRUE)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(reshape2)

# root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analyses_dir <- file.path(root_dir, 'analyses', 'bretigea')
output_dir <- file.path(analyses_dir, 'output')

# source plotting theme
source(file.path(root_dir, 'util', 'pubTheme.R'))

# function to compute and plot Brain Cell Type Specific Gene Expression scores
plot_bretigea <- function(expr, meta, title, method){
  
  expr <- readRDS(expr)
  meta <- readRDS(meta)
  ct_res <- brainCells(expr, nMarker = 50)
  ct_res <- melt(as.matrix(ct_res), varnames = c("sample", "cell_type"), value.name = "bretigea_score")
  ct_res <- ct_res %>%
    inner_join(meta, by = 'sample')
  
  # shapiro test of normality
  cell_types <- unique(ct_res$cell_type)
  for(i in 1:length(cell_types)){
    print(cell_types[i])
    print(shapiro.test(x = ct_res[ct_res$type == "pos" & ct_res$cell_type == cell_types[i], 'bretigea_score'])) 
    print(shapiro.test(x = ct_res[ct_res$type == "neg" & ct_res$cell_type == cell_types[i], 'bretigea_score'])) 
  }
  
  # change factor for variable
  levels(ct_res$cell_type) <- c("Astrocytes", "Endothelial Cells", "Microglia",
                                "Neurons", "Oligodendrocytes", "Oligodendrocyte PC") 
  
  my_comparisons <- list(c("pos", "neg"))
  p <- ggplot(ct_res, aes(x = type, y = bretigea_score, fill = type)) + 
    facet_wrap(~cell_type) +
    stat_boxplot(geom ='errorbar', width = 0.2) +
    geom_boxplot(lwd = 0.5, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) +
    scale_x_discrete(labels = c("pos" = "clk1_splice_pos", 
                                "neg" = "clk1_splice_neg")) +
    scale_fill_manual(values = c("pos" = "#F8766D",
                                 "neg" = "#00BFC4")) +
    theme_Publication() + xlab("") + ylab("BRETIGEA score") +
    guides(fill = FALSE) + ggtitle(title) +
    stat_compare_means(label.y = 0.9, color = "darkred", paired = FALSE, method = method)
  return(p)
}

p <- plot_bretigea(expr = file.path(root_dir, 'data', 'pbta-hgat-dx-prog-pm-gene-expression-subset-rsem-tpm-uncorrected.rds'), 
                   meta = file.path(root_dir, 'data', 'pbta-hgat-dx-prog-pm-gene-expression-subset-metadata.rds'), 
                   title = "Brain cell type proportions (nSamples = 8)", method = "t.test")
ggsave(filename = file.path(output_dir, 'bretigea_plot.pdf'), plot = p, device = "pdf", width = 10, height = 7)


