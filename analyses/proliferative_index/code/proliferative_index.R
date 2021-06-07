# Author: Komal S. Rathi
# GSVA using proliferative genes
# cell cycle genes obtained from table S2 of PMID: 30041684

# load libraries
library(GSVA)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(reshape2)

# root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analyses_dir <- file.path(root_dir, 'analyses', 'proliferative_index')
output_dir <- file.path(analyses_dir, 'output')

# source plotting theme
source(file.path(root_dir, 'util', 'pubTheme.R'))

# read cell cycle genes
cell_cycle_sig <- file.path(analyses_dir, 'input', 'cell-cycle-genes.txt')
cell_cycle_sig <- read.delim(cell_cycle_sig, header = F, stringsAsFactors = F)
cell_cycle_sig <- list(cell_cycle = cell_cycle_sig$V1)

proliferative_index <- function(expr, meta, sig = cell_cycle_sig, title, method){
  expr <- readRDS(expr)
  meta <- readRDS(meta)
  
  # log2 FPKM data
  expr <- log2(expr+1)
  
  # GSVA
  ssgsea_output <- gsva(as.matrix(expr), 
                        gset.idx.list = sig, 
                        method = "gsva",
                        min.sz = 1, max.sz = 1500,
                        mx.diff = TRUE)
  ssgsea_output <- melt(ssgsea_output, varnames = c("molecular_function", "sample"), value.name = "proliferative_index")
  ssgsea_output <- ssgsea_output %>%
    inner_join(meta, by = 'sample')
  
  # shapiro test of normality
  print(shapiro.test(x = ssgsea_output[ssgsea_output$type == "pos", 'proliferative_index'])) 
  print(shapiro.test(x = ssgsea_output[ssgsea_output$type == "neg", 'proliferative_index'])) 
  
  # plot
  my_comparisons <- list(c("pos", "neg"))
  p <- ggplot(ssgsea_output, aes(x = type, y = proliferative_index, fill = type)) + 
    facet_wrap(~molecular_function) +
    stat_boxplot(geom ='errorbar', width = 0.2) +
    geom_boxplot(lwd = 0.5, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) +
    scale_x_discrete(labels = c("pos" = "clk1_splice_pos", 
                                "neg" = "clk1_splice_neg")) +
    scale_fill_manual(values = c("pos" = "#F8766D",
                                 "neg" = "#00BFC4")) +
    theme_Publication() + xlab("") + ylab("Proliferative Index") +
    guides(fill = FALSE) + ggtitle(title) +
    stat_compare_means(label.y = 0.9, color = "darkred", paired = FALSE, method = method)
  return(p)
}

p <- proliferative_index(expr = file.path(root_dir, 'data', 'pbta-hgat-dx-prog-pm-gene-expression-subset-rsem-tpm-uncorrected.rds'), 
                         meta = file.path(root_dir, 'data', 'pbta-hgat-dx-prog-pm-gene-expression-subset-metadata.rds'), 
                         sig = cell_cycle_sig,
                         title = "Proliferative Index (nSamples = 8)", method = "t.test") 
ggsave(filename = file.path(output_dir, 'proliferative_index_plot.pdf'), device = "pdf", plot = p, width = 10, height = 7)

