# function to compute and plot Brain Cell Type Specific Gene Expression scores
suppressPackageStartupMessages({
  library(BRETIGEA, quietly = TRUE)
  library(tidyverse)
  library(ggpubr)
  library(ggplot2)
  library(reshape2)
})

brain_cells_prop <- function(expr, meta, title, method){
  
  expr <- readRDS(expr)
  meta <- readRDS(meta)
  ct_res <- brainCells(expr, nMarker = 50)
  ct_res <- melt(as.matrix(ct_res), varnames = c("Kids_First_Biospecimen_ID", "cell_type"), value.name = "bretigea_score")
  ct_res <- ct_res %>%
    inner_join(meta, by = 'Kids_First_Biospecimen_ID')
  
  # change factor for variable
  levels(ct_res$cell_type) <- c("Astrocytes", "Endothelial Cells", "Microglia",
                                "Neurons", "Oligodendrocytes", "Oligodendrocyte PC") 
  
  # plot
  p <- ggplot(ct_res, aes(x = type, y = bretigea_score, fill = type)) +
    facet_wrap(~cell_type) + 
    stat_boxplot(geom ='errorbar', width = 0.2) +
    geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.color = "white",
                 lwd = 0.5, fatten = 0.7, width = 0.5) + 
    ggtitle('Brain cells proportion') +
    theme_bw() + theme_pubr(base_size = 10) + ylab('Brain cells proportion') + xlab("") + 
    guides(fill = F) + 
    geom_jitter(position=position_jitter(width=.01), shape = 21) +
    scale_x_discrete(labels = c("high" = "CLK1 high exon inclusion", 
                                "low" = "CLK1 low exon inclusion")) +
    scale_fill_manual(values = c("high" = "#F8766D",
                                 "low" = "#00BFC4")) +
    stat_compare_means(label.y = 1.4, color = "darkred", paired = FALSE, method = method) +
    labs(fill = "Type")
  return(p)
}