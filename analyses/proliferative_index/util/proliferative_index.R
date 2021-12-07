# function to calculate proliferative index and create boxplot
suppressPackageStartupMessages({
  library(ggpubr)
  library(ggplot2)
  library(GSVA)
  library(reshape2)
})

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
  ssgsea_output <- melt(ssgsea_output, varnames = c("molecular_function", "Kids_First_Biospecimen_ID"), value.name = "proliferative_index")
  ssgsea_output <- ssgsea_output %>%
    inner_join(meta, by = 'Kids_First_Biospecimen_ID')
  
  # plot
  p <- ggplot(ssgsea_output, aes(x = type, y = proliferative_index, fill = type)) +
    stat_boxplot(geom ='errorbar', width = 0.2) +
    geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.color = "white",
                 lwd = 0.5, fatten = 0.7, width = 0.5) + 
    ggtitle('Proliferative index') +
    theme_bw() + theme_pubr(base_size = 10) + ylab('Proliferative index') + xlab("") + 
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