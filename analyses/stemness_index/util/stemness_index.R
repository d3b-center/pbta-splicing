# function to calculate stemness index and create boxplot
suppressPackageStartupMessages({
  library(ggpubr)
  library(ggplot2)
})
stemness_index <- function(expr, meta, fnSig, fnOut, method = c('t.test')){
  
  # use predict function
  main.predict(fnSig = fnSig, expr = expr, fnOut = fnOut)
  
  # load data
  meta <- readRDS(meta)
  
  # merge with metadata
  predicted_stemness_scores <- read.delim(fnOut, header = F)
  colnames(predicted_stemness_scores) <- c("Kids_First_Biospecimen_ID", "stemness_score")
  predicted_stemness_scores <- predicted_stemness_scores %>%
    inner_join(meta, by = "Kids_First_Biospecimen_ID")
  
  # plot
  p <- ggplot(predicted_stemness_scores, aes(x = type, y = stemness_score, fill = type)) +
    stat_boxplot(geom ='errorbar', width = 0.2) +
    geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.color = "white",
                 lwd = 0.5, fatten = 0.7, width = 0.5) + 
    ggtitle('OCLR predicted Stemness index') +
    theme_bw() + theme_pubr(base_size = 10) + ylab('Stemness Index') + xlab("") + 
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
