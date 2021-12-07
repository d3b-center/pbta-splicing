# function to calculate oxidative stress index and create boxplot
suppressPackageStartupMessages({
  library(msigdbr)
  library(GSVA)
  library(tidyverse)
  library(ggpubr)
  library(ggplot2)
  library(reshape2)
})

oxidative_stress_index <- function(expr, meta, method){
  expr <- readRDS(expr)
  meta <- readRDS(meta)
  
  # log2 FPKM data
  expr <- log2(expr+1)
  
  # GSVA
  hallmark_geneset <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
  oxidative_stress_geneset <- hallmark_geneset %>% 
    filter(gs_name %in% c("GOBP_RESPONSE_TO_OXIDATIVE_STRESS",
                          "GOBP_RESPONSE_TO_HYDROGEN_PEROXIDE",
                          "GOBP_CELLULAR_RESPONSE_TO_HYDROGEN_PEROXIDE"))
  oxidative_stress_geneset <- unstack(oxidative_stress_geneset[,c("human_gene_symbol", "gs_name")])
  
  ssgsea_output <- gsva(as.matrix(expr), 
                        gset.idx.list = oxidative_stress_geneset, 
                        method = "gsva",
                        min.sz = 1, max.sz = 1500,
                        mx.diff = TRUE)
  ssgsea_output <- melt(ssgsea_output, varnames = c("molecular_function", "Kids_First_Biospecimen_ID"), value.name = "oxidative_stress_index")
  ssgsea_output <- ssgsea_output %>%
    inner_join(meta, by = 'Kids_First_Biospecimen_ID')
  
  # plot
  p <- ggplot(ssgsea_output, aes(x = type, y = oxidative_stress_index, fill = type)) +
    stat_boxplot(geom ='errorbar', width = 0.2) +
    geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.color = "white",
                 lwd = 0.5, fatten = 0.7, width = 0.5) + 
    ggtitle('Oxidative stress index') +
    theme_bw() + theme_pubr(base_size = 10) + ylab('Oxidative stress index') + xlab("") + 
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