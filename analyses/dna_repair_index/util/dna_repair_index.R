# function to calculate DNA repair index and create boxplot
suppressPackageStartupMessages({
  library(msigdbr)
  library(GSVA)
  library(tidyverse)
  library(ggpubr)
  library(ggplot2)
  library(reshape2)
})

dna_repair_index <- function(expr, meta, method){
  expr <- readRDS(expr)
  meta <- readRDS(meta)
  
  # log2 FPKM data
  expr <- log2(expr+1)
  
  # GSVA
  dna_repair_geneset <- msigdbr::msigdbr(species = "Homo sapiens", 
                                         category = "C2", 
                                         subcategory = "CP:KEGG")
  dna_repair_geneset <- dna_repair_geneset %>%
    filter(gs_name %in% c("KEGG_MISMATCH_REPAIR",
                          "KEGG_HOMOLOGOUS_RECOMBINATION",
                          "KEGG_NON_HOMOLOGOUS_END_JOINING",
                          "KEGG_BASE_EXCISION_REPAIR"))
  dna_repair_geneset <- unstack(dna_repair_geneset[,c("human_gene_symbol", "gs_name")])
  
  ssgsea_output <- gsva(as.matrix(expr), 
                        gset.idx.list = dna_repair_geneset, 
                        method = "gsva",
                        min.sz = 1, max.sz = 1500,
                        mx.diff = TRUE)
  ssgsea_output <- melt(ssgsea_output, varnames = c("molecular_function", "Kids_First_Biospecimen_ID"), value.name = "dna_repair_index")
  ssgsea_output <- ssgsea_output %>%
    inner_join(meta, by = 'Kids_First_Biospecimen_ID')
  
  # plot
  p <- ggplot(ssgsea_output, aes(x = type, y = dna_repair_index, fill = type)) +
    stat_boxplot(geom ='errorbar', width = 0.2) +
    geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.color = "white",
                 lwd = 0.5, fatten = 0.7, width = 0.5) + 
    ggtitle('DNA repair index') +
    theme_bw() + theme_pubr(base_size = 10) + ylab('DNA repair index') + xlab("") + 
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