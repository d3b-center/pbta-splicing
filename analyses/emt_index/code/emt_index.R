# Author: Komal S. Rathi
# GSVA using EMT gene sets

# load libraries
library(msigdbr)
library(GSVA)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(reshape2)

# root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analyses_dir <- file.path(root_dir, 'analyses', 'emt_index')
output_dir <- file.path(analyses_dir, 'output')

# source plotting theme
source(file.path(root_dir, 'util', 'pubTheme.R'))

emt_index <- function(expr, meta, title, method){
  expr <- readRDS(expr)
  meta <- readRDS(meta)
  
  # log2 FPKM data
  expr <- log2(expr+1)
  
  # GSVA
  hallmark_geneset <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  emt_geneset <- hallmark_geneset %>% filter(gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") %>%
    .$human_gene_symbol 
  emt_geneset <- list(emt_geneset = emt_geneset)
  
  ssgsea_output <- gsva(as.matrix(expr), 
                        gset.idx.list = emt_geneset, 
                        method = "gsva",
                        min.sz = 1, max.sz = 1500,
                        mx.diff = TRUE)
  ssgsea_output <- melt(ssgsea_output, varnames = c("molecular_function", "sample"), value.name = "emt_index")
  ssgsea_output <- ssgsea_output %>%
    inner_join(meta, by = 'sample')
  
  # shapiro test of normality
  molecular_functions <- unique(ssgsea_output$molecular_function)
  for(i in 1:length(molecular_functions)){
    print(molecular_functions[i])
    print(shapiro.test(x = ssgsea_output[ssgsea_output$type == "pos" & ssgsea_output$molecular_function == molecular_functions[i], 'emt_index'])) 
    print(shapiro.test(x = ssgsea_output[ssgsea_output$type == "neg" & ssgsea_output$molecular_function == molecular_functions[i], 'emt_index'])) 
  }
  
  
  # plot
  my_comparisons <- list(c("pos", "neg"))
  p <- ggplot(ssgsea_output, aes(x = type, y = emt_index, fill = type)) + 
    facet_wrap(~molecular_function) +
    stat_boxplot(geom ='errorbar', width = 0.2) +
    geom_boxplot(lwd = 0.5, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) +
    scale_x_discrete(labels = c("pos" = "clk1_splice_pos", 
                                "neg" = "clk1_splice_neg")) +
    scale_fill_manual(values = c("pos" = "#F8766D",
                                 "neg" = "#00BFC4")) +
    theme_Publication() + xlab("") + ylab("EMT Index") +
    guides(fill = FALSE) + ggtitle(title) +
    stat_compare_means(label.y = 0.9, color = "darkred", paired = FALSE, method = method)
  return(p)
}

p <- emt_index(expr = file.path(root_dir, 'data', 'pbta-hgat-dx-prog-pm-gene-expression-subset-rsem-tpm-uncorrected.rds'), 
               meta = file.path(root_dir, 'data', 'pbta-hgat-dx-prog-pm-gene-expression-subset-metadata.rds'), 
               title = "EMT Index (nSamples = 8)", 
               method = "t.test") 
ggsave(filename = file.path(output_dir, 'emt_index_plot.pdf'), device = "pdf", plot = p,width = 10, height = 7)



