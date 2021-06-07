# Author: Komal S. Rathi
# GSVA using Oxidative Stress genes

# load libraries
library(msigdbr)
library(GSVA)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(reshape2)

# root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analyses_dir <- file.path(root_dir, 'analyses', 'oxidative_stress_index')
output_dir <- file.path(analyses_dir, 'output')

# source plotting theme
source(file.path(root_dir, 'util', 'pubTheme.R'))

oxidative_stress_index <- function(expr, meta, title, method){
  expr <- readRDS(expr)
  meta <- readRDS(meta)
  
  # log2 FPKM data
  expr <- log2(expr+1)

  # GSVA
  hallmark_geneset <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
  oxidative_stress_geneset <- hallmark_geneset %>% 
    filter(gs_name %in% c("GO_RESPONSE_TO_OXIDATIVE_STRESS",
                          "GO_RESPONSE_TO_HYDROGEN_PEROXIDE",
                          "GO_CELLULAR_RESPONSE_TO_HYDROGEN_PEROXIDE"))
  oxidative_stress_geneset <- unstack(oxidative_stress_geneset[,c("human_gene_symbol", "gs_name")])
  
  ssgsea_output <- gsva(as.matrix(expr), 
                        gset.idx.list = oxidative_stress_geneset, 
                        method = "gsva",
                        min.sz = 1, max.sz = 1500,
                        mx.diff = TRUE)
  ssgsea_output <- melt(ssgsea_output, varnames = c("molecular_function", "sample"), value.name = "oxidative_stress_index")
  ssgsea_output <- ssgsea_output %>%
    inner_join(meta, by = 'sample')
  
  # shapiro test of normality
  molecular_functions <- unique(ssgsea_output$molecular_function)
  for(i in 1:length(molecular_functions)){
    print(shapiro.test(x = ssgsea_output[ssgsea_output$type == "pos" & ssgsea_output$molecular_function == molecular_functions[i], 'oxidative_stress_index'])) 
    print(shapiro.test(x = ssgsea_output[ssgsea_output$type == "neg" & ssgsea_output$molecular_function == molecular_functions[i], 'oxidative_stress_index'])) 
  }
  
  
  # plot
  my_comparisons <- list(c("pos", "neg"))
  p <- ggplot(ssgsea_output, aes(x = type, y = oxidative_stress_index, fill = type)) + 
    facet_wrap(~molecular_function) +
    stat_boxplot(geom ='errorbar', width = 0.2) +
    geom_boxplot(lwd = 0.5, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) +
    scale_x_discrete(labels = c("pos" = "clk1_splice_pos", 
                                "neg" = "clk1_splice_neg")) +
    scale_fill_manual(values = c("pos" = "#F8766D",
                                 "neg" = "#00BFC4")) +
    theme_Publication() + xlab("") + ylab("Oxidative Stress Index") +
    guides(fill = FALSE) + ggtitle(title) +
    stat_compare_means(label.y = 0.9, color = "darkred", paired = FALSE, method = method)
  return(p)
}

p <- oxidative_stress_index(expr = file.path(root_dir, 'data', 'pbta-hgat-dx-prog-pm-gene-expression-subset-rsem-tpm-uncorrected.rds'), 
                            meta = file.path(root_dir, 'data', 'pbta-hgat-dx-prog-pm-gene-expression-subset-metadata.rds'), 
                            title = "Oxidative Stress Index (nSamples = 8)", method = "t.test") 
ggsave(filename = file.path(output_dir, 'oxidative_stress_index_plot.pdf'), device = "pdf", plot = p,width = 10, height = 7)
