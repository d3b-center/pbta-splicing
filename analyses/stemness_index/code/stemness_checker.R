# Author: Komal S. Rathi
# Function: Run OCLR based Stemness profiling
# References: 
# Paper: https://www.cell.com/cell/pdf/S0092-8674(18)30358-1.pdf
# Tutorial: http://tcgabiolinks.fmrp.usp.br/PanCanStem/mRNAsi.html

# load libraries
# use synapser instead of synapseClient
library(synapser) 
library(gelnet)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(reshape2)

# Synapse login required to get PCBC data
readRenviron('~/.Renviron')
SYNAPSE_ID <- Sys.getenv("SYNAPSE_ID")
SYNAPSE_PWD <- Sys.getenv("SYNAPSE_PWD")
synLogin(email = SYNAPSE_ID, password = SYNAPSE_PWD) # your id and password from https://www.synapse.org/

# root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analyses_dir <- file.path(root_dir, 'analyses', 'stemness_index')
output_dir <- file.path(analyses_dir, 'output')

# source plotting theme
source(file.path(root_dir, 'util', 'pubTheme.R'))

# source scripts for training and prediction
source(file.path(analyses_dir, 'util', 'genes2hugo.R')) # convert entrez ids to hugo symbols
source(file.path(analyses_dir, 'util', 'main.train.R')) # train using PCBC dataset
source(file.path(analyses_dir, 'util', 'main.predict.R')) # predict using vector generated in main.train.R

# run training set (the output of this should be vector of length 78 containing only 1s)
# stem/progenitor cells from the Progenitor Cell Biology Consortium
fname <- file.path(analyses_dir, 'input', 'pcbc-stemsig.tsv')
if(file.exists(fname)){
  print("Training set ready")
} else {
  print("Run training set")
  auc <- main.train(fnOut = fname, fnGenes = NULL)
}

# function to create stemness index boxplot
stemness_plot <- function(expr, meta, fnSig, fnOut, facet = NULL, method = c('t.test')){
  
  # use predict function
  main.predict(fnSig = fnSig, expr = expr, fnOut = fnOut)
  
  # load data
  meta <- readRDS(meta)
  
  # visualization
  predicted.stemness.scores <- read.delim(fnOut, header = F)
  colnames(predicted.stemness.scores) <- c("sample_name", "stemness_score")
  predicted.stemness.scores <- predicted.stemness.scores %>%
    inner_join(meta, by = c("sample_name" = "sample"))
  
  # shapiro test of normality
  print(shapiro.test(x = predicted.stemness.scores[predicted.stemness.scores$type == "pos", 'stemness_score']))
  print(shapiro.test(x = predicted.stemness.scores[predicted.stemness.scores$type == "neg", 'stemness_score']))
  
  # let's use Kruskal Wallis here as type a is not normally distributed
  predicted.stemness.scores$type <- factor(predicted.stemness.scores$type, levels = c("pos", "neg"))
  my_comparisons <- list(c("pos", "neg"))
  if(!is.null(facet)){
    p <- ggplot(predicted.stemness.scores, aes(x = type, y = stemness_score, fill = type)) +
      facet_wrap(~library_type) +
      stat_boxplot(geom ='errorbar', width = 0.2) +
      geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.color = "white",
                   lwd = 0.5, fatten = 0.7, width = 0.5) + 
      ggtitle('OCLR predicted Stemness index\nCLK1 splice-positive vs splice-negative') +
      theme_bw() + theme_Publication(base_size = 10) + ylab('Stemness Index') + xlab("") + 
      geom_text(aes(label = sample_name), size = 3, hjust = 0.7, vjust = -0.5) +
      guides(fill = F) + 
      geom_jitter(position=position_jitter(width=.01), shape = 21) +
      scale_x_discrete(labels = c("pos" = "clk1_splice_pos", 
                                  "neg" = "clk1_splice_neg")) +
      scale_fill_manual(values = c("pos" = "#F8766D",
                                   "neg" = "#00BFC4")) +
      stat_compare_means(comparisons = my_comparisons, color = "darkred", size = 3, paired = FALSE, method = method) +
      stat_compare_means(label.y = 1.4, color = "darkred", paired = FALSE, method = method) +
      labs(fill = "Type")
  } else {
    p <- ggplot(predicted.stemness.scores, aes(x = type, y = stemness_score, fill = type)) +
      stat_boxplot(geom ='errorbar', width = 0.2) +
      geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.color = "white",
                   lwd = 0.5, fatten = 0.7, width = 0.5) + 
      ggtitle('OCLR predicted Stemness index (nSamples = 8)') +
      theme_bw() + theme_Publication(base_size = 10) + ylab('Stemness Index') + xlab("") + 
      geom_text(aes(label = sample_name), size = 3, hjust = 0.7, vjust = -0.5) +
      guides(fill = F) + 
      geom_jitter(position=position_jitter(width=.01), shape = 21) +
      scale_x_discrete(labels = c("pos" = "clk1_splice_pos", 
                                  "neg" = "clk1_splice_neg")) +
      scale_fill_manual(values = c("pos" = "#F8766D",
                                   "neg" = "#00BFC4")) +
      stat_compare_means(comparisons = my_comparisons, color = "darkred", size = 3, paired = FALSE, method = method) +
      stat_compare_means(label.y = 1.4, color = "darkred", paired = FALSE, method = method) +
      labs(fill = "Type")
  }
  return(p)
}

# list (n = 8)
p <- stemness_plot(expr = file.path(root_dir, 'data', 'pbta-hgat-dx-prog-pm-gene-expression-subset-rsem-tpm-uncorrected.rds'), 
                   meta = file.path(root_dir, 'data', 'pbta-hgat-dx-prog-pm-gene-expression-subset-metadata.rds'), 
                   fnSig = fname,
                   fnOut = file.path(output_dir, 'stemness_scores.tsv'), 
                   facet = NULL, method = "t.test")
ggsave(filename = file.path(output_dir, 'stemness_index_plot.pdf'), device = "pdf", plot = p,width = 10, height = 7)
