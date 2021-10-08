################################################################################
# splicing_functional_sites.R
# script that plots a table of splicing events overlapping uniprot sites
# written by Ammar Naqvi
#
# usage: Rscript splicing_functional_sites.R <file>
################################################################################

suppressPackageStartupMessages(library(ggthemes))
suppressPackageStartupMessages(library(extrafont))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggstatsplot))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## get command line arg -- tsv 
args <- commandArgs(trailing = TRUE)
file_pos = args[1]
file_neg = args[2]

#print (file_neg)
## hardcoded for now
## table with splicing events and associated uniprot-annotaed site 
## table with splicing events and associated uniprot-annotaed site dvided by negative and positive
#file_pos = "/Users/naqvia/Desktop/pbta-splicing/analyses/splicing_events_functional_sites/results/splicing_events.total.pos.intersectUnip.ggplot.txt" 
#file_neg = "/Users/naqvia/Desktop/pbta-splicing/analyses/splicing_events_functional_sites/results/splicing_events.total.neg.intersectUnip.ggplot.txt" 

dpsi_unip_pos <- read.table(file_pos, header=TRUE,sep = "\t")
dpsi_unip_neg <- read.table(file_neg, header=TRUE,sep = "\t")

set.seed(123)
plot1 <- ggstatsplot::ggbetweenstats(
  data = dpsi_unip_neg, 
  x = Uniprot, 
  y = dPSI,
  k = 3,
  nboot = 15,
  outlier.label = SpliceID, # label to attach to outlier values
  outlier.label.args = list(color = "red"), # outlier point label color
  notch = TRUE,
  mean.ci = TRUE,
  outlier.tagging = TRUE,
  #title = "Tumor supressors",
  type = "robust",
  xlab = "Unipro-defined Site",
  pairwise.comparisons = TRUE,
  messages = FALSE
) + theme_Publication()

set.seed(123)
plot2 <- ggstatsplot::ggbetweenstats(
  data = dpsi_unip_pos, 
  x = Uniprot, 
  y = dPSI,
  k = 3,
  nboot = 15,
  outlier.label = SpliceID, # label to attach to outlier values
  outlier.label.args = list(color = "red"), # outlier point label color
  notch = TRUE,
  mean.ci = TRUE,
  #outlier.tagging = TRUE,
  type = "robust",
  xlab = "Unipro-defined Site",
  pairwise.comparisons = TRUE,
  messages = FALSE
) + theme_Publication()

setwd("/Users/naqvia/Desktop/pbta-splicing/analyses/splicing_events_functional_sites")
ggsave("plots/dPSI_across_functional_sites_pos.pdf", width = 15, height = 5)
ggsave("plots/dPSI_across_functional_sites_neg.pdf", width = 15, height = 5)


# ##fitlered by given gene lists
#kinase
gene_list_kinase <- read.delim("/Users/naqvia/Desktop/create_oncoprint/data/kinase_list.txt", sep = "\t",
                               header = F, as.is = T)

#brain cancer genes
gene_list_brain <- read.delim("/Users/naqvia/Desktop/create_oncoprint/data/brain-goi-list-new.txt", sep = "\t",
                              header = F, as.is = T) 

#splicing factors
gene_list_sfs <- read.delim("/Users/naqvia/Desktop/create_oncoprint/data/splicing_factors.txt", sep = "\t",
                            header = F, as.is = T)

#SWI/SNF genes
gene_list_swisnf <- read.delim("/Users/naqvia/Desktop/create_oncoprint/data/swi_snf_genes.txt", sep = "\t",
                               header = F, as.is = T)

#epigenetic related genes
gene_list_epi <- read.delim("/Users/naqvia/Desktop/create_oncoprint/data/epi_genes_list.grep.txt", sep = "\t",
                            header = F, as.is = T) 

#tf related genes
gene_list_tf <- read.delim("/Users/naqvia/Desktop/create_oncoprint/data/tf_gene_list.txt", sep = "\t",
                           header = F, as.is = T) 
# ## brain cancer
# dpsi_unip_neg_brain_goi_filter = dpsi_unip_neg[dpsi_unip_neg$Gene %in% gene_list_brain$V1,]
# dpsi_unip_pos_brain_goi_filter = dpsi_unip_pos[dpsi_unip_pos$Gene %in% gene_list_brain$V1,]
# 
# ##kinases
# dpsi_unip_neg_kinase_filter = dpsi_unip_neg[dpsi_unip_neg$Gene %in% gene_list_kinase$V1,]
# dpsi_unip_pos_kinase_filter = dpsi_unip_pos[dpsi_unip_pos$Gene %in% gene_list_kinase$V1,]
# 
# ##epigenetic
# dpsi_unip_neg_epi_filter = dpsi_unip_neg[dpsi_unip_neg$Gene %in% gene_list_epi$V1,]
# dpsi_unip_pos_epi_filter = dpsi_unip_pos[dpsi_unip_pos$Gene %in% gene_list_epi$V1,]
# 
# ##tfs
# dpsi_unip_neg_tf_filter = dpsi_unip_neg[dpsi_unip_neg$Gene %in% gene_list_tf$V1,]
# dpsi_unip_pos_tf_filter = dpsi_unip_pos[dpsi_unip_pos$Gene %in% gene_list_tf$V1,]
# 
# ##sfs
# dpsi_unip_neg_sf_filter = dpsi_unip_neg[dpsi_unip_neg$Gene %in% gene_list_sfs$V1,]
# dpsi_unip_pos_sf_filter = dpsi_unip_pos[dpsi_unip_pos$Gene %in% gene_list_sfs$V1,]
# 
# ##plots
# set.seed(123)
# ggstatsplot::ggbetweenstats(
#   data = dpsi_unip_neg_brain_goi_filter, 
#   x = Uniprot, 
#   y = dPSI,
#   k = 3,
#   nboot = 15,
#   outlier.label = SpliceID, # label to attach to outlier values
#   outlier.label.args = list(color = "red"), # outlier point label color
#   notch = TRUE,
#   mean.ci = TRUE,
#   outlier.tagging = TRUE,
#   title = "Cancer Genes Splicing",
#   type = "robust",
#   xlab = "Unipro-defined Site",
#   pairwise.comparisons = TRUE,
#   messages = FALSE
# ) + theme_Publication()
# 
# set.seed(123)
# ggstatsplot::ggbetweenstats(
#   data = dpsi_unip_pos_brain_goi_filter, 
#   x = Uniprot, 
#   y = dPSI,
#   k = 3,
#   nboot = 15,
#   outlier.label = SpliceID, # label to attach to outlier values
#   outlier.label.args = list(color = "red"), # outlier point label color
#   notch = TRUE,
#   mean.ci = TRUE,
#   outlier.tagging = TRUE,
#   title = "Cancer Genes Splicing",
#   type = "robust",
#   xlab = "Unipro-defined Site",
#   pairwise.comparisons = TRUE,
#   messages = FALSE
# ) + theme_Publication()
# 
# ##
# set.seed(123)
# ggstatsplot::ggbetweenstats(
#   data = dpsi_unip_pos_epi_filter, 
#   x = Uniprot, 
#   y = dPSI,
#   k = 3,
#   nboot = 15,
#   outlier.label = SpliceID, # label to attach to outlier values
#   outlier.label.args = list(color = "red"), # outlier point label color
#   notch = TRUE,
#   mean.ci = TRUE,
#   outlier.tagging = TRUE,
#   title = "Cancer Genes Splicing",
#   type = "robust",
#   xlab = "Unipro-defined Site",
#   pairwise.comparisons = TRUE,
#   messages = FALSE
# ) + theme_Publication()
# 
# ## kinases
# set.seed(123)
# ggstatsplot::ggbetweenstats(
#   data = dpsi_unip_pos_kinase_filter, 
#   x = Uniprot, 
#   y = dPSI,
#   k = 3,
#   nboot = 15,
#   outlier.label = SpliceID, # label to attach to outlier values
#   outlier.label.args = list(color = "red"), # outlier point label color
#   #notch = TRUE,
#   #mean.ci = TRUE,
#   #outlier.tagging = TRUE,
#   title = "Kinase Splicing",
#   type = "robust",
#   xlab = "Unipro-defined Site",
#   pairwise.comparisons = TRUE,
#   messages = FALSE
# ) + theme_Publication()
# 
# JitterPlot = qplot(Uniprot, dPSI, data=dpsi_unip_neg_kinase_filter, geom=c("boxplot", "jitter"), 
#                    main="Kinase Aberrant Splicing",xlab="Type", ylab="dPSI") 
# 
# JitterPlot
# 
# set.seed(123)
# ggstatsplot::ggbetweenstats(
#   data = dpsi_unip_neg_kinase_filter, 
#   x = Uniprot, 
#   y = dPSI,
#   k = 3,
#   nboot = 15,
#   outlier.label = SpliceID, # label to attach to outlier values
#   outlier.label.args = list(color = "red"), # outlier point label color
#   #notch = TRUE,
#   #mean.ci = TRUE,
#   #outlier.tagging = TRUE,
#   title = "Kinase Splicing",
#   type = "robust",
#   xlab = "Unipro-defined Site",
#   pairwise.comparisons = TRUE,
#   messages = FALSE
# ) + theme_Publication()
# 


