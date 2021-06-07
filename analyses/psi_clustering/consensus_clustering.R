################################################################################
# consensus_clustering.R
# script that takes in PSI data file and performs consensus clustering
# written by Ammar Naqvi
#
# usage: Rscript consensus_clustering.R
################################################################################

library("randomcoloR")
library("pheatmap")
library("ConsensusClusterPlus")
library("ggplot2")
library("dplyr")
library("vroom")
library("ggplot2")

suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("rlist"))

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`


## input psi matrix (removing duplicates, cell lines, and second malignancies)
dataDir = "/Users/naqvia/Desktop/pbta-splicing/analyses/psi_clustering/input/"
file_psi <- "/Users/naqvia/Desktop/pbta-splicing/analyses/psi_clustering/results/pan_cancer_splicing.thr10.report_select.remDup.v2.txt"

## temp for gene expr
#file_psi <- "/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/results/pan_cancer_expr_tpm_report_select.remDup.txt"
psi_tab  = read.delim(file_psi, sep = "\t", row.names=1, header=TRUE)

rnames <- psi_tab[,1]
row.names(psi_tab) <- psi_tab$Splice_ID
mat_hm <- data.matrix(psi_tab[,2:ncol(psi_tab)])

d=mat_hm

## reduce the dataset to the top 5% most variable genes, measured by median absolute deviation
mads=apply(d,1,mad)
d=d[rev(order(mads))[1:5417],] ## top 5% .05*108352

## the default settings of the agglomerative hierarchical clustering algorithm using Pearson correlation distance, so it is appropriate to gene median center d using
d = sweep(d,1, apply(d,1,median,na.rm=T))

## remove NAs
is.na(d) <- sapply(d, is.infinite)
d[is.na(d)] <- 0
d[is.nan(d)] <- 0

## k= 3 clusters pam+spearman
results = ConsensusClusterPlus((d),maxK=10,reps=100,pItem=0.8,
                     title="clustering",clusterAlg="pam",distance="spearman",seed=123,innerLinkage = "average", finalLinkage = "average")

# choose a cluster that seems best and assign to n_cluster
CC_group <- results[[3]]$consensusClass %>%
  as.data.frame()
colnames(CC_group) <- "CC"

write.table(CC_group, file = "/Users/naqvia/Desktop/pbta-splicing/analyses/psi_clustering/results/CC_groups_remDup.v2.txt", row.names=TRUE, sep="\t")



## run script to add clustering info to histology files // input/pbta-histologies_w_clusters.tsv // not working
#system("/Users/naqvia/Desktop/AS-DMG/analyses/psi_clustering/combine_clin_cluster.pl /Users/naqvia/Desktop/AS-DMG/analyses/psi_clustering/results/CC_groups_remDup.txt /Users/naqvia/Desktop/AS-DMG/analyses/psi_clustering/input/pbta-histologies.RNA-Seq.initial.tsv")

# read in consensus clustering matrix
CC_consensus_mat <- results[[3]]$consensusMatrix
colnames(CC_consensus_mat) <- rownames(CC_group)
rownames(CC_consensus_mat) <- rownames(CC_group)

clin_file = "/Users/naqvia/Desktop/pbta-splicing/analyses/psi_clustering/input/pbta-histologies_w_clusters.v2.tsv"
cluster_mem = read.delim(clin_file, sep = "\t", header=TRUE)

hist_sample <- cbind(data.frame(clin_tab$Kids_First_Biospecimen_ID), data.frame(clin_tab$short_histology), data.frame(clin_tab$molecular_subtype), data.frame(clin_tab$Cluster))

## specificy colors 
anno_palette <- cluster_mem$Color
hist_sample$clin_tab.short_histology <- as.factor(hist_sample$clin_tab.short_histology )
rownames(hist_sample)<- hist_sample$clin_tab.Kids_First_Biospecimen_ID

##remove colum
hist_sample = subset(hist_sample, select = -c(clin_tab.Kids_First_Biospecimen_ID))

pheatmap::pheatmap(
  CC_consensus_mat,
  #color = cluster_mem$Color,
  annotation_col=hist_sample,
  #annotation_col= cluster_mem$Color,
  #annotation_colors=anno_palette,
  cluster_rows = results[[3]]$consensusTree,
  cluster_cols = results[[3]]$consensusTree,
  show_rownames = F,
  show_colnames = F
)




theme_Publication <- function(base_size=15, base_family="Helvetica") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size, base_family=base_family)
      + theme(plot.title = element_text(face = "bold",
                                        size = rel(1.2), hjust = 0.5),
              text = element_text(),
              panel.background = element_rect(colour = NA),
              plot.background = element_rect(colour = NA),
              panel.border = element_rect(colour = NA),
              axis.title = element_text(face = "bold",size = rel(1)),
              axis.title.y = element_text(angle=90,vjust =2),
              axis.title.x = element_text(vjust = -0.2),
              axis.text = element_text(),
              axis.line = element_line(colour="black"),
              axis.ticks = element_line(),
              panel.grid.major = element_line(colour="#f0f0f0"),
              panel.grid.minor = element_blank(),
              legend.key = element_rect(colour = NA),
              legend.position = "right",
              legend.direction = "vertical",
              legend.key.size= unit(0.5, "cm"),
              # legend.margin = unit(0.5, "cm"),
              legend.margin = margin(5,5,5,5),
              legend.title = element_text(face="bold"),
              #plot.margin=unit(c(10,5,5,5),"mm"),
              plot.margin=unit(c(10,5,5,10),"mm"),
              strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
              strip.text = element_text(face="bold")
      ))
  }
