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
dataDir = "/Users/naqvia/Desktop/pbta-splicing/analyses/psi_clustering/results/"
file_psi <- "/Users/naqvia/Desktop/pbta-splicing/analyses/psi_clustering/results/pan_cancer_splicing.thr10.report_select.remDup.v2.txt"

df = read.table(args[1], header=TRUE)

print(df)


## temp for gene expr
#file_psi <- "/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/results/pan_cancer_expr_tpm_report_select.remDup.txt"
psi_tab  = read.delim(file_psi, sep = "\t", header=TRUE)

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
colnames(CC_group) <- "Cluster"


## write cluster file for vtest
cluster_tab <- rownames_to_column(CC_group,var="Kids_First_Biospecimen_ID")
write.table(cluster_tab, file = "/Users/naqvia/Desktop/pbta-splicing/analyses/psi_clustering/results/CC_groups.txt", quote=FALSE,row.names=FALSE, sep="\t")

## run script to add clustering info to histology files // input/pbta-histologies_w_clusters.tsv // not working
# system("/Users/naqvia/Desktop/AS-DMG/analyses/psi_clustering/combine_clin_cluster.pl /Users/naqvia/Desktop/AS-DMG/analyses/psi_clustering/results/CC_groups_remDup.txt /Users/naqvia/Desktop/AS-DMG/analyses/psi_clustering/input/pbta-histologies.RNA-Seq.initial.tsv")

# read in consensus clustering matrix
CC_consensus_mat <- results[[3]]$consensusMatrix
colnames(CC_consensus_mat) <- rownames(CC_group)
rownames(CC_consensus_mat) <- rownames(CC_group)

clin_file = "/Users/naqvia/Desktop/pbta-splicing/analyses/psi_clustering/input/pbta-histologies.RNA-Seq.initial.tsv"
clin_tab = read.delim(clin_file, sep = "\t", header=TRUE)

#clin_file = "/Users/naqvia/Desktop/AS-DMG/analyses/psi_clustering/input/pbta-histologies.tsv"
#clin_tab = read.delim(clin_file, sep = "\t", header=TRUE)

## add cluster membership info for BS IDS
clin_tab <- clin_tab %>% left_join(rownames_to_column(CC_group,var="Kids_First_Biospecimen_ID"),by="Kids_First_Biospecimen_ID")

## make table with ID, Short histology and cluster info
hist_sample <- cbind(data.frame(clin_tab$Kids_First_Biospecimen_ID), data.frame(clin_tab$short_histology),data.frame(clin_tab$Cluster))

## specificy colors // not working
#anno_palette <- cluster_mem$Color

## convert to factors
hist_sample$clin_tab.short_histology <- as.factor(hist_sample$clin_tab.short_histology )
hist_sample$clin_tab.Cluster <- as.factor(hist_sample$clin_tab.Cluster )

rownames(hist_sample)<- hist_sample$clin_tab.Kids_First_Biospecimen_ID

##remove colum
hist_sample = subset(hist_sample, select = -c(clin_tab.Kids_First_Biospecimen_ID))

setwd("/Users/naqvia/Desktop/pbta-splicing/analyses/psi_clustering")
pheatmap::pheatmap(
  CC_consensus_mat,
  annotation_col=hist_sample,
  annotation_colors=anno_palette,
  cluster_rows = results[[3]]$consensusTree,
  cluster_cols = results[[3]]$consensusTree,
  show_rownames = F,
  show_colnames = F,
  filename = "plots/CC_heatmap.png"
)