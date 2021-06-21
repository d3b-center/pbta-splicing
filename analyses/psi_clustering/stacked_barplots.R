################################################################################
# stacked_barplots.R
# script that takes in cluster info file and plots stacked barplot
# written by Ammar Naqvi
#
# usage: Rscript stacked_barplots.R
################################################################################

library("pheatmap")
library("ConsensusClusterPlus")
library("ggplot2")
library("dplyr")
library("vroom")
library("ggplot2")
library("tidyverse")

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## get command line arg -- file
args <- commandArgs(trailing = TRUE)

##cluster file
#CC_file = "/Users/naqvia/Desktop/pbta-splicing/analyses/psi_clustering/results/CC_groups.txt"
#cluster_mem      = read.delim(CC_file, sep = "\t",header=TRUE)

## clinical histology file
#clin_file = "/Users/naqvia/Desktop/pbta-splicing/analyses/psi_clustering/input/pbta-histologies.RNA-Seq.initial.tsv"
#clin_tab = read.delim(clin_file, sep = "\t", header=TRUE)

## get cluster and histology files and make tables
setwd(".")
CC_file   <- args[1]
clin_file <- args[2]

cluster_mem  = read.delim(CC_file, sep = "\t",header=TRUE)
clin_tab     = read.delim(clin_file, sep = "\t", header=TRUE)

## combine histologies and cluster info
clin_tab    <- clin_tab %>% left_join(cluster_mem, by="Kids_First_Biospecimen_ID")
cluster_mem <- cbind(data.frame(clin_tab$short_histology), data.frame(clin_tab$Cluster))

##make plot
ggplot(cluster_mem, aes(fill=clin_tab.short_histology, x=clin_tab.Cluster)) +
  geom_bar(stat="count", position="stack") # + scale_fill_manual(values=c(cbbPalette)) + theme_Publication()

plot_barplot <- ggplot(cluster_mem, aes(fill=clin_tab.short_histology, x=clin_tab.Cluster)) +
  geom_bar(stat="count", position="stack") #+ scale_fill_manual(values=c(cbbPalette)) + theme_Publication()

png("plots/stacked_barplot_clusters.v2.png")
print(plot_barplot)
dev.off()
