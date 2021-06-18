################################################################################
# stacked_barplots.R
# script that takes in cluster info file and plots stacked barplot
# written by Ammar Naqvi
#
# usage: Rscript stacked_barplots.R
################################################################################

library("randomcoloR")
library("pheatmap")
library("ConsensusClusterPlus")
library("ggplot2")
library("dplyr")
library("vroom")
library("ggplot2")

file = "~/Desktop/cluster_hist_members.RO3.txt"
cluster_mem      = read.delim(file, sep = "\t",header=TRUE)


##stacked barplot of cluster memberships
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


cluster_mem_file = "/Users/naqvia/Desktop/pbta-splicing/analyses/psi_clustering/input/cluster_hist_members.txt"
cluster_mem      = read.delim(cluster_mem_file, sep = "\t",header=TRUE)

ggplot(cluster_mem, aes(fill=short_histology, x=Cluster)) +
  geom_bar(stat="count", position="stack") + scale_fill_manual(values=c(cbbPalette)) + theme_Publication()

##stacked barplot of cluster memberships of HGATs
cluster_mem_hgat_file = "/Users/naqvia/Desktop/pbta-splicing/analyses/psi_clustering/results/cluster_hist_members.hgats.txt"
cluster_mem_hgat      = read.delim(cluster_mem_hgat_file, sep = "\t",header=TRUE)

# Stacked
cbbPalette_hgg <- c("#800000", "#8B0000", "#A52A2A", "#B22222", "#DC143C", "#FF0000", "#FF6347", "#FF7F50","#CD5C5C","#F08080","#E9967A","#FA8072");

ggplot(cluster_mem_hgat, aes(fill=Subtype, x=Cluster)) +
  geom_bar(stat="count", position="stack")  + theme_Publication()
#+ scale_fill_manual(values=c(cbbPalette_hgg))

cluster_mem_hgat %>%
  ggplot() +
  geom_bar(aes(x = Cluster, fill = Subtype)) +
  scale_fill_manual(values=c(cbbPalette_hgg)) +
  theme_Publication()


cluster_mem_lgat_file = "/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/results/cluster_hist_makeup.lgats.v3.txt"
cluster_mem_lgat      = read.delim(cluster_mem_lgat_file, sep = "\t",header=TRUE)

# Stacked
cbbPalette_lgg <- c("#191970", "#000080", "#00008B", "#0000CD", "#0000FF", "#4169E1", "#87CEFA", "#87CEEB","#ADD8E6","#1E90FF","#00BFFF","#6495ED", "#4682B4", "#5F9EA0");

ggplot(cluster_mem_lgat, aes(fill=Subtype, x=Cluster)) +
  geom_bar(position="stack", stat="count") + scale_fill_manual(values=c(cbbPalette_lgg)) + theme_Publication()

cluster_mem_lgat %>%
  ggplot() +
  geom_bar(aes(x = Cluster, fill = Subtype)) +
  scale_fill_manual(values=c(cbbPalette_lgg)) +
  theme_Publication()
