################################################################################
# upsetR_of_hist.R
# script that takes in "results/CC_memberships_expr" and CC_memberships_psi data 
#files and computes cluster membership/groupings overlap
# written by Ammar Naqvi
#
# usage: Rscript upsetR_of_clusters.R
################################################################################

library(UpSetR)
library(ggplot2)
library(gplots)

list_1 = read.table("/Users/naqvia/Desktop/pbta-splicing/analyses/histology_specific/results/perc_hist_as.30prev.ATRT.tsv", sep = "\t")
list_2 = read.table("/Users/naqvia/Desktop/pbta-splicing/analyses/histology_specific/results/perc_hist_as.30prev.Cran.tsv", sep = "\t")
list_3 = read.table("/Users/naqvia/Desktop/pbta-splicing/analyses/histology_specific/results/perc_hist_as.30prev.Epend.tsv", sep = "\t")
list_4 = read.table("/Users/naqvia/Desktop/pbta-splicing/analyses/histology_specific/results/perc_hist_as.30prev.HGAT.tsv", sep = "\t")
list_5 = read.table("/Users/naqvia/Desktop/pbta-splicing/analyses/histology_specific/results/perc_hist_as.30prev.LGAT.tsv", sep = "\t")
list_6 = read.table("/Users/naqvia/Desktop/pbta-splicing/analyses/histology_specific/results/perc_hist_as.30prev.medul.tsv", sep = "\t")
list_7 = read.table("/Users/naqvia/Desktop/pbta-splicing/analyses/histology_specific/results/perc_hist_as.30prev.Gang.tsv",sep = "\t")

listInput <- list("ATRT" =list_1$V2, 
                  "Cran" =list_2$V2,
                  "Epend" =list_3$V2,
                  "HGAT" =list_4$V2,
                  "LGAT" =list_5$V2, 
                  "Med" =list_6$V2,
                  "Gang" =list_7$V2)

upset(fromList(listInput), 
      mainbar.y.label = "", sets.x.label = "Clusters", order.by = "freq",
      mb.ratio = c(0.5,0.50), text.scale = c(1.3, 1.3, 1.3, 1.3, 2, 1.4),point.size = 2, line.size = 1.5, nsets = 7)
