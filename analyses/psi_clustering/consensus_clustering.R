################################################################################
# consensus_clustering.R
# script that takes in PSI data file and performs consensus clustering
# written by Ammar Naqvi
#
# usage: Rscript consensus_clustering.R
################################################################################

suppressPackageStartupMessages({
  library("pheatmap")
  library("ConsensusClusterPlus")
  library("ggplot2")
  library("dplyr")
  library("tidyverse")
  library("optparse")
})

##directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "psi_clustering")

results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")
input_dir <- file.path(analysis_dir, "input")

plots_dir <- file.path(analysis_dir, "plots")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

##download file from https://figshare.com/s/47bd539da8e9887143c8 ## too large for github
file_psi <- "pan_cancer_splicing_SE.txt"
psi_tab <- readr::read_tsv(file.path(input_dir, file_psi))
d <- psi_tab %>%
  as.data.frame() %>%
  tibble::column_to_rownames("Splice_ID")

## reduce the dataset to the top 5% most variable genes, measured by median absolute deviation
mads <- apply(d,1,mad)
top_5_perc <- round(length(mads)*0.05) ## top 5% .05*108352 (vs. 199134)
d <- d[rev(order(mads))[1:top_5_perc],] 

## the default settings of the agglomerative hierarchical clustering algorithm using Pearson correlation distance, so it is appropriate to gene median center d using
d <- sweep(d,1, apply(d,1,median,na.rm=T))

## remove NAs
is.na(d) <- sapply(d, is.infinite)
d[is.na(d)] <- 0

## k= 3 clusters pam+spearman after visual inspection
results <- ConsensusClusterPlus(as.matrix(d),
                                maxK=10,
                                reps=100,
                                pItem=0.8,
                                title="clustering",
                                clusterAlg="pam",
                                distance="spearman",
                                seed=123,
                                innerLinkage = "average", 
                                finalLinkage = "average")

# choose a cluster that seems best and assign to n_cluster
CC_group <- results[[3]]$consensusClass %>%
  as.data.frame()
colnames(CC_group) <- "Cluster"

## write cluster file for vtest 
cluster_tab <- tibble::rownames_to_column(CC_group, 
                                          var="Kids_First_Biospecimen_ID")
readr::write_tsv(cluster_tab, file = file.path(results_dir, 
                                               "CC_groups.tsv"))

# read in consensus clustering matrix
CC_consensus_mat <- results[[3]]$consensusMatrix
colnames(CC_consensus_mat) <- rownames(CC_group)
rownames(CC_consensus_mat) <- rownames(CC_group)

# read in clinical file
clin_tab <- readr::read_tsv(file.path(data_dir, "v19_plus_20210311_pnoc_rna.tsv")) %>%
  ## add cluster membership info for BS IDS
  dplyr::left_join(tibble::rownames_to_column(CC_group,var="Kids_First_Biospecimen_ID"),
                   by="Kids_First_Biospecimen_ID")

## make table with ID, Short histology and cluster info
hist_sample <- cbind(data.frame(clin_tab$Kids_First_Biospecimen_ID), 
                     data.frame(clin_tab$short_histology),
                     data.frame(clin_tab$Cluster))

## convert to factors
hist_sample$clin_tab.short_histology <- as.factor(hist_sample$clin_tab.short_histology )
hist_sample$clin_tab.Cluster <- as.factor(hist_sample$clin_tab.Cluster )

rownames(hist_sample)<- hist_sample$clin_tab.Kids_First_Biospecimen_ID

##remove column
hist_sample = subset(hist_sample, select = -c(clin_tab.Kids_First_Biospecimen_ID))

pheatmap::pheatmap(
  CC_consensus_mat,
  annotation_col=hist_sample,
  #annotation_colors=anno_palette,
  cluster_rows = results[[3]]$consensusTree,
  cluster_cols = results[[3]]$consensusTree,
  show_rownames = F,
  show_colnames = F,
  filename = file.path(plots_dir, "CC_heatmap.png")
)
