################################################################################
# consensus_clustering.R
# script that takes in PSI data file and performs consensus clustering
# written by Ammar Naqvi
#
# usage: Rscript consensus_clustering.R
################################################################################

## do all of the samples // pan_cancer_splicing.thr10.report_all.txt
dataDir = "~/Desktop/AS-DMG/data/"
file_fus <- "pan_cancer_splicing.thr10.report_all.txt"
psi_tab  = read.delim(paste0(dataDir, file_fus), sep = "\t", header=TRUE)

rnames <- psi_tab[,1]
row.names(psi_tab) <- psi_tab$Splice_ID
mat_hm <- data.matrix(psi_tab[,2:ncol(psi_tab)])

d=mat_hm

## reduce the dataset to the top 1,000 most variable genes, measured by median absolute deviation
mads=apply(d,1,mad)
d=d[rev(order(mads))[1:1000],]

## the default settings of the agglomerative hierarchical clustering algorithm using Pearson correlation distance, so it is appropriate to gene median center d using
d = sweep(d,1, apply(d,1,median,na.rm=T))


## remove NAs
is.na(d) <- sapply(d, is.infinite)

d[is.na(d)] <- 0
d[is.nan(d)] <- 0

## k= 7 clusters
results = ConsensusClusterPlus((d),maxK=7,reps=100,pItem=0.8,
                     title="clustering",clusterAlg="hc",distance="spearman",seed=123,innerLinkage = "average", finalLinkage = "average")

##save to file
#title="/Users/naqvia/Desktop/pan_cancer_rmats"
#results= ConsensusClusterPlus((d),maxK=7,reps=100,pItem=0.8,
#                     title="Total Samples PSI Clustering",clusterAlg="hc",distance="spearman",seed=123,innerLinkage = "average", finalLinkage = "average", plot="png")



# choose a cluster that seems best and assign to n_cluster
CC_group <- results[[7]]$consensusClass %>%
  as.data.frame()
colnames(CC_group) <- "CC"

# read in consensus clustering matrix
CC_consensus_mat <- results[[7]]$consensusMatrix
colnames(CC_consensus_mat) <- rownames(CC_group)
rownames(CC_consensus_mat) <- rownames(CC_group)


clin_file = "~/Desktop/AS-DMG/data/pbta-histologies.RNA-Seq.tsv"
clin_tab = read.delim(clin_file, sep = "\t", header=TRUE)

hist_sample <- cbind(data.frame(clin_tab$Kids_First_Biospecimen_ID), data.frame(clin_tab$short_histology))
library(randomcoloR)

## n is total num of histology (short)
anno_palette <- distinctColorPalette(43)

names(anno_palette)<- unique(clin_tab$short_histology)

hist_sample$clin_tab.short_histology <- as.factor(hist_sample$clin_tab.short_histology )
names(anno_palette)<- levels(hist_sample$clin_tab.short_histology)
rownames(hist_sample)<- hist_sample$clin_tab.Kids_First_Biospecimen_ID

##remove colum
hist_sample = subset(hist_sample, select = -c(clin_tab.Kids_First_Biospecimen_ID))

pheatmap::pheatmap(
  CC_consensus_mat,
  #color = colorRampPalette(brewer.pal(n = 9, name = "Blues"))(9),
  annotation_col=hist_sample, 
  #annotation_colors=anno_palette,
  cluster_rows = results[[7]]$consensusTree,
  cluster_cols = results[[7]]$consensusTree, 
  show_rownames = F,
  show_colnames = F
  #filename = "cc_test/png",
  #fontsize = 10, cellwidth = 10, cellheight = 10
)

# Generate annotations for rows and columns
annotation_col = data.frame(
  CellType = factor(rep(c("CT1", "CT2"), 5)), 
  Time = 1:5
)

