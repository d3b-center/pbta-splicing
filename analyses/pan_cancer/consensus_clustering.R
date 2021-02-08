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
ConsensusClusterPlus((d),maxK=7,reps=100,pItem=0.8,
                     title="clustering",clusterAlg="hc",distance="spearman",seed=123,innerLinkage = "average", finalLinkage = "average")

##save to file
#title="/Users/naqvia/Desktop/pan_cancer_rmats"
#results= ConsensusClusterPlus((d),maxK=7,reps=100,pItem=0.8,
#                     title="Total Samples PSI Clustering",clusterAlg="hc",distance="spearman",seed=123,innerLinkage = "average", finalLinkage = "average", plot="png")
