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

suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("rlist"))


# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## do all of the samples // pan_cancer_splicing.thr10.report_all.txt
dataDir = "~/Desktop/AS-DMG/analyses/pan_cancer/results/"

##results/pan_cancer_splicing.thr10.report_select.txt
file_psi <- "pan_cancer_splicing.thr10.report_select.txt"

##removing duplicates, cell lines, and second malignancies
file_psi <- "pan_cancer_splicing.thr10.report_select.remDup.txt"
psi_tab  = read.delim(paste0(dataDir, file_psi), sep = "\t", row.names=1, header=TRUE)

rnames <- psi_tab[,1]
row.names(psi_tab) <- psi_tab$Splice_ID
mat_hm <- data.matrix(psi_tab[,2:ncol(psi_tab)])

d=mat_hm
d<-data.matrix(psi_tab)


#psi_tab$Splice_ID
means <- rowMeans(psi_tab)
vars <- apply(psi_tab, 1, var)
varorder <- order(vars, decreasing = T)
psi_tab <- psi_tab[varorder, ] %>% head(10000)

# tumorData mean and sd
tumorData_means <- rowMeans(psi_tab, na.rm = TRUE)
tumorData_sd <- apply(psi_tab, 1, sd, na.rm = TRUE)

# subtract mean
GeneExpzscored <- sweep(psi_tab, 1, tumorData_means, FUN = "-")

# divide by SD remove NAs and Inf values from zscore for genes with 0
GeneExpzscored <- sweep(GeneExpzscored, 1, tumorData_sd, FUN = "/") %>%
  na_if(Inf) %>%
  na.omit()

## reduce the dataset to the top 10,000 most variable genes, measured by median absolute deviation
mads=apply(d,1,mad)
d=d[rev(order(mads))[1:10000],]

## the default settings of the agglomerative hierarchical clustering algorithm using Pearson correlation distance, so it is appropriate to gene median center d using
d = sweep(d,1, apply(d,1,median,na.rm=T))

## remove NAs
is.na(d) <- sapply(d, is.infinite)

d[is.na(d)] <- 0
d[is.nan(d)] <- 0

## k= 6 clusters
results = ConsensusClusterPlus((d),maxK=10,reps=100,pItem=0.8,
                     title="clustering",clusterAlg="hc",distance="spearman",seed=123,innerLinkage = "average", finalLinkage = "average")
#GeneExpzscored
##save to file
#title="/Users/naqvia/Desktop/pan_cancer_rmats"
#results= ConsensusClusterPlus((d),maxK=10,reps=100,pItem=0.8,
#                     title="Total Samples PSI Clustering",clusterAlg="hc",distance="spearman",seed=123,innerLinkage = "average", finalLinkage = "average", plot="png")



# choose a cluster that seems best and assign to n_cluster
CC_group <- results[[6]]$consensusClass %>%
  as.data.frame()
colnames(CC_group) <- "CC"

write.table(CC_group, file = "/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/results/CC_groups_remDup.txt", row.names=TRUE, sep="\t")

# read in consensus clustering matrix
CC_consensus_mat <- results[[6]]$consensusMatrix
colnames(CC_consensus_mat) <- rownames(CC_group)
rownames(CC_consensus_mat) <- rownames(CC_group)



clin_file = "~/Desktop/AS-DMG/data/pbta-histologies.RNA-Seq.addedCluster.tsv"
clin_tab = read.delim(clin_file, sep = "\t", header=TRUE)

clin_file = "/Users/naqvia/Desktop/AS-DMG/data/compiled_molecular_subtypes_with_clinical_pathology_feedback.tsv"
clin_tab = read.delim(clin_file, sep = "\t", header=TRUE)

#cluster_mem = "results/CC_select_members.txt"
cluster_file <- "CC_select_members.txt"
cluster_mem  = read.delim(paste0(dataDir, cluster_file), sep = "\t", header=TRUE)

hist_sample <- cbind(data.frame(clin_tab$Kids_First_Biospecimen_ID), data.frame(clin_tab$short_histology), data.frame(clin_tab$molecular_subtype), data.frame(clin_tab$Cluster))


#library(randomcoloR)

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
  cluster_rows = results[[6]]$consensusTree,
  cluster_cols = results[[6]]$consensusTree, 
  show_rownames = F,
  show_colnames = F
  #filename = "cc_test/png",
  #fontsize = 10, cellwidth = 10, cellheight = 10
)


##pie chart of clusters
cluster_members <- "cc_k7_memberships.txt"
cluster_members  = read.delim(paste0(dataDir, cluster_members), sep = "\t", header=TRUE)


ggplot(cluster_members, aes("", Counts, fill = Hist)) +
  geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
  coord_polar("y") +
  geom_text(aes(label = paste0(round(Counts), "%")), 
            position = position_stack(vjust = 0.5)) +
  labs(x = NULL, y = NULL, fill = NULL, 
       title = "market share") +
  guides(fill = guide_legend(reverse = TRUE)) +
  #scale_fill_manual(values = c("#ffd700", "#bcbcbc", "#ffa500", "#254290")) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "#666666"))

# Factominer PCA
# merge multi-genomics data as dataframe of most variant feature sets
res.pca <- FactoMineR::PCA(t(psi_tab), graph = T)
res.hcpc <- FactoMineR::HCPC(res.pca, nb.clust = -1, graph = T, metric = "euclidean", method = "average", graph.scale = "intertia")

##stacked barplot of cluster memberships

cluster_mem_file = "/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/results/cluster_hist_makeup.txt"
cluster_mem      = read.delim(cluster_mem_file, sep = "\t",header=TRUE)

# Stacked
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(cluster_mem, aes(fill=Hist, y=Counts, x=Cluster)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_brewer(palette="Accent") + theme_Publication()

##stacked barplot of cluster memberships

cluster_mem_hgat_file = "/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/results/cluster_hist_makeup.hgats.v2.txt"
cluster_mem_hgat      = read.delim(cluster_mem_hgat_file, sep = "\t",header=TRUE)

# Stacked
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(cluster_mem_hgat, aes(fill=Subtype, y=Counts, x=Cluster)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_brewer(palette="Spectral") + theme_Publication()


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

