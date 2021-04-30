library("EnhancedVolcano")
library("DESeq2")

full_path = "/Users/naqvia/Desktop/DIPG/gene_counts_tpm.sfactors.rounded_ceil.v2.tsv"
countTable <- read.table(full_path,header=TRUE,row.names=1)
head(countTable)

##construct metadata
design = data.frame(row.names = colnames( countTable ),
                    condition = c(rep("norm",5), rep("DIPG",48) ),
                    libType   = c(rep("paired-end",53)))

singleSamples = design$libType == "paired-end"
new_countTable = countTable[ , singleSamples ]
condition = design$condition[ singleSamples ]

cds = newCountDataSet( new_countTable, condition )
cds = estimateSizeFactors( cds )
sizeFactors( cds )
cds = estimateDispersions( cds )

res = nbinomTest( cds, "norm", "DIPG")
res$Significant <- ifelse(res$pval< 0.05, "P-val < 0.05", "Not Sig")

ggplot(res, aes(x = log2FoldChange, y = -log10(pval))) +
  geom_point(aes(color = Significant)) +
  scale_color_manual(values = c("grey", "red")) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom") +
  geom_text_repel(
    data = subset(res, pval <= 0.05 & abs(res$log2FoldChange)>=1),
    aes(label = id),
    size = 2,
    box.padding = unit(.5, "lines"),
    point.padding = unit(.5, "lines"))

##enhanced
EnhancedVolcano(res,
                lab = (res$id),
                x = 'log2FoldChange',
                y = 'pval',
                xlim = c(-5,5),
                ylim = c(0,14),
                title = 'BS vs DMG (SFs)',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.0,
                labSize = 3.0)


## with batch tpm_ordered_by_polyA_dmg_vs_5norms.SFs_only.csv
library("EnhancedVolcano")
library("DESeq2")
library("sva")




## volcano of HGGs vs DMGs

## expression of HGGs + DMGs
gene_tpms <- read.csv("/Users/naqvia/Desktop/DIPG/gene_counts_tpm.by_status_HGGs.txt",row.names=1, header=TRUE)
gene_tpms <- read.csv("/Users/naqvia/Desktop/DIPG/gene_counts_tpm.by_status_HGGs.SFs.txt",row.names=1, header=TRUE)
countTable <- gene_tpms
head(countTable)

filtered.counts <- countTable[rowSums(countTable>=10) > 10, ]
countTable <- filtered.counts

# x=plotMDS(log(gene_tpms), cex.lab= 1, cex = 1, col = c( c(rep("darkgreen",11)),c(rep("red",36)), c(rep("blue",137)) ), main = "PCA of HGG Expr", pch=c(rep(19,184)))

##construct metadata
design = data.frame(row.names = colnames( countTable ),
                    condition = c(rep("DMG",47), rep("non-DMG",135) ),
                    libType   = c(rep("paired-end",182)))

singleSamples = design$libType == "paired-end"
new_countTable = countTable[ , singleSamples ]
condition = design$condition[ singleSamples ]

cds = newCountDataSet( new_countTable, condition )
cds = estimateSizeFactors( cds )
sizeFactors( cds )
cds = estimateDispersions( cds )

res = nbinomTest( cds, "DMG", "non-DMG")
res$Significant <- ifelse(res$pval< 0.05, "P-val < 0.05", "Not Sig")


#x=plotMDS(log(gene_tpms), cex.lab= 1, cex = 1, col = c( c(rep("darkgreen",11)),c(rep("red",36)), c(rep("blue",137)) ), main = "PCA of HGG Expr", pch=c(rep(19,184)))

##enhanced
EnhancedVolcano(res,
                lab = (res$id),
                x = 'log2FoldChange',
                y = 'pval',
                xlim = c(-6,6),
                ylim = c(0,20),
                title = 'DMG vs non-DMG',
                pCutoff = 0.005,
                FCcutoff = 2,
                pointSize = 1.0,
                labSize = 3.0)

## with batch correction ##
library("DESeq2")
library("EnhancedVolcano")

library("sva")
library("EnhancedVolcano")
library("DESeq2")

gene_tpms <- read.csv("/Users/naqvia/Desktop/DIPG/tpm_ordered_by_polyA_dmg_only.csv",row.names=1, header=TRUE)
batch <- c(rep(1, 33), rep(2, 22))
filtered.counts <- countTable[rowSums(corrected_mat>=10) > 10, ]
countTable <- filtered.counts


adjusted <- ComBat_seq(as.matrix(log2(gene_tpms+1)), batch=batch, group=NULL)
corrected_mat <- 2^(adjusted)

corrected_df <- as.data.frame(corrected_mat)
corrected_df <- sapply( corrected_df, as.integer )

## construct metadata
design = data.frame(row.names = colnames(corrected_df),
                    condition = c(rep("H3-WT",5), rep("H3K28",28), rep("H3-WT", 6), rep("H3K28", 16) ),
                    libType   = c(rep("paired-end",55)))

singleSamples = design$libType == "paired-end"
new_countTable = (corrected_df[ , singleSamples ])
condition = design$condition[ singleSamples ]

cds = newCountDataSet( (new_countTable), condition )
cds = estimateSizeFactors( cds )
sizeFactors( cds )

cds = estimateDispersions(cds)

res = nbinomTest( cds, "H3-WT", "H3K28")
res$Significant <- ifelse(res$pval< 0.05, "P-val < 0.05", "Not Sig")


#x=plotMDS(log(corrected_mat), cex.lab= 1, cex = 1, col = c( c(rep("blue", 19)),c(rep("red",5)), c(rep("darkred",28 )),c(rep("blue",101)), c(rep("red",6)), c(rep("darkred", 16))), main = "PCA of HGGs vs DMGs", pch=c(rep(16,52), c(rep(17,123))))
#legend("bottomright", legend=c("HGGs-H3 WT PolyA","DMG-H3 WT PolyA","DMG-H3 K28 PolyA","HGGs-H3 WT Non-PolyA","DMG-H3 WT Non-PolyA","DMG-H3 K28 Non-PolyA" ),col = c("blue","red","darkred","blue","red","darkred"), pch = c(16,16,16,17,17,17), horiz=TRUE, cex=0.5)


EnhancedVolcano(res,
                lab = (res$id),
                x = 'log2FoldChange',
                y = 'pval',
                ylim = c(0,2.5),
                title = 'H3-WT versus H3K27',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.0,
                labSize = 2.0)





##theme
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
