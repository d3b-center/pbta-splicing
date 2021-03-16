library("sva")
library("EnhancedVolcano")
library("DESeq2")
library("ggplot2")
library("Hmisc")
library("corrplot")
library("plyr")
library("dplyr")
library("VennDiagram")
library(Hmisc)
library(corrplot)



gene_tpms <- read.delim("/Users/naqvia/Desktop/DIPG/normals_vs_DMG.txt", sep = ",", row.names=1,header=TRUE)
batch <- c(rep(1, 33), rep(2, 22), rep(1,9))

filtered.counts <- gene_tpms[rowSums(gene_tpms>=100) > 10, ]
countTable <- filtered.counts ## filter out  genes with low read counts

## batch correction
adjusted <- ComBat_seq(as.matrix(log2(gene_tpms+1)), batch=batch, group=NULL)
corrected_mat <- 2^(adjusted) 

## convert to dataframe and integer for DESeq2
corrected_df <- as.data.frame(corrected_mat)
corrected_df[] <- sapply( corrected_df, as.integer)

##check to see if batch correction worked?
x=plotMDS(log(corrected_mat), cex.lab= 1, cex = 1, col = c( c(rep("red", 33)),c(rep("darkred",22)),  c(rep("blue",9))), main = "PCA of  DMGs vs Norms Expr", pch=c( rep(16,33), rep(17,22), rep(16,9)))
          
## construct metadata
design = data.frame(row.names = colnames(corrected_df),
                    condition = c(rep("DMG",55), rep("Norm",9)),  
                    libType   = c(rep("paired-end",64)))

singleSamples = design$libType == "paired-end"
new_countTable = (corrected_df[ , singleSamples ])
condition = design$condition[ singleSamples ]

cds = newCountDataSet( (new_countTable), condition )
cds = estimateSizeFactors( cds )
sizeFactors( cds )

cds = estimateDispersions(cds) 

res = nbinomTest( cds, "DMG", "Norm")
res$Significant <- ifelse(res$pval< 0.05, "P-val < 0.05", "Not Sig")

## plot 
EnhancedVolcano(res,
                lab = (res$id),
                x = 'log2FoldChange',
                y = 'pval',
                ylim = c(0,30),
                xlim =c(-10,10),
                title = 'DMG versus Healthy',
                pCutoff = 0.005,
                FCcutoff = 2,
                pointSize = 2.0,
                labSize = 3.0)

write_delim(res, "/Users/naqvia/Desktop/DIPG/res_dmg_healthy_total.tsv", delim= "\t");


##correlation between TAF1 vs its targets
gene_tpm_psi <- read.csv("/Users/naqvia/Desktop/DIPG/NFIC_targets_sign_tpm_vs_psi_tab.tsv", sep = ",", header=TRUE) 

## remove columns that have stdv of 0
tab <- gene_tpm_psi[, sapply(gene_tpm_psi, function(x) { sd(x) != 0} )]

res <- cor(tab)
corrplot(res, type = "upper", order = "hclust", tl.col = "black", tl.cex = 0.5,tl.srt = 100)

res2<-rcorr(as.matrix(tab))

## Plot corr. matrix with insignificant correlations are leaved blank
corrplot(res2$r, type="lower", order="hclust", 
         p.mat = res2$P, sig.level = 0.05, insig = "blank", tl.col = "black", tl.cex = 0.7,tl.srt = 1)

##all genes

gene_tpms <- read.csv("/Users/naqvia/Desktop/DIPG/tpm_ordered_by_polyA_dmg_vs_5norms.SFs_only.csv",row.names=1, header=TRUE)
batch <- c(rep(1, 33), rep(2, 22), rep(1,5))

filtered.counts <- gene_tpms[rowSums(gene_tpms>=10) > 10, ]
countTable <- filtered.counts ## filter out  genes with low read counts

## batch correction
adjusted <- ComBat_seq(as.matrix(log2(gene_tpms+1)), batch=batch, group=NULL)
corrected_mat <- 2^(adjusted) 

## convert to dataframe and integer for DESeq2
corrected_df <- as.data.frame(corrected_mat)
corrected_df[] <- sapply( corrected_df, as.integer )

##check to see if batch correction worked?
x=plotMDS(log(corrected_mat), cex.lab= 1, cex = 1, col = c( c(rep("red", 33)),c(rep("darkred",22)),  c(rep("blue",5))), main = "PCA of  DMGs vs Norms Expr", pch=c( rep(16,33), rep(17,22), rep(16,5)))

## construct metadata
design = data.frame(row.names = colnames(corrected_df),
                    condition = c(rep("DMG",55), rep("Norm",5)),  
                    libType   = c(rep("paired-end",60)))

singleSamples = design$libType == "paired-end"
new_countTable = (corrected_df[ , singleSamples ])
condition = design$condition[ singleSamples ]

cds = newCountDataSet( (new_countTable), condition )
cds = estimateSizeFactors( cds )
sizeFactors( cds )

cds = estimateDispersions(cds) 

res = nbinomTest( cds, "DMG", "Norm")
res$Significant <- ifelse(res$pval< 0.05, "P-val < 0.05", "Not Sig")

## plot 
EnhancedVolcano(res,
                lab = (res$id),
                x = 'log2FoldChange',
                y = 'pval',
                ylim = c(0,20),
                xlim =c(-5,5),
                title = 'DMG versus BS',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 4.0)

