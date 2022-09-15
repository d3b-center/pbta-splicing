################################################################################
# highExon4_vs_lowExon4_diffExpr_volcano
# written by Ammar Naqvi
#
# usage: Rscript highExon4_vs_lowExon4_diffExpr_volcano.R
################################################################################

suppressPackageStartupMessages({
  library("ggplot2")
  library("dplyr")
  library("tidyverse")
  library("EnhancedVolcano")
  library("DESeq2")
})

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1_specific")

input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

## get data files and make table
input      = file.path(data_dir,"stranded_gene_counts_tab.tsv")
tab_rsem <- read.delim(input, header=TRUE, row.names=1)


## extract relevant samples only (top 4 and lowest 4)
tab_rsem = tab_rsem %>%
  select(c('gene_id', 'BS_Q13FQ8FV', 'BS_ZV1P6W9C','BS_WH8G4VFB','BS_NNPEC7W1','BS_PZVHMSYN', 
           'BS_DRY58DTF','BS_GXTFW99H','BS_E60JZ9Z3','BS_9CA93S6D'))

countTable <- tab_rsem[rowSums(tab_rsem>=10) > 9, ]
rownames(countTable) <- countTable$gene_id
countTable <- (countTable[-c(1)])

## construct metadata
design = data.frame(row.names = colnames(countTable),
                    condition = c(rep("High",5), rep("Low",4)))
condition = c(rep("High",5), rep("Low",4) )

cds = DESeqDataSetFromMatrix(countData=round(countTable),
                             colData=design,
                             design= ~ condition)

cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)
cds <- DESeq(cds)

res <- results(cds)

res$Significant <- ifelse(res$pvalue< 0.05, "P-val < 0.05", "Not Sig")

EnhancedVolcano(res,
                lab = gsub("ENSG[1234567890]+[.][1234567890]+_", "",row.names(res)), ## remove ensembleid portion
                x = 'log2FoldChange',
                y = 'pvalue',
                ylim = c(0,8),
                xlim = c(-4,4),
                title = 'High Exon 4  versus Low Exon 4 Tumors',
                pCutoff = 0.005,
                FCcutoff = 2,
                pointSize = 1,
                labSize = 3,
                typeConnectors = "closed",
                #drawConnectors = TRUE,
                widthConnectors = 0.15,
                colConnectors = 'black')

plot_file = file.path(plots_dir,"highExon4_vs_lowExon4.volano.pdf") 
ggsave(plot_file, width = 20, height = 15)


__DATA__
BS_Q13FQ8FV
BS_ZV1P6W9C
BS_WH8G4VFB
BS_NNPEC7W1
BS_PZVHMSYN

CLK1 low exon inclusion (group 2)
BS_XM1AHBDJ
BS_DRY58DTF
BS_GXTFW99H
BS_E60JZ9Z3
BS_9CA93S6D