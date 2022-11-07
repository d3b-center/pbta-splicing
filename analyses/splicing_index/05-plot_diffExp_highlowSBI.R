################################################################################
# diffExp_highlowSBI.R
#
# generates volcanto plot of differentail expression between high vs low 
# splicing burden tumors using output from generate_splicing_index_tab_using_tumors.pl
#
# written by Ammar Naqvi
#
# usage: Rscript diffExp_highlowSBI.R
################################################################################

## libraries needed
suppressPackageStartupMessages({
  library("dplyr")
  library("EnhancedVolcano")
  library("DESeq2")
  })


# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "splicing_index")

input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

## theme for all plots
# source function for theme for plots survival
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## output files for final plots
file_volc_hgat_plot <- file.path(analysis_dir, "plots", "enhancedVolcano_hgat_sbi.png")
file_volc_non_hgat_plot <- file.path(analysis_dir, "plots", "enhancedVolcano_nonhgat_sbi.png")

#count table for HGAT 
input      = file.path(input_dir,"tab_rsem.str.sbi.hgat.txt")
tab_rsem_hgat <- read.delim(input, header=TRUE, row.names=1)

#count table for non-HGAT  
input      = file.path(input_dir,"tab_rsem.str.sbi.non-hgat.txt")
tab_rsem_non_hgat <- read.delim(input, header=TRUE, row.names=1)

## HGAT differential gene expression analysis
# remove low expression genes
filtered.counts <- tab_rsem_hgat[rowSums(tab_rsem_hgat>=10) >= 38, ]
countTable <- filtered.counts

## construct metadata
design = data.frame(row.names = colnames(countTable),
                    condition = c(rep("High",36), rep("Low",40)))

condition = c(rep("High",36), rep("Low",40) )

cds = DESeqDataSetFromMatrix(countData=round(countTable),
                             colData=design,
                             design= ~ condition)

cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)
cds <- DESeq(cds)

res <- results(cds)

res$Significant <- ifelse(res$pvalue< 0.05, "P-val < 0.05", "Not Sig")

volc_hgat_plot <- EnhancedVolcano(res,
                lab = gsub("ENSG[1234567890]+[.][1234567890]+_", "",row.names(res)), ## remove ensembleid portion
                x = 'log2FoldChange',
                y = 'pvalue',
                ylim = c(0,8),
                xlim = c(-2,2),
                title = 'Low vs High SBI in HGATs',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1,
                labSize = 3,
                typeConnectors = "closed",
                #drawConnectors = TRUE,
                widthConnectors = 0.15,
                colConnectors = 'black')


# Save plot as PNG
png(file_volc_hgat_plot, 
    res = 800, width = 8, height = 8, units = "in")
volc_hgat_plot
dev.off()

##non HGAT differential gene expression analysis
#filter low expressed genes
filtered.counts <- tab_rsem_non_hgat[rowSums(tab_rsem_non_hgat>=10) >= 145, ]
countTable <- filtered.counts

## construct metadata
design = data.frame(row.names = colnames(countTable),
                    condition = c(rep("High",145), rep("Low",145)))

condition = c(rep("High",145), rep("Low",145) )

cds = DESeqDataSetFromMatrix(countData=round(countTable),
                             colData=design,
                             design= ~ condition)

cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)
cds <- DESeq(cds)

res <- results(cds)

res$Significant <- ifelse(res$pvalue< 0.05, "P-val < 0.05", "Not Sig")

volc_non_hgat_plot <- EnhancedVolcano(res,
                lab = gsub("ENSG[1234567890]+[.][1234567890]+_", "",row.names(res)), ## remove ensembleid portion
                x = 'log2FoldChange',
                y = 'pvalue',
                ylim = c(0,30),
                xlim = c(-2,2),
                title = 'Low vs High SBI in non-HGATs',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1,
                labSize = 3,
                typeConnectors = "closed",
                #drawConnectors = TRUE,
                widthConnectors = 0.15,
                colConnectors = 'black')

# Save plot as PNG
png(file_volc_non_hgat_plot, 
    res = 800, width = 8, height = 8, units = "in")
volc_non_hgat_plot
dev.off()
