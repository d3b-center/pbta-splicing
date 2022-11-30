################################################################################
# 03-plot_highExon4_vs_lowExon4_diffExpr_volcano.R
# written by Ammar Naqvi
#
# usage: Rscript 03-plot_highExon4_vs_lowExon4_diffExpr_volcano.R
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
file_gene_counts <- file.path(data_dir, "gene-counts-rsem-expected_count-collapsed.rds") 

## extract relevant samples only (top 4 exon 4 inclusion vs and lowest ex4 inclusion)
tab_rsem <- readRDS(file_gene_counts) %>%  
  select(c('BS_Q13FQ8FV', 'BS_ZV1P6W9C','BS_WH8G4VFB','BS_NNPEC7W1',
           'BS_PZVHMSYN',  'BS_DRY58DTF','BS_GXTFW99H','BS_E60JZ9Z3',
           'BS_9CA93S6D'))

# remove genes which do not have TPM >=10 in at least 8 samples
countTable <- tab_rsem[rowSums(tab_rsem>=10) > 8, ]

## construct metadata
design = data.frame(row.names = colnames(countTable),
                    condition = c(rep("High",5), rep("Low",4)),
                    libType   = c(rep("paired-end",9)))

cds = DESeqDataSetFromMatrix(countData=round(countTable),
                             colData=design,
                             design= ~ condition)

cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)
cds <- DESeq(cds)

res <- results(cds)

# add significance column 
# res$Significant <- ifelse(res$padj< 0.05, "P-val < 0.05", "Not Sig")

volc_plot <- EnhancedVolcano(res,
                lab = gsub("ENSG[1234567890]+[.][1234567890]+_", "",row.names(res)), ## remove ensembleid portion
                x = 'log2FoldChange',
                y = 'padj',
                ylim = c(0,8),
                xlim = c(-4,4),
                title = 'Tumors with low vs. high CLK1 exon 4 inclusion',
                pCutoff = 0.05,
                subtitle = NULL,
                FCcutoff = 2,
                pointSize = 1,
                labSize = 3,
                typeConnectors = "closed",
                widthConnectors = 0.15,
                colConnectors = 'black',
                caption = NULL,
                ylab = bquote(bold(~-log[10] ~ p[adj])),
                xlab = bquote(bold(~log[2] ~ "Fold Change")),
                legendLabels = c("NS", expression(log[2] ~ FC), "adj p-value", expression("adj p-value and" ~ log[2] ~ FC)))

plot_file = file.path(plots_dir,"high_v_low_clk1_exon4_inclusion.volano.pdf") 

# Save plot as PDF
pdf(plot_file, 
    width = 10, height = 10)
volc_plot
dev.off()







