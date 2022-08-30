################################################################################
# volcano_plot_mRNA.R
# written by Ammar Naqvi
#
# usage: Rscript volcano_plot.R
################################################################################


suppressPackageStartupMessages({
  library("sva")
  library("EnhancedVolcano")
  library("DESeq2")
})


# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "mRNA_diff_expr")

input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")


## output files for final plots
file_volc_hgat_SF_plot <- file.path(analysis_dir, "plots", "enhancedVolcano_ctrl_hgat_SFs.png")
file_volc_H3K28_SF_plot <- file.path(analysis_dir, "plots", "enhancedVolcano_ctrl_h3k28_SFs.png")

## get splicing index table
file <- "/tpm_norm_vs_tumor.hgg.SFs.txt"
gene_counts  <-  read.csv(paste0(input_dir, file), row.names = 1,  header=TRUE)

batch <- c(rep(1, 61), rep(2, 137))
filtered.counts <- gene_tpms[rowSums(gene_counts>=2) >= 1, ]

adjusted <- ComBat_seq(as.matrix(log2(filtered.counts+1)), batch=batch, group=NULL)
corrected_mat <- 2^(adjusted)
corrected_df <- as.data.frame(corrected_mat)
corrected_df <- sapply( corrected_df, as.integer )

## construct metadata
design = data.frame(row.names = colnames(corrected_df),
                    condition = c(rep("Healthy",10), rep("Tumor",188) ),
                    libType   = c(rep("paired-end",198)))

singleSamples = design$libType == "paired-end"
new_countTable = (gene_tpms[ , singleSamples ])
condition = design$condition[ singleSamples ]



cds = DESeqDataSetFromMatrix(countData=round(corrected_df),
                             colData=design,
                             design= ~ condition)


cds = estimateSizeFactors( cds )
cds = estimateDispersions(cds)

## run deseq function to compute pvalues
cds <- DESeq(cds)
res <- results(cds)

## label anything below <0.05 as signficant
res$Significant <- ifelse(res$pvalue< 0.05, "P-val < 0.05", "Not Sig")



EnhancedVolcano(res,
                lab = gsub("ENSG[1234567890]+_", "",row.names(filtered.counts)), ## remove ensembleid portion
                x = 'log2FoldChange',
                y = 'pvalue',
                ylim = c(0,33),
                xlim = c(-3,3),
                title = 'Healthy versus HGG',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = ,
                labSize = 3)



ggsave(
  file_volc_hgat_SF_plot,
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width =6.73,
  height = 10.38,
  units = "in",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)


## get count table index table
file <- "/tpm_norm_vs_tumor.h3k28.SFs.v2.txt"
gene_counts  <-  read.csv(paste0(input_dir, file), row.names = 1,  header=TRUE)
batch <- c(rep(1, 45), rep(2, 38))

filtered.counts <- gene_tpms[rowSums(gene_counts>=2) >= 1, ]

adjusted <- ComBat_seq(as.matrix(log2(filtered.counts+1)), batch=batch, group=NULL)
corrected_mat <- 2^(adjusted)
corrected_df <- as.data.frame(corrected_mat)
corrected_df <- sapply( corrected_df, as.integer )


## construct metadata
design = data.frame(row.names = colnames(corrected_df),
                    condition = c(rep("Healthy",10), rep("Tumor",73) ),
                    libType   = c(rep("paired-end",83)))

singleSamples = design$libType == "paired-end"
new_countTable = (gene_counts[ , singleSamples ])
condition = design$condition[ singleSamples ]

cds = DESeqDataSetFromMatrix(countData=round(corrected_df),
                             colData=design,
                             design= ~ condition)

cds = estimateSizeFactors( cds )
sizeFactors( cds )

cds = estimateDispersions(cds)

## run deseq function to compute pvalues
cds <- DESeq(cds)
res <- results(cds)

res$Significant <- ifelse(res$pvalue< 0.05, "P-val < 0.05", "Not Sig")

EnhancedVolcano(res,
                lab = gsub("ENSG[1234567890]+_", "",row.names(filtered.counts)), ## remove ensembleid portion
                x = 'log2FoldChange',
                y = 'pvalue',
                ylim = c(0,27),
                xlim = c(-3,3),
                title = 'Healthy versus H3K28',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2,
                labSize = 3)



ggsave(
  file_volc_H3K28_SF_plot,
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width =6.73,
  height = 10.38,
  units = "in",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)


