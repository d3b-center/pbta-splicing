################################################################################
# 01-volcano_plot_mRNA.R
# written by Ammar Naqvi
#
# This script uses DESeq2 to perform differential gene expression on select 
# samples annotaed by the histology files. THe results are plotted in a volcano 
# plot and splicing factor genes are then selected based on that for downstream 
# analyses. 
#
# usage: Rscript volcano_plot_mRNA.R
################################################################################


suppressPackageStartupMessages({
  library("EnhancedVolcano")
  library("DESeq2")
  library("tidyverse")
  library("dplyr")
})


# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "splicing-factor_dysregulation")

input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")

if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive=TRUE)
}

## output files for final plots
file_volc_hgg_SF_plot <- file.path(analysis_dir, "plots", 
                                    "enhancedVolcano_ctrl_hggs_SFs.pdf")

## input files
sf_file <- file.path(input_dir,"splicing_factors.txt")
clin_file <- file.path(data_dir,"histologies.tsv")
file_gene_counts = file.path(data_dir,"gene-counts-rsem-expected_count-collapsed.rds")
file_nontumor_count = file.path(input_dir,"rsem_counts.non_tumor.tsv")


## get splicing factor list to subset later
sf_list <- read.csv(sf_file, header=FALSE)

## get clinical histlogy file filtered by HGG samples
clin_tab <- read.delim(clin_file, sep = "\t", header=TRUE) %>% 
  filter(short_histology == 'HGAT') %>% filter(RNA_library == 'stranded') %>%
                                          filter(cohort == 'PBTA') %>%
                                          filter(CNS_region == 'Midline')

## get gene count table with splicing factors and midline HGGs filter
count_data <- readRDS(file_gene_counts)

#filter for only splicing factors from above SF list
count_data_sf <- count_data[rownames(count_data) %in% sf_list$V1, ]  %>% 
                select(any_of(clin_tab$Kids_First_Biospecimen_ID))       

#filter for HGG midline samples
count_data_sf <- cbind(gene = rownames(count_data_sf), count_data_sf)
rownames(count_data_sf) <- NULL

## get corresponding non-brain tumor samples
gene_counts_nontumor  <-  read.delim(file_nontumor_count,header=TRUE, sep = "\t")

gene_counts_combined <- inner_join(gene_counts_nontumor,count_data_sf, by = 'gene')
filtered.counts <- gene_counts_combined[rowSums(gene_counts_combined>=2) >= 1, ]

## construct metadata
design = data.frame(row.names = colnames(filtered.counts$gene),
                    condition = c(rep("Healthy",10), rep("Tumor",53) ),
                    libType   = c(rep("paired-end",63)))

## remove first column
filtered.counts_removed <- select(filtered.counts, -gene)

cds = DESeqDataSetFromMatrix(countData=round(filtered.counts_removed),
                             colData=design,
                             design= ~ condition)

## run deseq function to compute pvalues
cds <- DESeq(cds)
res <- results(cds)

## label anything below <0.05 as signficant
res$Significant <- ifelse(res$pvalue< 0.05, "P-val < 0.05", "Not Sig")

EnhancedVolcano(res,
                lab = filtered.counts$gene, ## remove ensembleid portion
                x = 'log2FoldChange',
                y = 'pvalue',
                ylim = c(0,21),
                xlim = c(-3,3),
                title = 'non-Tumor versus HGG',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3,
                labSize = 5)


ggsave(
  file_volc_hgg_SF_plot,
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width =6.73,
  height = 10.38,
  units = "in",
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)

## write significant genes to table for subsequent correlation analyses
gene_sign_list <- as.data.frame(res) %>% mutate(gene = filtered.counts$gene) %>% filter(padj < 0.05) %>% filter(abs(log2FoldChange) > 1)  %>% select(gene) 
write_delim(gene_sign_list,file.path(results_dir,"sign_genes.txt"), delim = "\t")

