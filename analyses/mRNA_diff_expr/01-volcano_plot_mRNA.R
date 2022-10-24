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
analysis_dir <- file.path(root_dir, "analyses", "mRNA_diff_expr")

input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

## output files for final plots
file_volc_hgat_SF_plot <- file.path(analysis_dir, "plots", 
                                    "enhancedVolcano_ctrl_hgat_SFs.png")
file_volc_H3K28_SF_plot <- file.path(analysis_dir, "plots", 
                                     "enhancedVolcano_ctrl_h3k28_SFs.png")

## get splicing factor list to subset later
sf_file = "splicing_factors.txt"
sf_list <- read.csv(paste0(input_dir, "/", sf_file),  header=FALSE)

## get cliincal histlogy file fitlered by HGG samples
clin_file = "v1/histologies.tsv"
clin_tab <- read.delim(paste0(data_dir,"/",clin_file), sep = "\t", header=TRUE) %>% 
  filter(short_histology == 'HGAT') %>% filter(RNA_library == 'stranded') %>%
                                          filter(cohort == 'PBTA') %>%
                                          filter(CNS_region == 'Midline')

## get gene count table with splicing factors and midline HGGs filter
file_gene_counts = "gene-counts-rsem-expected_count-collapsed.rds" 
count_data <- readRDS(paste0(data_dir, "/v1/", file_gene_counts)) 

#filter for only splicing factors from above SF list
count_data_sf <- count_data[rownames(count_data) %in% sf_list$V1, ]  %>% 
                select(any_of(clin_tab$Kids_First_Biospecimen_ID))       

#filter for HGG midline samples
count_data_sf <- cbind(gene = rownames(count_data_sf), count_data_sf)
rownames(count_data_sf) <- NULL

## get corresponding non-brain tumor samples
file_nontumor_count = "rsem_counts.non_tumor.tsv"
gene_counts_nontumor  <-  read.delim(paste0(input_dir, "/",file_nontumor_count),header=TRUE, sep = "\t")

gene_counts_combined <- inner_join(gene_counts_nontumor,count_data_sf, by = 'gene')
filtered.counts <- gene_counts_combined[rowSums(gene_counts_combined>=2) >= 1, ]

## construct metadata
design = data.frame(row.names = colnames(filtered.counts$gene),
                    condition = c(rep("Healthy",10), rep("Tumor",53) ),
                    libType   = c(rep("paired-end",63)))

singleSamples = design$libType == "paired-end"
new_countTable = (filtered.counts[ , singleSamples ])
condition = design$condition[ singleSamples ]

## remove first column
filtered.counts_removed <- select(filtered.counts, -gene)

cds = DESeqDataSetFromMatrix(countData=round(filtered.counts_removed),
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
                lab = filtered.counts$gene, ## remove ensembleid portion
                x = 'log2FoldChange',
                y = 'pvalue',
                ylim = c(0,21),
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


## compute gene expression SF between Midline HGGs vs H3K28M
clin_file = "v1/histologies.tsv"

## get cliincal histlogy file fitlered by HGG samples
clin_tab_hgg <- read.delim(paste0(data_dir,"/",clin_file), sep = "\t", 
                           header=TRUE) %>% filter(short_histology == 'HGAT') %>% 
                                            filter(RNA_library == 'stranded') %>%
                                            filter(cohort == 'PBTA') %>%
                                            filter(CNS_region == 'Midline')  %>%
                                            filter(pathology_diagnosis == 'High-grade glioma/astrocytoma (WHO grade III/IV)')
## get DIPG samples
clin_tab_dipg <- read.delim(paste0(data_dir,"/",clin_file), sep = "\t", 
                            header=TRUE) %>% filter(short_histology == 'HGAT') %>% 
                                             filter(RNA_library == 'stranded') %>%
                                             filter(cohort == 'PBTA') %>%
                                             filter(CNS_region == 'Midline')  %>%
                                             filter(pathology_diagnosis == 'Brainstem glioma- Diffuse intrinsic pontine glioma')

file_gene_counts = "gene-counts-rsem-expected_count-collapsed.rds" 
count_data <- readRDS(paste0(data_dir, "/v1/", file_gene_counts)) 

## filter for only splicing factors from above SF list
count_data_sf <-count_data[rownames(count_data) %in% sf_list$V1, ] 

#filter for HGG and DIPG midline samples
count_data_sf_hgg <- select(count_data_sf, any_of(clin_tab_hgg$Kids_First_Biospecimen_ID) )
count_data_sf_dipg <- select(count_data_sf,any_of(clin_tab_dipg$Kids_First_Biospecimen_ID) )

count_data_sf_hgg <- cbind(gene = rownames(count_data_sf_hgg), count_data_sf_hgg)
rownames(count_data_sf_hgg) <- NULL

count_data_sf_dipg <- cbind(gene = rownames(count_data_sf_dipg), count_data_sf_dipg)
rownames(count_data_sf_dipg) <- NULL

## combine in count table and filter
gene_counts_combined <- inner_join(count_data_sf_hgg,count_data_sf_dipg, by = 'gene')
filtered.counts <- gene_counts_combined[rowSums(gene_counts_combined>=2) >= 10, ]

## construct metadata
design = data.frame(row.names = colnames(filtered.counts$gene),
                    condition = c(rep("Non-DIPG",40), rep("DIPG",13) ),
                    libType   = c(rep("paired-end",53)))

singleSamples = design$libType == "paired-end"
new_countTable = (filtered.counts[ , singleSamples ])
condition = design$condition[ singleSamples ]

## remove first column
filtered.counts_removed <- select(filtered.counts, -gene)

cds = DESeqDataSetFromMatrix(countData=round(filtered.counts_removed),
                             colData=design,
                             design= ~ condition)

cds = estimateSizeFactors( cds )
cds = estimateDispersions(cds)

## run deseq function to compute pvalues
cds <- DESeq(cds)
res <- results(cds)

## label anything below <0.05 as significant
res$Significant <- ifelse(res$pvalue< 0.05, "P-val < 0.05", "Not Sig")

EnhancedVolcano(res,
                lab = filtered.counts$gene, ## remove ensembleid portion
                x = 'log2FoldChange',
                y = 'pvalue',
                ylim = c(0,5),
                xlim = c(-3,3),
                title = 'Non-DIPG versus DIPG',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = ,
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


