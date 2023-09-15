################################################################################
# 01-diffExpr-ctrl_vs_morph.R
# Gene Set Enrichment Analysis using ClusterProfiler
# Author: Ammar Naqvi and Shehbeel Arif
# Script adapted from ALSF's DE Analysis Tutorial (https://alexslemonade.github.io/refinebio-examples/03-rnaseq/pathway-analysis_rnaseq_02_gsea.html)
################################################################################

suppressPackageStartupMessages({
  library("EnhancedVolcano")
  library("DESeq2")
  library("tidyverse")
  library("dplyr")
  library("vroom")
  library("clusterProfiler")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data/")
analysis_dir <- file.path(root_dir, "analyses", "CLK1_splicing_impact")

input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")

if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive=TRUE)
}

## ouput files
de_output = "ctrl_vs_treated.de.tsv"
file_volc_plot = "ctrl_vs_clk1-morp_volcano.tiff"

tpm_count_file <- "ctrl-vs-morpholino-gene-counts-rsem-expected_count.tsv"
count_data <- vroom(paste0(data_dir, tpm_count_file)) %>% 
               filter( (CTRL1 + CTRL2 + CTRL3 > 10) & (Treated1 + Treated2 + Treated3 > 10) )

## construct metadata
design = data.frame(row.names = colnames(count_data)[-1],
                    condition = c(rep("Ctrl",3), rep("Treated",3) ),
                    libType   = c(rep("paired-end",6)))


## remove first column
count_data_removed <- dplyr::select(count_data, -gene)

cds = DESeqDataSetFromMatrix(countData=round(count_data_removed),
                             colData=design,
                             design= ~ condition)


## run DeSeq function to compute pvalues
cds <- DESeq(cds)
res <- results(cds)

## label anything below <0.05 as signficant
res$Significant <- ifelse(res$pvalue< 0.05, "P-val < 0.05", "Not Sig")

EnhancedVolcano(res,
                lab = gsub("ENSG[1234567890]+[.][1234567890]+_", "",count_data$gene), ## remove ensembleid portion , ## remove ensembleid portion
                x = 'log2FoldChange',
                y = 'pvalue',
                #ylim = c(0,21),
                #xlim = c(-3,3),
                title = 'Treated vs Ctrl',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2,
                labSize = 4)
ggsave(
  paste0(plots_dir,"/", file_volc_plot),
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

de_results <- cbind(res, count_data) 

write_delim(
  as_tibble(de_results),
  paste0(results_dir,"/", de_output),
  delim = "\t"
)

## brain-goi list
brain_goi_file <- "brain-goi-list-new.txt"
brain_goi_data <- read.table(file.path(input_dir, brain_goi_file),header=TRUE)

## res data with filtered goi from above 
res2 <- as.data.frame(res) %>% dplyr::mutate(gene=gsub("ENSG[1234567890]+[.][1234567890]+_", "",count_data$gene)) %>% inner_join(brain_goi_data, by="gene")

EnhancedVolcano(res2,
                lab = res2$gene, ## remove ensembleid portion , ## remove ensembleid portion
                x = 'log2FoldChange',
                y = 'pvalue',
                #ylim = c(0,21),
                #xlim = c(-3,3),
                title = 'Ctrl vs Treated',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 2,
                labSize = 4)

ggsave(
  paste0(plots_dir,"/", "ctrl_vs_clk1-morp_volcano.brain-goi.tiff"),
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

unlink(file.path("Rplots.pdf"))