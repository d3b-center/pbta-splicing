suppressPackageStartupMessages({
  library("EnhancedVolcano")
  library("DESeq2")
  library("tidyverse")
  library("dplyr")
  library("vroom")
})


# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data/v5")
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

tpm_count_file <- "ctrl_vs_morpho.rsem.genes.results.tsv"
count_data <- vroom(paste0(data_dir, "/v5/", tpm_count_file)) %>% 
               filter( (CTRL1 + CTRL2 + CTRL3 > 10) & (Treated1 + Treated2 + Treated3 > 10) )

## construct metadata
design = data.frame(row.names = colnames(count_data$gene),
                    condition = c(rep("Ctrl",3), rep("Treated",3) ),
                    libType   = c(rep("paired-end",6)))


## remove first column
count_data_removed <- select(count_data, -gene)

cds = DESeqDataSetFromMatrix(countData=round(count_data_removed),
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
                lab = gsub("ENSG[1234567890]+[.][1234567890]+_", "",count_data$gene), ## remove ensembleid portion , ## remove ensembleid portion
                x = 'log2FoldChange',
                y = 'pvalue',
                #ylim = c(0,21),
                #xlim = c(-3,3),
                title = 'Ctrl vs Treated',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2,
                labSize = 4)

de_results <- cbind(res, count_data) 

write_delim(
  as_tibble(de_results),
  paste0(results_dir, de_output),
  delim = "\t"
)
