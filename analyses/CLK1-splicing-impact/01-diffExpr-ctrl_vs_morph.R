################################################################################
# 01-diffExpr-ctrl_vs_morph.R
# Differential gene expression analysis and plots
# Author: Ammar Naqvi 
# usage: Rscript 01-diffExpr-ctrl_vs_morph.R
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
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing-impact")

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
de_reg_output = "ctrl_vs_treated.de.regl.tsv"
file_volc_plot = "ctrl_vs_clk1-morp_volcano.pdf"


## input files
# count data
tpm_count_file <- "ctrl-vs-morpholino-gene-counts-rsem-expected_count.tsv" 

# gene lists
known_rbp_file <- file.path(input_dir,'RBP_known.txt')
known_kinase_file <- file.path(input_dir,'kinase_known.txt')
known_tf_file <- file.path(input_dir,'tf_known.txt')
known_kinase_file <- file.path(input_dir,'kinase_known.txt')
known_epi_file <- file.path(input_dir,'epi_known.txt')
known_ts_file <- file.path(input_dir,'ts_known.txt')
known_onco_file <- file.path(input_dir,'onco_known.txt')
brain_goi_file <- file.path(input_dir,'brain-goi-list-new.txt')

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
                drawConnectors = TRUE,
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

## subset with gene lists
brain_goi_data <- read.table(brain_goi_file, header=TRUE)
known_rbp <- read.table(known_rbp_file,header=FALSE) %>% 
  rename('gene' = V1)
known_epi <- read.table(known_epi_file,header=FALSE) %>% 
  rename('gene' = V1)
known_kinase <- read.table(known_kinase_file,header=FALSE) %>% 
  rename('gene' = V1)
known_tf <- read.table(known_tf_file,header=FALSE) %>% 
  rename('gene' = V1)
known_onco <- read.table(known_onco_file,header=FALSE) %>% 
  rename('gene' = V1)
known_ts <- read.table(known_ts_file,header=FALSE) %>% 
  rename('gene' = V1)

res_rbp <- as.data.frame(res) %>% 
  dplyr::mutate(gene=gsub("ENSG[1234567890]+[.][1234567890]+_", "",count_data$gene)) %>% 
  inner_join(known_rbp, by="gene") %>% mutate(Class="RBP")
res_tf <- as.data.frame(res) %>% 
  dplyr::mutate(gene=gsub("ENSG[1234567890]+[.][1234567890]+_", "",count_data$gene)) %>% 
  inner_join(known_tf, by="gene") %>% mutate(Class="TF")
res_kinase <- as.data.frame(res) %>% 
  dplyr::mutate(gene=gsub("ENSG[1234567890]+[.][1234567890]+_", "",count_data$gene)) %>% 
  inner_join(known_kinase, by="gene") %>% mutate(Class="kinase")
res_epi <- as.data.frame(res) %>% 
  dplyr::mutate(gene=gsub("ENSG[1234567890]+[.][1234567890]+_", "",count_data$gene)) %>% 
  inner_join(known_epi, by="gene") %>% mutate(Class="Epi")
res_ts <- as.data.frame(res) %>% 
  dplyr::mutate(gene=gsub("ENSG[1234567890]+[.][1234567890]+_", "",count_data$gene)) %>% 
  inner_join(known_ts, by="gene") %>% mutate(Class="TS")
res_onco <- as.data.frame(res) %>% 
  dplyr::mutate(gene=gsub("ENSG[1234567890]+[.][1234567890]+_", "",count_data$gene)) %>% 
  inner_join(known_onco, by="gene") %>% mutate(Class="Onco")
res_brain <- as.data.frame(res) %>% 
  dplyr::mutate(gene=gsub("ENSG[1234567890]+[.][1234567890]+_", "",count_data$gene)) %>% 
  inner_join(brain_goi_data, by="gene") %>% mutate(Class="Brain_GOI")

res_all_regl_df <- rbind(res_rbp, res_tf, res_kinase, res_epi, res_brain)
EnhancedVolcano(res_all_regl_df,
                lab = res_all_regl_df$gene, 
                x = 'log2FoldChange',
                y = 'pvalue',
                drawConnectors = TRUE,
                #ylim = c(0,21),
                #xlim = c(-3,3),
                title = 'Treated vs Ctrl',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3,
                labSize = 3)

ggsave(
  paste0(plots_dir,"/", "ctrl_vs_clk1-morp_volcano.goi.pdf"),
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

## write to file
write_delim(
  as_tibble(res_all_regl_df),
  paste0(results_dir,"/", de_reg_output),
  delim = "\t"
)


unlink(file.path("Rplots.pdf"))