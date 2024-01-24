################################################################################
# 01-diffExpr-ctrl_vs_morph.R
# Performs gene expression analysis and plots on cells untreated vs treated with
# morpholinos targeting CLK1
#
# Author: Ammar Naqvi 
# usage: Rscript --vanilla 01-diffExpr-ctrl_vs_morph.R
################################################################################

suppressPackageStartupMessages({
  library("EnhancedVolcano")
  library("DESeq2")
  library("tidyverse")
  library("dplyr")
  library("vroom")
  library("clusterProfiler")
  library("annoFuseData")
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


## theme for all plots
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## output files
de_output = "ctrl_vs_treated.de.tsv"
de_reg_output = "ctrl_vs_treated.de.regl.tsv"
file_volc_plot = "ctrl_vs_clk1-morp_volcano.pdf"
file_gene_family_plot = file.path(plots_dir,"gene-fam-DE-plot.pdf")

## input files
# count data
tpm_count_file <- "ctrl-vs-morpholino-gene-counts-rsem-expected_count.tsv" 

# gene list files
known_rbp_file <- file.path(input_dir,'RBP_known.txt')
known_epi_file <- file.path(input_dir,'epi_known.txt')

## hgnc file for liftover
hgnc_file <- file.path(input_dir,"hgnc_complete_set.txt")

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
genelist_ref_df <-read.delim(system.file("extdata", "genelistreference.txt", package = "annoFuseData")) 

known_kinase <- genelist_ref_df %>% 
  mutate(Class = if_else(grepl("Kinase", type), 'Kinase', NA)) %>%
  na.omit() %>%
  dplyr::select(Gene_Symbol,Class) 

known_tf<- genelist_ref_df %>% 
  mutate(Class = if_else(grepl("TranscriptionFactor", type), "TF", NA)) %>%
  na.omit() %>%
  dplyr::select(Gene_Symbol,Class) 
  
known_ts <- genelist_ref_df %>% 
    mutate(Class = if_else(grepl("TumorSuppressorGene", type), "TS", NA)) %>%
    na.omit() %>%
  dplyr::select(Gene_Symbol,Class) 
    
known_onco <- genelist_ref_df %>% 
  mutate(Class = if_else(grepl("Oncogene", type), "Onco", NA)) %>%
  na.omit() %>%
  dplyr::select(Gene_Symbol,Class) 

## other non-annoFuse gene lists, RBPs and epigenetic genes
known_rbp <- read.table(known_rbp_file,header=FALSE) %>% 
  dplyr::rename('Gene_Symbol' = V1)
known_epi <- read.table(known_epi_file,header=FALSE) %>% 
  dplyr::rename('Gene_Symbol' = V1)

res_rbp <- as.data.frame(res) %>% 
  dplyr::mutate(gene=gsub("ENSG[1234567890]+[.][1234567890]+_", "",count_data$gene)) %>% 
  dplyr::rename('Gene_Symbol'=gene) %>% 
  inner_join(known_rbp, by="Gene_Symbol") %>%
  mutate(Class="RBP")


res_tf <- as.data.frame(res) %>% 
  dplyr::mutate(gene=gsub("ENSG[1234567890]+[.][1234567890]+_", "",count_data$gene)) %>% 
  dplyr::rename('Gene_Symbol'=gene) %>% 
  inner_join(known_tf, by="Gene_Symbol") %>% 
  mutate(Class="TranscriptionFactor")

res_kinase <- as.data.frame(res) %>% 
  dplyr::mutate(gene=gsub("ENSG[1234567890]+[.][1234567890]+_", "",count_data$gene)) %>% 
  dplyr::rename('Gene_Symbol'=gene) %>% 
  inner_join(known_kinase, by="Gene_Symbol") %>% 
  mutate(Class="Kinase")

res_epi <- as.data.frame(res) %>% 
  dplyr::mutate(gene=gsub("ENSG[1234567890]+[.][1234567890]+_", "",count_data$gene)) %>% 
  dplyr::rename('Gene_Symbol'=gene) %>% 
  inner_join(known_epi, by="Gene_Symbol") %>% 
  mutate(Class="Epigenetic")

res_ts <- as.data.frame(res) %>% 
  dplyr::mutate(gene=gsub("ENSG[1234567890]+[.][1234567890]+_", "",count_data$gene)) %>% 
  dplyr::rename('Gene_Symbol'=gene) %>% 
  inner_join(known_ts, by="Gene_Symbol") %>% 
  mutate(Class="TumorSuppressor")

res_onco <- as.data.frame(res) %>% 
  dplyr::mutate(gene=gsub("ENSG[1234567890]+[.][1234567890]+_", "",count_data$gene)) %>% 
  dplyr::rename('Gene_Symbol'=gene) %>% 
  inner_join(known_onco, by="Gene_Symbol") %>% 
  mutate(Class="Oncogene")

## combine results and liftover geneSymbols to new names
liftover_df <- vroom(hgnc_file) %>%
  dplyr::select(symbol, prev_symbol)

res_all_regl_df <- rbind(res_rbp, res_tf, res_kinase, res_epi, res_ts,res_onco) %>%
  ## liftover geneSymbols
  left_join(liftover_df, by=c('Gene_Symbol'='prev_symbol')) %>% 
  mutate(Gene_Symbol = case_when(
    is.na(symbol) ~ Gene_Symbol,
    TRUE ~ symbol
  )) 

## write significant genes to table for subsequent correlation analyses
sign_regl_gene_df <- res_all_regl_df %>%
  as.data.frame() %>%
  dplyr::filter(padj < 0.05,
         abs(log2FoldChange) >= 1) %>%
  mutate(Direction= case_when(log2FoldChange<1 ~ 'Down',
                              log2FoldChange>1 ~ 'Up'))


## plot num of hits per gene fam
plot_barplot_family <- ggplot(sign_regl_gene_df, aes(x = fct_rev(fct_infreq(Class)), fill= Direction)) +
                       geom_bar(stat="count", position='dodge', color="black") + 
                       xlab("Gene Family")     + 
                       ylab("Number of Signficantly DE Genes") + 
                       scale_fill_manual(name = "Direction",
                                         values=c("#FFC20A","#0C7BDC")) + 
                       geom_text(stat='count',aes(label=after_stat(count)), 
                                 position = position_dodge(width = 1),
                                 hjust = -0.5, size = 4) +
                      theme_Publication() +
                      coord_flip() 
  

# print and save plot
pdf(file_gene_family_plot, height = 10.38, width = 6.73, useDingbats = FALSE) 
print(plot_barplot_family)
dev.off()