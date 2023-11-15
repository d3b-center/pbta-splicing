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
                                    "midline_hggs_v_ctrl_SFs_volcano.pdf")

gene_sign_list_file <- file.path(results_dir,"midline_hggs_v_ctrl_SFs_sig_genes.txt")

## input files
sf_file <- file.path(input_dir,"splicing_factors.txt")
clin_file <- file.path(data_dir,"histologies.tsv")
file_gene_counts <- file.path(data_dir,"gene-counts-rsem-expected_count-collapsed.rds")
file_nontumor_count <- file.path(input_dir,"gene_counts_normals_final.csv")


## get splicing factor list to subset later
sf_list <- read_lines(sf_file)

## get clinical histlogy file filtered by HGG samples
clin_tab <- read_tsv(clin_file, guess_max = 100000) %>% 
  filter(short_histology == 'HGAT',
         RNA_library == 'stranded',
         cohort == 'PBTA',
         CNS_region == 'Midline')

## get gene count table with splicing factors and midline HGGs filter
count_data <- readRDS(file_gene_counts) %>%
  #filter for HGG midline samples
  select(any_of(clin_tab$Kids_First_Biospecimen_ID))

##  filter for only splicing factors from above SF list
count_data_sf <- count_data[rownames(count_data) %in% sf_list, ] 
count_data_sf <- count_data_sf %>%
  mutate(gene = rownames(count_data_sf)) %>%
  # rearrange
  select(gene, any_of(clin_tab$Kids_First_Biospecimen_ID))

## get corresponding non-brain tumor samples
gene_counts_nontumor  <-  read_csv(file_nontumor_count) %>% 
  mutate(gene=str_replace(gene, "ENSG[1234567890]+_", "") )

gene_counts_combined <- inner_join(count_data_sf, gene_counts_nontumor, by = 'gene')

# filter for rows with >= 10 count
filtered_counts <- gene_counts_combined[rowSums(gene_counts_combined>2) == 62, ]

filtered_counts <- gene_counts_combined %>% 
  filter(sum(c_across(where(is.numeric))) >= 620) %>%
  ungroup

## construct metadata
design <- data.frame(row.names = colnames(filtered_counts$gene),
                    condition = c(rep("HGG",53), rep("Control",9)),
                    libType   = c(rep("paired-end",62))) %>%
  # the DESEQ analysis will automatically use alphabetical order, so relevel to get the correct comparison
  mutate(condition = fct_relevel(condition, c("HGG", "Control")))

## remove first column
filtered_counts_gene_rm <- select(filtered_counts, -gene)

cds = DESeqDataSetFromMatrix(countData=round(filtered_counts_gene_rm),
                             colData=design,
                             design= ~ condition)

## run deseq function to compute pvalues
cds <- DESeq(cds)

## check levels so they are in the right order with HGG first
levels(cds$condition)

# get results
res <- results(cds)

## label anything below <0.05 as signficant
res$Significant <- ifelse(res$padj < 0.05, "P-val < 0.05", "Not Sig")
res$gene <- filtered_counts$gene

  volc <- EnhancedVolcano(res,
                  lab = res$gene, # Use the new label column
                  subtitle = "",
                  x = 'log2FoldChange',
                  y = 'pvalue',
                  xlab = expression(bold("log"[2]*" Fold Change")),
                  ylab = expression(bold("-log"[10]*" p-value")),
                 # ylim = c(0,21),
                # xlim = c(-3,3),
                  title = 'Midline HGG vs. Bainstem control',
                  drawConnectors = TRUE,
                  pCutoff = 0.05,
                  FCcutoff = 1,
                  pointSize = 2,
                  labSize = 4) 
# print plot
pdf(file_volc_hgg_SF_plot, height = 7, width = 7, useDingbats = FALSE)
print(volc)
dev.off()

## write significant genes to table for subsequent correlation analyses
gene_sign_list <- res %>%
  as.data.frame() %>%
  filter(padj < 0.05,
         abs(log2FoldChange) > 1) %>%
  select(gene, everything(res)) %>%
write_tsv(gene_sign_list_file)

