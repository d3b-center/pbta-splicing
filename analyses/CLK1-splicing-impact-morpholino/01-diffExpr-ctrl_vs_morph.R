################################################################################
# 01-diffExpr-ctrl_vs_morph.R
# Performs gene expression analysis and plots on cells untreated vs treated with
# morpholinos targeting CLK1
#
# Author: Ammar Naqvi, Jo Lynne Rokita
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
  library(ggplot2)
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data/")
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing-impact-morpholino")

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
de_output = file.path(results_dir, "ctrl_vs_treated.de.tsv")
file_volc_plot = file.path(plots_dir, "ctrl_vs_clk1-morp_volcano.pdf")
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
res$Significant <- ifelse(res$padj< 0.05, "P-val < 0.05", "Not Sig")

volc <- EnhancedVolcano(res,
                lab = gsub("ENSG[1234567890]+[.][1234567890]+_", "",count_data$gene), ## remove ensembleid portion
                x = 'log2FoldChange',
                y = 'padj',
                drawConnectors = TRUE,
                #ylim = c(0,21),
               # xlim = c(-6,6),
                title = 'CLK1 Exon 4 Morpholino vs Non-targeting Morpholino',
                caption = "",
                subtitle = "",
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2,
                labSize = 4) 
  # Attempt to override axis titles post-hoc
  volc <- volc + labs(x = expression(bold(Log[2] * " Fold Change")), 
                      y = expression(bold("-Log"[10] * " p-value")))

pdf(file_volc_plot, height = 7, width = 8)
print(volc)
dev.off()

de_results <- as_tibble(cbind(res, count_data)) %>%
  extract(col = gene, 
          into = c("ENS_ID", "Gene_Symbol"), 
          regex = "^(ENSG[0-9]+\\.[0-9]+)_(.+)$",
          remove = FALSE) %>%
  select(-gene)

write_tsv(de_results, de_output)

## subset with gene lists
# gene list file
oncokb_gene_file <-file.path(input_dir,'cancerGeneList.tsv')

## oncoKB gene list
oncokb_gene_ref <- vroom(oncokb_gene_file) %>% 
  dplyr::rename('gene' = `Hugo Symbol`) %>%
  mutate(
    type = case_when(
      `Is Oncogene` == "Yes" & `Is Tumor Suppressor Gene` == "Yes" ~ "Both",
      `Is Oncogene` == "Yes" ~ "Oncogene",
      `Is Tumor Suppressor Gene` == "Yes" ~ "Tumor Suppressor",
      TRUE ~ NA_character_
    )) %>%
  select(gene,type) %>% 
  filter(!is.na(type))

# combine with the DE results
genes_to_plot <- as.data.frame(res) %>% 
  dplyr::mutate(gene=gsub("ENSG[1234567890]+[.][1234567890]+_", "",count_data$gene)) %>% 
  inner_join(oncokb_gene_ref, by="gene")                          

## use only significant results for the plot
sign_regl_gene_df <- genes_to_plot %>%
  as.data.frame() %>%
  dplyr::filter(padj < 0.05,
         abs(log2FoldChange) >= 1) %>%
  mutate(Direction= case_when(log2FoldChange < -1 ~ 'Up',
                              log2FoldChange > 1 ~ 'Down')) 


## plot num of hits per gene fam
plot_barplot_family <- ggplot(sign_regl_gene_df, aes(x = fct_rev(fct_infreq(type)), fill= Direction)) +
                       geom_bar(stat="count", position='dodge', color="black") + 
                       facet_wrap(~type, scales = "free_y", ncol = 1) +
                       xlab("Cancer Gene Type")     + 
                       ylab("Number of Genes Signficantly DE") + 
                       scale_fill_manual(name = "Direction (CLK1 exon 4 high)",
                                         values=c("#FFC20A","#0C7BDC")) + 
                       geom_text(stat='count',aes(label=after_stat(count)), 
                                 position = position_dodge(width = 1),
                                 hjust = -0.5, size = 3.5) +
                      theme_Publication() +
                      theme(legend.position = "top", legend.direction = "horizontal") +
                      coord_flip() 

  

# print and save plot
pdf(file_gene_family_plot, height = 6, width = 8, useDingbats = FALSE) 
print(plot_barplot_family)
dev.off()
