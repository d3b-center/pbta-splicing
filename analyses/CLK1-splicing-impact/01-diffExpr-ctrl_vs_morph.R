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
res$Significant <- ifelse(res$pvalue< 0.05, "P-val < 0.05", "Not Sig")

volc <- EnhancedVolcano(res,
                lab = gsub("ENSG[1234567890]+[.][1234567890]+_", "",count_data$gene), ## remove ensembleid portion , ## remove ensembleid portion
                x = 'log2FoldChange',
                y = 'pvalue',
                drawConnectors = TRUE,
                #ylim = c(0,21),
                #xlim = c(-3,3),
                title = 'CLK1 Exon 4 Morpholino vs Non-targeting Morpholino',
                caption = "",
                subtitle = "",
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2,
                labSize = 4)

ggsave(file_volc_plot,
  plot = volc,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8,
  height = 10,
  units = "in",
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)

de_results <- cbind(res, count_data) 

write_tsv(as_tibble(de_results), de_output)

## subset with gene lists
genelist_ref_df <-read.delim(system.file("extdata", "genelistreference.txt", package = "annoFuseData"))  %>%
  # remove cosmic census
  filter(type != "CosmicCensus") %>%
  mutate(type = gsub("CosmicCensus, |, CosmicCensus", "", type))

## other non-annoFuse gene lists, RBPs and epigenetic genes
known_rbp_inlist <- read.table(known_rbp_file,header=FALSE) %>% 
  dplyr::rename('Gene_Symbol' = V1) %>%
  filter(Gene_Symbol %in% genelist_ref_df$Gene_Symbol)

known_rbp_not_inlist <- read.table(known_rbp_file,header=FALSE) %>% 
  dplyr::rename('Gene_Symbol' = V1) %>%
  filter(!Gene_Symbol %in% genelist_ref_df$Gene_Symbol) %>%
  mutate(plot_type = "Other",
         plot_subtype = "RNA Binding Protein")

known_epi_inlist <- read.table(known_epi_file,header=FALSE) %>% 
  dplyr::rename('Gene_Symbol' = V1) %>%
  filter(Gene_Symbol %in% genelist_ref_df$Gene_Symbol)

known_epi_not_inlist <- read.table(known_epi_file,header=FALSE) %>% 
  dplyr::rename('Gene_Symbol' = V1) %>%
  filter(!Gene_Symbol %in% genelist_ref_df$Gene_Symbol) %>%
  mutate(plot_type = "Other",
         plot_subtype = "Epigenetic")


genelist_cat <- genelist_ref_df %>%
  mutate(type = ifelse(Gene_Symbol %in% known_rbp_inlist$Gene_Symbol, paste(type, "RNA Binding Protein", sep = ", "), type),
         type = ifelse(Gene_Symbol %in% known_epi_inlist$Gene_Symbol, paste(type, "Epigenetic", sep = ", "), type)) %>%
  # add epi/rbp for those which are already in list and recategorize
  mutate(plot_type = case_when(grepl("TumorSuppressorGene, Oncogene", type) ~ "Oncogene or Tumor Suppressor",
                               type %in% c("TumorSuppressorGene, TranscriptionFactor", "TumorSuppressorGene, Kinase", "TumorSuppressorGene", 
                                           "Kinase, TumorSuppressorGene", "TumorSuppressorGene, RNA Binding Protein",
                                           "TumorSuppressorGene, TranscriptionFactor, RNA Binding Protein, Epigenetic",
                                           "TumorSuppressorGene, RNA Binding Protein, Epigenetic", 
                                           "TumorSuppressorGene, TranscriptionFactor, Epigenetic", "TumorSuppressorGene, Epigenetic",
                                           "TumorSuppressorGene, Kinase, Epigenetic", "Kinase, TumorSuppressorGene, Epigenetic") ~ "Tumor Suppressor",
                               type %in% c("Kinase, Oncogene", "Oncogene", "Oncogene, Kinase", "Oncogene, TranscriptionFactor",
                                           "Oncogene, RNA Binding Protein", "Oncogene, TranscriptionFactor, RNA Binding Protein",
                                           "Oncogene, TranscriptionFactor, Epigenetic", "Oncogene, Kinase, Epigenetic",
                                           "Oncogene, Epigenetic", "Kinase, Oncogene, Epigenetic", "Oncogene, RNA Binding Protein, Epigenetic",
                                           "TranscriptionFactor, Oncogene") ~ "Oncogene",
                               TRUE ~ "Other"),
         plot_subtype = case_when(type %in% c("Kinase", "Kinase, TranscriptionFactor", "Kinase, Oncogene", "TumorSuppressorGene, Kinase", 
                                              "Kinase, TumorSuppressorGene, Oncogene", "Oncogene, Kinase", "Kinase, RNA Binding Protein",
                                              "TumorSuppressorGene, Oncogene, Kinase, Epigenetic", "TumorSuppressorGene, Kinase, Epigenetic",
                                              "Oncogene, Kinase, Epigenetic", "Kinase, TumorSuppressorGene, Oncogene, Epigenetic",
                                              "Kinase, TumorSuppressorGene, Epigenetic", "Kinase, Oncogene, Epigenetic", "Kinase, Epigenetic",
                                              "Kinase, TumorSuppressorGene", "TumorSuppressorGene, Oncogene, Kinase") ~ "Kinase",
                                  type %in% c("Oncogene, TranscriptionFactor", "TumorSuppressorGene, Oncogene, TranscriptionFactor",
                                              "TumorSuppressorGene, TranscriptionFactor", "TranscriptionFactor, Oncogene") ~ "Transcription Factor",
                                  type %in% c("Oncogene, RNA Binding Protein", "TranscriptionFactor, RNA Binding Protein", 
                                              "TumorSuppressorGene, Oncogene, RNA Binding Protein", "TranscriptionFactor, RNA Binding Protein, Epigenetic",
                                              "TumorSuppressorGene, RNA Binding Protein, Epigenetic", "TumorSuppressorGene, TranscriptionFactor, RNA Binding Protein",
                                              "TumorSuppressorGene, TranscriptionFactor, RNA Binding Protein, Epigenetic",
                                              "Oncogene, RNA Binding Protein, Epigenetic", "Oncogene, TranscriptionFactor, RNA Binding Protein",
                                              "TumorSuppressorGene, RNA Binding Protein") ~ "RNA Binding Protein",
                                  type %in% c("TumorSuppressorGene, TranscriptionFactor, Epigenetic", "TumorSuppressorGene, Oncogene, TranscriptionFactor, Epigenetic",
                                              "TumorSuppressorGene, Oncogene, Epigenetic", "TumorSuppressorGene, Epigenetic",
                                              "TranscriptionFactor, Epigenetic", "Oncogene, TranscriptionFactor, Epigenetic",
                                              "Oncogene, Epigenetic") ~ "Epigenetic",
                                  type == "TumorSuppressorGene" ~ "Other Tumor Suppressor",
                                  type == "Oncogene" ~ "Other Oncogene",
                                  type == "TumorSuppressorGene, Oncogene" ~ "Oncogene or Tumor Suppressor",
                                  type == "TranscriptionFactor" ~ "Transcription Factor",
                                  TRUE ~ type)) %>%
  dplyr::select(Gene_Symbol, plot_type, plot_subtype) %>%
  # add other RBP, Epi not already in the list
  bind_rows(known_rbp_not_inlist, known_epi_not_inlist)

# check that all of the subgroups look right
table(genelist_cat$plot_subtype, genelist_cat$plot_type)

# combine with the DE results
genes_to_plot <- as.data.frame(res) %>% 
  dplyr::mutate(gene=gsub("ENSG[1234567890]+[.][1234567890]+_", "",count_data$gene)) %>% 
  dplyr::rename('Gene_Symbol'=gene) %>% 
  inner_join(genelist_cat, by="Gene_Symbol")                          

## use only significant results for the plot
sign_regl_gene_df <- genes_to_plot %>%
  as.data.frame() %>%
  dplyr::filter(padj < 0.05,
         abs(log2FoldChange) >= 1) %>%
  mutate(Direction= case_when(log2FoldChange<1 ~ 'Down',
                              log2FoldChange>1 ~ 'Up')) 

# relevel the plot_type and subtype
sign_regl_gene_df$plot_type <- factor(sign_regl_gene_df$plot_type, levels = c("Oncogene", "Tumor Suppressor",
                                                                              "Oncogene or Tumor Suppressor", "Other"))
sign_regl_gene_df$plot_subtype <- factor(sign_regl_gene_df$plot_subtype, levels = c("Kinase", "Oncogene or Tumor Suppressor",
                                                                              "Transcription Factor", "RNA Binding Protein",
                                                                              "Epigenetic", "Other Oncogene", "Other Tumor Suppressor"))
unique(sign_regl_gene_df$plot_subtype)
  
## plot num of hits per gene fam
plot_barplot_family <- ggplot(sign_regl_gene_df, aes(x = fct_rev(fct_infreq(plot_subtype)), fill= Direction)) +
                       geom_bar(stat="count", position='dodge', color="black") + 
                       facet_wrap(~plot_type, scales = "free_y", ncol = 1) +
                       xlab("Gene Family")     + 
                       ylab("Number of Genes Signficantly Differentially Expressed") + 
                       scale_fill_manual(name = "Direction",
                                         values=c("#FFC20A","#0C7BDC")) + 
                       geom_text(stat='count',aes(label=after_stat(count)), 
                                 position = position_dodge(width = 1),
                                 hjust = -0.5, size = 3.5) +
                      theme_Publication() +
                      coord_flip() 

  

# print and save plot
pdf(file_gene_family_plot, height = 8, width = 8, useDingbats = FALSE) 
print(plot_barplot_family)
dev.off()
