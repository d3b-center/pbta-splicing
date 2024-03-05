################################################################################
# 01-plot-diffExp_highlowSBI.R
#
# generates volcanto plot of differentail expression between high vs low 
# splicing burden HGG tumors 
#
# written by Ammar Naqvi
#
# usage: Rscript 01-plot-diffExp_highlowSBI.R
################################################################################

## libraries needed
suppressPackageStartupMessages({
  library("dplyr")
  library("EnhancedVolcano")
  library("DESeq2")
  library(ggplot2)
})


# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "splicing-factor_dysregulation")

input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

## theme for all plots
# source function for theme for plots survival
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## output files for final plots and results
file_volc_hgg_plot <- file.path(analysis_dir, "plots", "enhancedVolcano_hgg_sbi.pdf")
file_barplot_SFs_plot <- file.path(analysis_dir, "plots", "barplot_hgg_SFs.pdf")

gene_sign_list_file <- file.path(results_dir,"diffSFs_sig_genes.txt")

## get and setup input files
sbi_coding_file  <- file.path(analysis_dir, "splicing_index/results/splicing_index.SE.txt")
clin_file <- file.path(data_dir, "histologies.tsv")
file_gene_counts <- file.path(data_dir,"gene-counts-rsem-expected_count-collapsed.rds")

# get splicing factor list to subset later
sf_file <- file.path(analysis_dir, "splicing-factor_dysregulation/input/splicing_factors.txt")
sf_list <- readr::read_lines(sf_file)

## get clinical histlogy file filtered by HGG samples
clin_tab <- readr::read_tsv(clin_file, guess_max = 100000) %>% 
  filter(short_histology == 'HGAT',
         RNA_library == 'stranded',
         cohort == 'PBTA')

# read in files, join palette with sbi file
sbi_coding_hgg_df  <-  readr::read_tsv(sbi_coding_file, comment = "#") %>%
  dplyr::rename(Kids_First_Biospecimen_ID = Sample) %>%
  filter(Kids_First_Biospecimen_ID %in% clin_tab$Kids_First_Biospecimen_ID)

## compute quartiles to define high vs low SBI tumors
quartiles_sbi <- quantile(sbi_coding_df$SI, probs=c(.25, .75), na.rm = FALSE)
lower_sbi <- quartiles_sbi[1]
upper_sbi <- quartiles_sbi[2] 

# annotate as high/low SBI
sbi_coding_hgg_df <- sbi_coding_hgg_df %>%
  mutate(SI_level = case_when(SI > upper_sbi ~ "High",
                              SI < lower_sbi ~ "Low",
                              TRUE ~ NA_character_)) %>%
  filter(!is.na(SI_level))
sbi_coding_hgg_df$SI_level <- factor(sbi_coding_hgg_df$SI_level, levels = c("Low", "High"))


## get gene count table with midline HGGs filter
count_data <- readRDS(file_gene_counts) %>% 
  #filter for HGG midline samples stranded and high sbi
  dplyr::select(any_of(sbi_coding_hgg_df$Kids_First_Biospecimen_ID)) 

# Add gene names as a column to count_data
count_data <- count_data %>%
  mutate(gene = rownames(.))

# Filter count_data based on sf_list, then select specified columns
count_data_sf <- count_data %>%
  filter(gene %in% sf_list) %>%
  select(gene, any_of(clin_tab$Kids_First_Biospecimen_ID)) %>%
  rowwise() %>%  # Ensure you use parentheses here
  filter(sum(c_across(where(is.numeric))) >= 1) %>%
  ungroup()

## remove first column
filtered_counts_gene_rm <- dplyr::select(count_data_sf, -gene)

## construct metadata
design = data.frame(row.names = sbi_coding_hgg_df$Kids_First_Biospecimen_ID,
                    condition = sbi_coding_hgg_df$SI_level)

condition = sbi_coding_hgg_df$SI_level

cds = DESeqDataSetFromMatrix(countData=round(filtered_counts_gene_rm),
                             colData=design,
                             design= ~ condition)

cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)
cds <- DESeq(cds)

res <- results(cds)
res$Significant <- ifelse(res$padj< 0.05, "P-Adj < 0.05", "Not Sig") 
res <- as_tibble(res) %>% 
  tibble::add_column(gene=count_data_sf$gene) 

volc_hgg_plot <- EnhancedVolcano(res,
                                 lab = gsub("ENSG[1234567890]+[.][1234567890]+_", "",count_data_sf$gene), ## remove ensembleid portion
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 xlim = c(-4, 6.5),
                                 title = 'High vs Low SBI HGGs',
                                 subtitle = NULL,
                                 caption = NULL,
                                 pCutoff = 0.005,
                                 FCcutoff = .5,
                                 pointSize = 4,
                                 labSize = 5,
                                 typeConnectors = "closed",
                                 drawConnectors = TRUE,
                                 widthConnectors = 1,
                                 colConnectors = 'black')

# Attempt to override axis titles post-hoc
volc_hgg_plot <- volc_hgg_plot + labs(x = expression(bold(Log[2] * " Fold Change")), 
                    y = expression(bold("-Log"[10] * " p-value")))


## write significant genes to table for subsequent correlation analyses
gene_sign_list <- res %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>%
  dplyr::select(gene, everything(res)) %>%
  readr::write_tsv(gene_sign_list_file)

## plot and focus on the two major splicing factor families (well known control exon-splicing)
plot_df <- gene_sign_list %>% filter(grepl("SRSF|HNRNP", gene)) %>% 
  mutate(Direction= case_when(log2FoldChange<0 ~ '-',
                              log2FoldChange>0 ~ '+')) 

plot_barplot_family <-  ggplot(plot_df, aes(x = reorder(gene,-padj), y = -log2(padj))) + 
  geom_bar(stat="identity", colour="black", fill="red") + 
  theme_Publication() + 
  xlab("Splicing Factor") + ylab("-log2 (padj)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  coord_flip() + 
  geom_text(aes(label =paste(Direction),ymax=0), 
            hjust = -0.5, size = 4) +
  ylim(c(0,45))

# Save plots as PDF
pdf(file_volc_hgg_plot, 
    width = 8, height = 8)
volc_hgg_plot
dev.off()

# Save plot as PDF
pdf(file_barplot_SFs_plot, 
    width = 4, height = 6)
plot_barplot_family
dev.off()
