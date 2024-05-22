################################################################################
# 09-clk1-nf1-protein-correlations.R
# Generate heatmap of correlation coefficients between CLK1 features (PSI, Expr) and CLK1 and NF1 RNA, protein, and phosphoprotein expression
# written by Ryan Corbett, Jo Lynne Rokita
# Usage: Rscript 09-clk1-nf1-protein-correlations.R
################################################################################

## load libraries
suppressPackageStartupMessages({
  library("ggplot2")
  library("tidyverse")
  library("ComplexHeatmap")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing_correlations")
results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")
figures_dir <- file.path(root_dir, "figures")

if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

# theme for all plots
source(file.path(figures_dir, "theme_for_plots.R"))

## define input files
cohort_file <- file.path(root_dir, "analyses", "cohort_summary",
                         "results", "histologies-plot-group.tsv")

data_file <- file.path(results_dir, "hgg-dmg-clk-nf1-expression-phosphorylation.tsv")
rsem_transc_counts <- file.path(data_dir,"rna-isoform-expression-rsem-tpm.rds")
rsem_tpm_counts <- file.path(data_dir,"gene-expression-rsem-tpm-collapsed.rds")


# output file
heatmap_output_file <- file.path(plots_dir, "clk1-nf1-single-sample-exp-protein-heatmap.pdf")

## read in histology, cohort, and independent specimens file
cohort_df <- read_tsv(cohort_file)


# define lists of nf1 and clk genes
goi <- c("CLK1", "NF1")
clk1_trans_list <- c("CLK1-201", "Total CLK1")
clk1_splice_list <- "CLK1-Exon4_PSI"
nf1_trans_list <- c("Total NF1", "NF1-202", "NF1-215")
nf1_splice_list <- c("NF1-Exon23a_PSI", "NF1-215_PSI")

trans_lists <- c(clk1_trans_list, nf1_trans_list)


# read in data file, zscored protein and psi values, but only keep all protein samples
data_df <- read_tsv(data_file) %>%
  filter(!is.na(NF1)) %>%
  select(match_id, all_of(nf1_splice_list), NF1, `NF1-S864`, `NF1-S2796`) %>%
  # zscore psi values
  mutate(`NF1-Exon23a_PSI` = scale(`NF1-Exon23a_PSI`),
         `NF1-215_PSI` = scale(`NF1-215_PSI`))

samples_oi <- cohort_df %>%
  filter(experimental_strategy == "RNA-Seq",
         RNA_library == "stranded",
         match_id %in% data_df$match_id)

## tpm table 
gene_tpm <- readRDS(rsem_tpm_counts) %>%
  dplyr::select(any_of(samples_oi$Kids_First_Biospecimen_ID)) %>%
  rownames_to_column(var = "gene_symbol") %>%
  filter(gene_symbol %in% goi) %>%
  mutate(gene_symbol = paste0("Total ", gene_symbol)) %>%
  column_to_rownames(var = "gene_symbol") %>%
  t() %>%
  as.data.frame()

# transcript tpm table
transcript_tpm_df <- readRDS(rsem_transc_counts) %>% 
  filter(gene_symbol %in% trans_lists) %>% 
  select(gene_symbol, any_of(samples_oi$Kids_First_Biospecimen_ID)) %>%
  column_to_rownames(var = "gene_symbol") %>%
  t() %>%
  as.data.frame() %>%
  bind_cols(gene_tpm) %>%
  rownames_to_column("Kids_First_Biospecimen_ID") %>%
  left_join(samples_oi[,c("match_id", "Kids_First_Biospecimen_ID")]) %>%
  select(-Kids_First_Biospecimen_ID)
  
exp_z <- transcript_tpm_df %>%
  # zscore
  mutate(across(-match_id, scale))

# join data
full_data <- data_df %>%
  left_join(exp_z) %>%
  as.data.frame() %>%
  column_to_rownames(var = "match_id") %>%
  t() %>%
  as.data.frame()

mat <- full_data %>%
  as.matrix()

storage.mode(mat) <- "numeric"
class(mat)

anno_names <- names(full_data) %>%
  as.data.frame()
names(anno_names) <- "match_id"

# create single sample heatmap
row_annot <- anno_names %>%
  left_join(unique(cohort_df[,c("match_id", "plot_group", "plot_group_hex")])) %>%
  column_to_rownames("match_id") %>%
  select(plot_group) %>%
  dplyr::rename(Histology = plot_group) %>%
  as.data.frame()

# create anno colors
anno_col <- list(Histology = c("DIPG or DMG" = "#ff40d9", "Other high-grade glioma" = "#ffccf5"))

# Heatmap annotation
row_anno = rowAnnotation(df = row_annot,
                         col = anno_col, show_legend = TRUE)

ha = HeatmapAnnotation(
  Histology = row_annot$Histology, col = anno_col,
  annotation_name_side = "right")

# Make heatmap without legends
heat_plot <- Heatmap(mat,
                     name = "Value",
                     col = colorRamp2(c(-2, 0, 2), c("#E66100", "white", "#5D3A9B")),
                     cluster_rows = FALSE,
                     #row_split = row_annot$Assay,
                     column_gap = 0.5,
                     show_row_names = TRUE,
                     show_column_names = FALSE,
                     show_heatmap_legend=TRUE,
                     cluster_columns = TRUE,
                    # right_annotation = row_anno,
                     top_annotation = ha,
                     #na_col = "lightgrey",
                     #rect_gp = gpar(col = "white"),
                     row_title = NULL,
                     column_title = NULL,
                     column_title_side = "top")

pdf(heatmap_output_file, width = 6, height = 4)
print(heat_plot)
dev.off()










