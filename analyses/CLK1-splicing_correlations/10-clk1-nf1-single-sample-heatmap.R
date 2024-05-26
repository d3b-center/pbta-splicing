################################################################################
# 10-clk1-nf1-single-sample-heatmap.R
# Generate heatmap of CLK1 features (PSI, Expr) and CLK1 and NF1 RNA, protein, and phosphoprotein expression
# written by Jo Lynne Rokita
# Usage: Rscript 10-clk1-nf1-single-sample-heatmap.R
################################################################################

## load libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
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
heatmap_output_file <- file.path(plots_dir, "clk1-nf1-single-sample-exp-protein-heatmap-dmg.pdf")

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
samples_oi <- cohort_df %>%
  filter(experimental_strategy == "RNA-Seq",
         RNA_library == "stranded",
         plot_group == "DIPG or DMG")

data_df <- read_tsv(data_file) %>%
  filter(!is.na(NF1)) %>%
  select(match_id, all_of(nf1_splice_list), all_of(clk1_splice_list), NF1, `NF1-S864`, `NF1-S2796`) %>%
  # zscore psi values
  mutate(`NF1-202 Exon23a PSI` = scale(`NF1-Exon23a_PSI`),
         `NF1-215 PSI` = scale(`NF1-215_PSI`),
         `CLK1-201 Exon4 PSI` = scale(`CLK1-Exon4_PSI`)) %>%
  select(-c(`CLK1-Exon4_PSI`, `NF1-Exon23a_PSI`, `NF1-215_PSI`)) %>%
  filter(match_id %in% samples_oi$match_id)

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
  # make labels pretty
  dplyr::rename(`Total NF1 protein` = NF1,
                `Total NF1 mRNA` = `Total NF1`,
              #  `NF1 pS2796` = `NF1-S2796`,
                `NF1 pS864` = `NF1-S864`) %>%
  select(-`Total CLK1`, -`NF1-S2796`) %>%
  # order by CLK1 exon 4 psi
  arrange(`CLK1-201 Exon4 PSI`) %>%
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
col_annot <- anno_names %>%
  left_join(unique(cohort_df[,c("match_id", "plot_group", "plot_group_hex")])) %>%
  column_to_rownames("match_id") %>%
  select(plot_group) %>%
  dplyr::rename(Histology = plot_group) %>%
  as.data.frame()

# create anno colors
col_anno_col <- list(Histology = c("DIPG or DMG" = "#ff40d9"))

ha = HeatmapAnnotation(
  Histology = col_annot$Histology, col = col_anno_col,
  annotation_name_side = "right")


# Heatmap annotation
row_annot <- as.data.frame(rownames(full_data)) %>%
  dplyr::rename(ID = `rownames(full_data)`) %>%
  mutate(Abundance = case_when(ID == "Total NF1 protein" ~ "Whole Cell Protein",
                           ID %in% c("NF1 pS864", "NF1 pS2796") ~ "Phosphoprotein",
                           ID %in% c("NF1-202", "CLK1-201", "NF1-215") ~ "Transcript",
                           ID == "Total NF1 mRNA" ~ "Gene",
                           ID %in% c("CLK1-201 Exon4 PSI", "NF1-202 Exon23a PSI", "NF1-215 PSI") ~ "mRNA splicing"))

row_anno_col <- list(Abundance = c("Transcript" = "#DC3220",
                                   "Gene" = "purple",
                                   "mRNA splicing" = "orange",
                                   "Phosphoprotein" = "#005AB5", 
                                   "Whole Cell Protein" = "#40B0A6"))

row_anno = rowAnnotation(Abundance = row_annot$Abundance,
                         col = row_anno_col, show_legend = TRUE)

# Make heatmap
heat_plot <- Heatmap(mat,
                     name = "Z-score",
                     col = colorRamp2(c(-2, 0, 2), c("darkblue", "white", "red")),
                     #col = colorRamp2(c(-2, 0, 2), c("#E66100", "white", "#5D3A9B")),
                     cluster_rows = TRUE,
                     column_gap = 0.5,
                     show_row_names = TRUE,
                     show_column_names = FALSE,
                     show_heatmap_legend=TRUE,
                     cluster_columns = FALSE,
                     clustering_distance_rows = "spearman",
                     right_annotation = row_anno,
                     top_annotation = ha,
                     na_col = "lightgrey",
                     rect_gp = gpar(col = "white"),
                     row_title = NULL,
                     column_title = NULL,
                     column_title_side = "top")

pdf(heatmap_output_file, width = 8, height = 4)
print(heat_plot)
dev.off()










