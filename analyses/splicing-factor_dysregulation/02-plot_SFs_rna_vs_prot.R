################################################################################
# 02-plot_SFs_rna_vs_prot.R
# written by Ammar Naqvi, Jo Lynne Rokita
#
# This script usesCPTAC data to generate and plot heatmap of RNA vs proteomics
# of select splicing factors that was previously identified in 
# 01-volcano_plot_mRNA.R script
#
# usage: Rscript 02-plot_SFs_rna_vs_prot.R
################################################################################

## load packages
suppressPackageStartupMessages({
  library(readxl)
  library(tidyverse)
  library(circlize)
  library(ComplexHeatmap)
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "splicing-factor_dysregulation")
data_dir   <- file.path(root_dir, "data")

input_dir   <- file.path(analysis_dir, "input")
plots_dir   <- file.path(analysis_dir, "plots")
results_dir   <- file.path(analysis_dir, "results")

## check and create plots dir
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

## file paths
heatmap_output_file <- file.path(plots_dir,"SF_RNA_vs_protein_levels_heatmap.pdf") 

## get CPTAC output table 
cptac_output_file <- file.path(input_dir,"CPTAC3-pbt_SFs.xls") 
de_file <- file.path(results_dir, "all_hgg-diffSFs_sig_genes.txt")

# Load dataset
de_genes <- read_tsv(de_file) %>%
  pull(gene)

cptac_data <- readxl::read_excel(cptac_output_file) %>%
  # select only rows with CLK1 or SRFs, remove muts
  filter(grepl("\\srna|\\spro", idx), 
               !`Data type` %in% c("mut", "cnv")) %>%
  # remove extra info from cols
  rename_with(~ gsub("X7316.", "7316-", .), everything()) %>%
  dplyr::rename(Assay = `Data type`) %>%
  filter(`Gene symbol` %in% de_genes) %>%
  # clean up naming for plotting
  mutate(Assay = case_when(Assay == "proteo" ~ "Whole Cell Proteomics",
                           Assay == "rna" ~ "RNA-Seq"),
         Assay = fct_relevel(Assay, c("RNA-Seq", "Whole Cell Proteomics")),
         # create new display name for phospho proteins to include phos site
         display_name = case_when(
                                  # add a space after the gene to trick into thinking the rownames are not duplicated
                                  Assay == "Whole Cell Proteomics" ~ paste(`Gene symbol`, " ", sep = " "),
                                  TRUE ~ `Gene symbol`)
  ) %>%
  dplyr::select(display_name, Assay, starts_with("7316")) %>%
  # remove NAs
  select_if(~ !any(is.na(.)))

# preserve gene names for rownames
rownames <- cptac_data$display_name

# convert to matrix and then add rownames
mat <- cptac_data %>%
  dplyr::select(3:(ncol(cptac_data))) %>%
  as.matrix()
storage.mode(mat) <- "numeric"
class(mat)


# set rownames, convert to matrix
rownames(mat) <- rownames

# select row annotations
row_annot <- cptac_data %>%
  dplyr::select(Assay) %>%
  as.data.frame()


# create anno colors
anno_col <- list(Assay = c("RNA-Seq" = "#DC3220", "Whole Cell Proteomics" = "#40B0A6"))

# Heatmap annotation
row_anno = rowAnnotation(df = row_annot,
                         col = anno_col, 
                         show_legend = TRUE, 
                         show_annotation_name = FALSE)

# Make heatmap without legends
heat_plot <- Heatmap(mat,
                     name = "Z-score",
                     col = colorRamp2(c(-2, 0, 2), c("#E66100", "white", "#5D3A9B")),
                     cluster_rows = FALSE,
                     row_split = row_annot$Assay, 
                     column_gap = 0.5,
                     show_row_names = TRUE,
                     show_column_names = FALSE,
                     show_heatmap_legend = TRUE,
                     cluster_columns = TRUE, 
                     right_annotation = row_anno,
                     row_title = NULL, 
                     column_title = NULL, 
                     row_names_gp = grid::gpar(fontsize = 9),
                     column_title_side = "top",
                     heatmap_legend_param = list(legend_direction = "horizontal", 
                                                 legend_position = "top"))

pdf(heatmap_output_file, width = 7, height = 5)
draw(heat_plot, heatmap_legend_side = "top")
dev.off()
