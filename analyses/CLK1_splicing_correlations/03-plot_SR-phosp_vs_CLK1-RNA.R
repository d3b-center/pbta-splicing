################################################################################
# 03-plot_SR-phosp_vs_CLK1-RNA.R
# written by Shehbeel Arif, Ammar Naqvi, Jo Lynne Rokita
#
# This script uses CPTAC data to generate and plot heatmap of RNA vs phosp of 
# CLK1 and SRSF phospho-levels
#
# usage: Rscript --vanilla 03-plot_SR-phosp_vs_CLK1-RNA.R
################################################################################

## load packages
suppressPackageStartupMessages({
  library(readxl)
  library(circlize)
  library(tidyverse)
  library(ComplexHeatmap)
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "CLK1_splicing_correlations")
data_dir   <- file.path(root_dir, "data")

input_dir   <- file.path(analysis_dir, "input")
plots_dir   <- file.path(analysis_dir, "plots")

## check and create plots dir
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

## file paths
heatmap_output_file <- file.path(plots_dir,"SR_phos_CLK1_exp_heatmap.pdf") 

## get CPTAC output table 
cptac_output_file <- file.path(input_dir,"CPTAC_SR_CLK1.xls") 
cptac_output_file2 <- file.path(input_dir,"CPTAC3-pbt.xls") 

# Load dataset
cptac_data <- readxl::read_excel(cptac_output_file2) %>%
  # select only rows with CLK1 or SRFs, remove muts
  filter(grepl("CLK1|SRSF", idx),
         !`Data type` %in% c("mut", "cnv")) %>%
  # remove extra info from cols
  rename_with(~ gsub("X7316.", "7316-", .), everything()) %>%
  dplyr::rename(Assay = `Data type`) %>%
  # create new display name for phospho proteins
  mutate(Assay = case_when(Assay == "proteo" ~ "Whole Cell Proteomics",
                           Assay == "phospho" ~ "Phospho-Proteomics",
                           Assay == "rna" ~ "RNA-Seq"),
         Assay = fct_relevel(Assay, c("RNA-Seq", "Whole Cell Proteomics", "Phospho-Proteomics")),
         phos_site = case_when(Assay == "Phospho-Proteomics" ~ str_split(idx, "phospho", simplify = TRUE)[, 2],
                               TRUE ~ NA_character_),
         display_name = case_when(Assay == "Phospho-Proteomics" ~ paste(`Gene symbol`, phos_site, sep = " "),
                                  # add a space after the gene to trick into thinking the rownames are not duplicated
                                  Assay == "Whole Cell Proteomics" ~ paste(`Gene symbol`, " ", sep = " "),
                                  TRUE ~ `Gene symbol`)
  ) %>%
  select(display_name, Assay, starts_with("7316")) %>%
  select_if(~ !any(is.na(.)))

# preserve gene names for rownames
rownames <- cptac_data$display_name

# convert to matrix and then add rownames
mat <- cptac_data %>%
  select(3:(ncol(cptac_data))) %>%
  as.matrix()
storage.mode(mat) <- "numeric"
class(mat)


# set rownames, convert to matrix
rownames(mat) <- rownames
# mat <- as.matrix(mat)

# select row annotations
row_annot <- cptac_data %>%
  select(Assay) %>%
  as.data.frame()

# add rownames
rownames(row_annot) <- rownames(mat)


# create anno colors
anno_col <- list(Assay = c("RNA-Seq" = "#DC3220", "Phospho-Proteomics" = "#005AB5", "Whole Cell Proteomics" = "#40B0A6"))

# Heatmap annotation
row_anno = rowAnnotation(df = row_annot,
                         col = anno_col, show_legend = TRUE)


# Make heatmap without legends
heat_plot <- Heatmap(mat,
                     name = "Z-score",
                     col = colorRamp2(c(-2, 0, 2), c("#E66100", "white", "#5D3A9B")),
                     cluster_rows = FALSE,
                     row_split = row_annot$Assay, 
                     column_gap = 0.5,
                     show_row_names = TRUE,
                     show_column_names = FALSE,
                     show_heatmap_legend=TRUE,
                     cluster_columns = TRUE, 
                     right_annotation = row_anno,
                     na_col = "lightgrey",
                     #rect_gp = gpar(col = "white"),
                     row_title = NULL, 
                     column_title = NULL, 
                     column_title_side = "top")

heat_plot

pdf(heatmap_output_file, width = 8, height = 6)
print(heat_plot)
dev.off()
