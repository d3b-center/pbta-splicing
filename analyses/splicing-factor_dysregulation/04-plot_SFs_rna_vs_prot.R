################################################################################
# 04-plot_SFs_rna_vs_prot.R
# written by Shehbeel Arif
#
# This script usesCPTAC data to generate and plot heatmap of RNA vs proteomics
# of select splicing factors that was previously identified in 
# 01-volcano_plot_mRNA.R script
#
# usage: Rscript 04-plot_SFs_rna_vs_prot.R
################################################################################

## load packages
suppressPackageStartupMessages({
  library("readxl")
  library("pheatmap")
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
plots_dir   <- file.path(analysis_dir, "plots")

## output files for final plots
heatmap_output_file <- file.path(plots_dir,"heatmap_SFs_rna_prot.tiff") 

## get CPTAC output table 
cptac_output_file <- file.path(input_dir,"CPTAC3-pbt.xls") 

# The problem with the dataframe was that it contained both 
# "NA" the string and NA for null, so I first made all "NA"s 
# as NA and then replaced them with 0.

# Load dataset
df <- read_excel(cptac_output_file) 
# Replace string "NA" with NA
df[df == "NA"] <- NA

# Convert dataframe to matrix and replace NAs with 0
df <- df %>% 
  as.matrix() %>% 
  replace(is.na(.), 0) %>% 
  as.data.frame() %>% 
  # Only look at RNA and Protein data
  filter(type %in% c("rna", "proteo"))

mat <- df %>%
  select(3:220) %>% 
  mutate_if(is.character, as.numeric) %>% 
  t() 
colnames(mat) <- paste0("row_", seq(ncol(mat)))

col_annot <- df %>%
  select("type", "gene")
rownames(col_annot) <- colnames(mat)

#col_annot$type <- factor(col_annot$type, levels = c("proteo", "rna"))
#col_annot$gene <- factor(col_annot$gene, levels = c("MSI1", "NOVA2", "RBM15", "SAMD4A"))

pheatmap(mat = mat,
         annotation_col = col_annot,
         gaps_col = c(2,4,6),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_colnames = FALSE,
         show_rownames = TRUE,
         height = 10,
         width = 7,
         filename = heatmap_output_file
)
################################################################################
## SAMPLES as rows
# Load dataset
df <- read_excel(cptac_output_file) 
# Replace string "NA" with NA
df[df == "NA"] <- NA

# Convert dataframe to matrix and replace NAs with 0
df <- df %>% 
  as.matrix() %>% 
  replace(is.na(.), 0) %>% 
  as.data.frame() %>% 
  # Only look at RNA and Protein data
  filter(type %in% c("rna", "proteo"))

mat <- df %>%
  select(3:220) %>% 
  mutate_if(is.character, as.numeric)

rownames(mat) <- paste0("col_", seq(nrow(mat)))

row_annot <- df %>%
  select("type", "gene")
rownames(row_annot) <- rownames(mat)

# Specify colors
annotation_colors = list(
  gene = c(MSI1="#994F00", NOVA2="#40B0A6", RBM15="#E66100", SAMD4A ="#5D3A9B"),
  type = c(rna="#FFC20A", proteo="#0C7BDC"))

pheatmap(mat = mat,
         annotation_row = row_annot,
         gaps_row = c(2,4,6),
         cluster_cols = TRUE,
         cluster_rows = FALSE,
         show_colnames = FALSE,
         show_rownames = FALSE,
         legend = TRUE,
         cellwidth = 8,
         fontsize = 30,
         height = 10,
         width = 30,
         #annotation = row_annot,
         annotation_colors  = annotation_colors,
         filename = heatmap_output_file
)
