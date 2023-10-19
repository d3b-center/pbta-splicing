################################################################################
# 03-plot_SR-phosp_vs_CLK1-RNA.R
# written by Shehbeel Arif, Ammar Naqvi
#
# This script uses CPTAC data to generate and plot heatmap of RNA vs phosp of 
# CLK1 and SRSF phospho-levels
#
# usage: Rscript 03-plot_SR-phosp_vs_CLK1-RNA.R
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
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing_correlations")

input_dir   <- file.path(analysis_dir, "input")
plots_dir   <- file.path(analysis_dir, "plots")

## check and create plots dir
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

## output files for final plots
heatmap_output_file <- file.path(plots_dir,"heatmap_SR-phosp_CLK1-RNA_rna_prot.tiff") 

## get CPTAC output table 
cptac_output_file <- file.path(input_dir,"CPTACT_SR_CLK1.xls") 

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
  filter(type %in% c("rna", "phospho"))

mat <- df %>%
  select(3:220) %>% 
  mutate_if(is.character, as.numeric)

rownames(mat) <- make.unique(as.character(c("CLK1","SRSF10","SRSF11","SRSF11","SRSF11","SRSF2","SRSF2","SRSF3","SRSF4","SRSF6","SRSF8","SRSF9")))


row_annot <- df %>%
  select("type")
rownames(row_annot) <- rownames(mat)

# Specify colors
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                                      "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
                                      
annotation_colors = list(type = c(rna="#FFC20A", phospho="#0C7BDC"))

pheatmap(mat = mat,
         annotation_row = row_annot,
         gaps_row = c(1,1,1),
         cluster_cols = TRUE,
         cluster_rows = FALSE,
         show_colnames = FALSE,
         show_rownames = TRUE,
         legend_labels = c("Z-Score"),
         legend = TRUE,
         cellwidth = 12,
         fontsize = 25,
         height = 20,
         width = 42,
         #annotation = row_annot,
         annotation_colors  = annotation_colors,
         filename = heatmap_output_file)