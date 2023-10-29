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
clin_file <- file.path(data_dir, "histologies.tsv")
heatmap_output_file <- file.path(plots_dir,"SR_phos_CLK1_exp_heatmap.pdf") 
tpm_file <- file.path(data_dir, "gene-expression-rsem-tpm-collapsed.rds")

# load clinical file
clin_rna <- read_tsv(clin_file, guess_max = 100000) %>%
  filter(cohort == "PBTA",
         experimental_strategy == "RNA-Seq")

# load TPM file for only cohort of interest
tpm_data <- readRDS(tpm_file) %>%
  select(clin_rna$Kids_First_Biospecimen_ID)

## get CPTAC output table 
cptac_output_file <- file.path(input_dir,"CPTAC_SR_CLK1.xls") 
cptac_output_file2 <- file.path(input_dir,"Gygi_TCMP_HMS_Phosphoproteome.phosphosite.tmt11_RearrangedCols.csv") 

# Load dataset
cptac_data <- read_excel(cptac_output_file, na = "NA")  %>%
  # Only look at RNA and Protein data
  filter(type %in% c("rna", "phospho")) %>%
  rename_with(~ gsub("X7316.", "7316-", .), everything()) %>%
  # remove any samples which have NA
  select_if(~ !any(is.na(.)))

cptac_data2 <- read_csv(cptac_output_file2)  %>%
  select(-c("Prot", "Peptide", "Organism")) %>%
  # remove extra info from cols
  rename_with(~ gsub(" Log Ratio|  Log Ratio", "", .), everything()) %>%
  filter(grepl("SRSF", Gene) & rowSums(!is.na(.[-(1:5)])) > 0) 
       
setdiff(names(cptac_data), names(cptac_data2))  
setdiff(names(cptac_data2), names(cptac_data)) 

cptac_data %>%
  filter(gene == "SRSF11") %>%
  select(gene, all_of(names(cptac_data[,3:8])))

cptac_data2 %>%
  filter(Gene == "SRSF11") %>%
  select(Gene, all_of(names(cptac_data[,3:8])))


# some genes have duplicate phos sites, so make them unique ids for later rownaming
cptac_data$suffix <- ave(cptac_data$gene, cptac_data$gene, FUN = function(x) 1 + seq_along(x) - 1)

cptac_data_suffixed <- cptac_data %>%
  mutate(gene = ifelse(duplicated(gene) | rev(duplicated(rev(gene))), glue::glue("{gene}-{suffix}"), gene))
# preserve gene names for rownames
rownames <- cptac_data_suffixed$gene

# convert to matrix and then add rownames
mat <- cptac_data_suffixed %>%
  select(3:(ncol(cptac_data_suffixed)-1))

# set rownames, convert to matrix
rownames(mat) <- rownames
mat <- as.matrix(mat)

# select row annotations
row_annot <- cptac_data_suffixed %>%
  select("type") %>%
  dplyr::rename(`Assay` = type) %>%
  as.data.frame()

# add rownames
rownames(row_annot) <- rownames(mat)
  
# create anno colors
anno_col <- list(`Assay` = c(rna = "#CD5C5C", phospho = "#0C7BDC"))

# Heatmap annotation
row_anno = rowAnnotation(df = row_annot,
                         col = anno_col, show_legend = TRUE)


# Make heatmap without legends
heat_plot <- Heatmap(mat,
                     name = "Z-score",
                     col = colorRamp2(c(-2, 0, 2), c("orange", "white", "purple")),
                     cluster_rows = FALSE,
                     row_split = row_annot$type, 
                     column_gap = 0.5,
                     show_row_names = TRUE,
                     show_column_names = FALSE,
                     show_heatmap_legend=TRUE,
                     cluster_columns = TRUE, 
                     right_annotation = row_anno,
                     #rect_gp = gpar(col = "white"),
                     row_title = NULL, 
                     column_title = NULL, 
                     column_title_side = "top")

heat_plot

pdf(heatmap_output_file, width = 7, height = 5)
print(heat_plot)
dev.off()

