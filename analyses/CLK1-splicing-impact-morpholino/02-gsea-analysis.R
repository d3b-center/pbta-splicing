################################################################################
# 02-gsea-analysis.R
# Gene Set Enrichment Analysis using ClusterProfiler
#
# Authors: Ammar Naqvi and Shehbeel Arif
# Script adapted from ALSF's DE Analysis Tutorial (https://alexslemonade.github.io/refinebio-examples/03-rnaseq/pathway-analysis_rnaseq_02_gsea.html)
################################################################################

## load libraries
suppressPackageStartupMessages({
  library("clusterProfiler")
  library("msigdbr")
  library("org.Hs.eg.db")
  library("ggplot2")
  library("DOSE")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`



## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
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
gsea_dotplot_path <- file.path(plots_dir, "CLK1_status_gsea_dotplot.pdf")

# input file paths
dge_results_file <- file.path(results_dir, "ctrl_vs_treated.de.formatted.tsv")
dge_df <- readr::read_tsv(dge_results_file) 

## specify MSigDB gene sets of interest
hs_msigdb_df <- msigdbr(species = "Homo sapiens")
pathway_df <- hs_msigdb_df %>%
  dplyr::filter(gs_cat == "H" | gs_subcat %in% c("CP:KEGG", "CP:BIOCARTA", "TFT:GTRD"))

## determine our pre-ranked genes list
## create a named vector ranked based on the log2 fold change values
lfc_vector <- dge_df$log2FoldChange
names(lfc_vector) <- dge_df$gene

## sort  log2 fold change values in descending order 
lfc_vector <- sort(lfc_vector, decreasing = TRUE)

## run GSEA using the GSEA() function
## set the seed so our results are reproducible:
set.seed(1234)

#run GSEA function
gsea_results <- GSEA(
  geneList = lfc_vector, # Ordered ranked gene list
  minGSSize = 10, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    pathway_df,
    gs_name,
    gene_symbol
  )
)

## convert GSEA results object to dataframe
gsea_result_df <- data.frame(gsea_results@result)

# save GSEA results
readr::write_csv(
  gsea_result_df,
  file.path(
    results_dir,
    "clk1_morph_gsea_results.csv"
  )
)

gsea_plot <- dotplot(gsea_results, showCategory=15, split=".sign") + facet_grid(.~.sign) + 
  theme_Publication() +
  scale_color_gradient(name = "Adjusted p-value", 
                       low = "orange", high = "#0C7BDC") +  # Modify color range
  labs(color = "B-H adj p-value")  # Modify legend title 


# Save GSEA dotplot as pdf
ggplot2::ggsave(gsea_dotplot_path,
                plot = gsea_plot,
                width=14,
                height=10,
                device="pdf")



