################################################################################
# 09-gsea-analysis.R
# Gene Set Enrichment Analysis using ClusterProfiler
# Author: Shehbeel Arif
# Script adapted from ALSF's DE Analysis Tutorial (https://alexslemonade.github.io/refinebio-examples/03-rnaseq/pathway-analysis_rnaseq_02_gsea.html)
################################################################################

## LOAD LIBRARIES
suppressPackageStartupMessages({
  # Library to perform GSEA and make nice visualizations
  library(clusterProfiler)
  # Package that contains MSigDB gene sets in tidy format
  library(msigdbr)
  # Human annotation package we'll use for gene identifier conversion
  library(org.Hs.eg.db)
  # We will need this so we can use the pipe: %>%
  library(magrittr)
  library(ggplot2)
  library(ggridges)
})

## SET DIRECTORIES
# Set input directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data", "v2")
analysis_dir <- file.path(root_dir, "analyses", "CLK1_specific")

# Set output directories
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

# Make output directories if they don't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Declare input file paths
dge_results_file <- file.path(results_dir, "clk1_group_diff_expr_results.tsv")

# Read in the contents of the differential expression results file
dge_df <- readr::read_tsv(dge_results_file) 

#######
# Specifying MSigDB gene sets of interest
hs_hallmark_sets <- msigdbr(
  species = "Homo sapiens", # Replace with species name relevant to your data
  category = "H"
)

#######
## PERFORM GSEA

# 1. Determine our pre-ranked genes list
# Check if there are any duplicate genes present
any(duplicated(dge_df$Gene)) # None

# Create a named vector ranked based on the log2 fold change values
lfc_vector <- dge_df$log2FoldChange
names(lfc_vector) <- dge_df$Gene

# Sort the log2 fold change values in descending order here
lfc_vector <- sort(lfc_vector, decreasing = TRUE)

# 2. Run GSEA using the GSEA() function
# Set the seed so our results are reproducible:
set.seed(1234)
# Run GSEA
gsea_results <- GSEA(
  geneList = lfc_vector, # Ordered ranked gene list
  minGSSize = 25, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    hs_hallmark_sets,
    gs_name,
    gene_symbol
  )
)


# We can access the results from our `gsea_results` object using `@result`
head(gsea_results@result)

# Convert GSEA results object to dataframe
gsea_result_df <- data.frame(gsea_results@result)

# Save GSEA results
readr::write_csv(
  gsea_result_df,
  file.path(
    results_dir,
    "clk1_status_gsea_results.csv"
  )
)

# 3. Visualize GSEA
# Look at the 3 gene sets with the most positive NES
gsea_result_df %>%
  # This returns the 3 rows with the largest NES values
  dplyr::slice_max(NES, n = 3)

# Make positive NES GSEA plot
most_positive_nes_plot <- enrichplot::gseaplot2(
  gsea_results,
  geneSetID = "HALLMARK_E2F_TARGETS",
  title = "HALLMARK_E2F_TARGETS"
)
most_positive_nes_plot

# Save GSEA enrichment plot as tiff
# Plot path
most_positive_nes_plot_path <- file.path(plots_dir, "CLK1_status_gsea_most_positive_plot.tiff")
# Save plot
ggplot2::ggsave(most_positive_nes_plot_path,
                width=10.7,
                height=5.33,
                device="tiff",
                dpi=300
)


# Look at the 3 gene sets with the most negative NES
gsea_result_df %>%
  # Return the 3 rows with the smallest (most negative) NES values
  dplyr::slice_min(NES, n = 3)

# Make negative NES GSEA plot
most_negative_nes_plot <- enrichplot::gseaplot2(
  gsea_results,
  geneSetID = "HALLMARK_KRAS_SIGNALING_DN",
  title = "HALLMARK_KRAS_SIGNALING_DOWN"
)
most_negative_nes_plot


# Save GSEA enrichment plot as tiff
# Plot path
most_negative_nes_plot_path <- file.path(plots_dir, "CLK1_status_gsea_most_negative_plot.tiff")
# Save plot
ggplot2::ggsave(most_negative_nes_plot_path,
                width=10.7,
                height=5.33,
                device="tiff",
                dpi=300
)


# Make GSEA enrichment plot for DNA Repair geneset
# Make DNA Repair NES GSEA plot
dna_repair_plot <- enrichplot::gseaplot2(
  gsea_results,
  geneSetID = "HALLMARK_DNA_REPAIR",
  title = "HALLMARK_DNA_REPAIR"
)
dna_repair_plot

# Save DNA Repair GSEA enrichment plot as tiff
# Plot path
dna_repair_plot_path <- file.path(plots_dir, "clk1_status_dna_repair_gsea_plot.tiff")
# Save plot
ggplot2::ggsave(dna_repair_plot_path,
                width=10.7,
                height=5.33,
                device="tiff",
                dpi=300
)


# Make dotplot of enriched gene sets
require(DOSE)
dotplot(gsea_results, showCategory=10, split=".sign") + facet_grid(.~.sign)

# Save GSEA dotplot as tiff
# Plot path
gsea_dotplot_path <- file.path(plots_dir, "CLK1_status_gsea_dotplot.tiff")
# Save plot
ggplot2::ggsave(gsea_dotplot_path,
                width=10.7,
                height=5.33,
                device="tiff",
                dpi=300
)


# Make ridgeplot of enriched gene sets
ridgeplot(gsea_results) + labs(x = "enrichment distribution")

# Save GSEA dotplot as tiff
# Plot path
gsea_ridgeplot_path <- file.path(plots_dir, "CLK1_status_gsea_ridgeplot.tiff")
# Save plot
ggplot2::ggsave(gsea_ridgeplot_path,
                width=10.7,
                height=8,
                device="tiff",
                dpi=300
)

dev.off()

