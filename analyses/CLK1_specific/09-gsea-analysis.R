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
  library(vroom)
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
dge_results_file <- file.path(data_dir, "gene-counts-rsem-expected_count-collapsed.rds")

# Read in the contents of the differential expression results file
dge_df <- readr::read_rds(dge_results_file) %>%
  as.data.frame()


# Load histology data
histology_df <- readr::read_tsv(file.path(data_dir, "histologies.tsv"), guess_max = 100000)

# Filter histology dataframe
histology_rna_df <- histology_df %>% 
  # Only include RNA-Seq samples 
  filter(experimental_strategy == "RNA-Seq") %>% 
  # Remove NAs
  filter(!is.na(RNA_library)) %>%
  # Only include samples from PBTA cohort
  filter(cohort == "PBTA") %>%
  # Include samples that are "stranded"
  filter(RNA_library == "stranded") %>% 
  # Include those samples originated from the Midline
  filter(CNS_region == 'Midline') %>% 
  # Include only HGAT sample histology
  filter(short_histology == 'HGAT') %>%
  # Filter for module-specific samples that was previously used in CLK1 splicing/comparisons
  filter(Kids_First_Biospecimen_ID %in% c('BS_Q13FQ8FV', 'BS_ZV1P6W9C', 'BS_WH8G4VFB', 
                                          'BS_NNPEC7W1', 'BS_PZVHMSYN', 'BS_DRY58DTF', 
                                          'BS_GXTFW99H', 'BS_E60JZ9Z3', 'BS_9CA93S6D'))

rmats_file <-  file.path(data_dir,"rMATS_merged.single.SE.tsv.gz")
rmats_df <-  vroom(rmats_file, comment = "#",delim="\t") %>%
  # select specific samples and extract CLK1 exon 4  
  dplyr::filter(geneSymbol=="CLK1") %>% 
  dplyr::filter(exonStart_0base=="200860124", exonEnd=="200860215") %>% 
  dplyr::select(sample, geneSymbol, IncLevel1) %>% 
  dplyr::inner_join(histology_rna_df, by=c('sample'='Kids_First_Biospecimen_ID')) %>%  
  dplyr::select(sample, geneSymbol, IncLevel1) 

# Subset expression data to include only samples from filtered histology filter dataset
dge_df <- dge_df %>%
  dplyr::select(histology_rna_df$Kids_First_Biospecimen_ID)

# For each type of the RNA library, we subset the expression matrix accordingly and run gsea scores for each RNA library 
rna_library_list <- histology_rna_df %>% 
  dplyr::pull(RNA_library) %>% 
  unique()

# Further subset to each cohort to deal with size issues
cohort_list <- histology_rna_df %>% 
  dplyr::pull(cohort) %>% 
  unique()



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













