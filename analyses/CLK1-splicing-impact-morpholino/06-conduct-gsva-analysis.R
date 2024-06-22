################################################################################
# 05-conduct-gsva-analysis.R
# Performs gsva analysis on cells treated with control morpholino 
# or morpholinos targeting CLK1
#
# Author: Jo Lynne Rokita
# usage: Rscript --vanilla 06-conduct-gsva-analysis.R
################################################################################

## Load and/or install libraries ##
library(tidyverse)
library(readr)
library(tibble)
library(msigdbr) ## Contains the hallmark data sets
library(GSVA)    ## Performs GSEA analysis

# Magrittr pipe
`%>%` <- dplyr::`%>%`

#### Set Up paths and file names --------------------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing-impact-morpholino")
input_dir <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

## input/output files
expression_data_file <- file.path(data_dir, "ctrl_vs_morpho.rsem.genes.results.tsv")
expression_collapsed_file <- file.path(results_dir, "ctrl_vs_morpho.rsem.genes.collapsed.rds")
de_results_file <- file.path(results_dir, "ctrl_vs_treated.de.tsv")
splice_func_file <- file.path(results_dir,"differential_splice_by_goi_category.tsv")

# dna repair gene lists
dna_all_file <- file.path(input_dir, "dna_repair_all.txt")
hr_file <- file.path(input_dir, "homologous_recombination.txt")
mmr_file <- file.path(input_dir, "mismatch_repair.txt")
ber_file <- file.path(input_dir, "base_excision_repair.txt")
ner_file <- file.path(input_dir, "nucleotide_excision_repair.txt")
nej_file <- file.path(input_dir, "nonhomologous_end_joining.txt")

#### Load input file and collapse genes-----------------------------------------
expression_data <- read_tsv(expression_data_file) %>%
  # split to get symbol
  extract(col = gene, 
          into = c("ENS_ID", "Gene_Symbol"), 
          regex = "^(ENSG[0-9]+\\.[0-9]+)_(.+)$",
          remove = FALSE) %>%
  dplyr::select(-gene, -ENS_ID) %>%
  # remove PAR_Y_
  mutate(Gene_Symbol = str_replace_all(Gene_Symbol, "PAR_Y_", "")) %>%
  unique() 


# read in significant DE genes
de_results<- read_tsv(de_results_file) %>%
  filter(padj <0.05)

# read in sig splice events
splice_res_func <- read_tsv(splice_func_file) %>%
  dplyr::rename(Gene_Symbol = gene) %>%
  select(Gene_Symbol, Uniprot, plot_type) %>%
  unique()

# filter for onco/tsg
splice_res_onco <- splice_res_func %>%
  # filter for onc/tsg
  filter(plot_type %in% c("Oncogene", "Oncogene or Tumor Suppressor", "Tumor Suppressor"))

# remove all genes with no expression
expr <- expression_data[which(rowSums(expression_data[,2:ncol(expression_data)]) > 0),] 

# take mean per row and use the max value for duplicated gene symbols
expr_collapsed <- expr %>% 
  mutate(means = rowMeans(dplyr::select(.,-Gene_Symbol))) %>% # take rowMeans
  arrange(desc(means)) %>% # arrange decreasing by means
  distinct(Gene_Symbol, .keep_all = TRUE) %>% # keep the ones with greatest mean value. If ties occur, keep the first occurencce
  dplyr::select(-means) %>%
  unique() %>%
  column_to_rownames("Gene_Symbol") %>%
  write_rds(expression_collapsed_file)

# create a second matrix - only DE genes
expr_collapsed_de <- expr_collapsed[de_results$Gene_Symbol,] %>%
  na.omit()

# create a third matrix - only sig splice genes
expr_splice <- expr_collapsed[splice_res_func$Gene_Symbol,] %>%
  na.omit()

# create a fourth matrix - only sig spliced onco or tsgs
expr_splice_onco <- expr_collapsed[splice_res_onco$Gene_Symbol,] %>%
  na.omit()

# load gene lists
human_hallmark  <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") ## human hallmark genes from `migsdbr` package. The loaded data is a tibble.
human_kegg  <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") ## human hallmark genes from `migsdbr` package. The loaded data is a tibble.

#### Prepare genesets Create a list of sets, each of which is a list of genes -----------------------------------------------
# hallmark
human_hallmark_twocols <- human_hallmark %>% dplyr::select(gs_name, human_gene_symbol)
human_hallmark_list    <- base::split(human_hallmark_twocols$human_gene_symbol, list(human_hallmark_twocols$gs_name))

# kegg
human_kegg_twocols <- human_kegg %>% dplyr::select(gs_name, human_gene_symbol)
human_kegg_list    <- base::split(human_kegg_twocols$human_gene_symbol, list(human_kegg_twocols$gs_name))

# DNA repair from Knijnenburg paper https://pubmed.ncbi.nlm.nih.gov/29617664/
# Read gene lists from files and directly assign them as character vectors
dna_all_genes <- read_lines(dna_all_file)
hr_genes <- read_lines(hr_file)
mmr_genes <- read_lines(mmr_file)
ber_genes <- read_lines(ber_file)
ner_genes <- read_lines(ner_file)
nej_genes <- read_lines(nej_file)

# Create a dna_repair_list with the desired structure
dna_repair_list <- list(
  "DNA Repair all genes" = dna_all_genes,
  "Homologous Recombination" = hr_genes,
  "Mismatch Repair" = mmr_genes,
  "Base Excision Repair" = ber_genes,
  "Nucleotide Excision Repair" = ner_genes,
  "Non-homologous End Joining" = nej_genes
)

#### Perform gene set enrichment analysis --------------------------------------------------------------------

# Prepare expression data: log2 transform re-cast as matrix
### Rownames are genes and column names are samples

# list the matrices
mat_list <- list(expr_collapsed = expr_collapsed, 
                 expr_collapsed_de = expr_collapsed_de, 
                 expr_splice = expr_splice, 
                 expr_splice_onco = expr_splice_onco)
gsea_scores_df_tidy <- data.frame()

for(i in names(mat_list)) {
  expression_data_each_log2_matrix <- as.matrix(log2(mat_list[[i]] + 1))
  
  # Remove genes with 0 variance
  keep <- apply(expression_data_each_log2_matrix, 1, function(x) var(x, na.rm = TRUE)) > 0 
  expression_data_each_log2_matrix_keep <- expression_data_each_log2_matrix[keep, ]

  genesets <- list(human_hallmark_list = human_hallmark_list, 
                 human_kegg_list = human_kegg_list, 
                 dna_repair_list = dna_repair_list)
  names(genesets) <- c("human_hallmark_list", "human_kegg_list", "dna_repair_list")
  
  for (geneset_name in names(genesets)){  
   geneset <- genesets[[geneset_name]]
   if (geneset_name == "human_hallmark_list"){ 
     out_file <- file.path(results_dir, paste0(i, "_clk1_ctrl_morpho_hallmark_gsva_scores.tsv"))
   }
     else if (geneset_name == "human_kegg_list"){
      out_file <- file.path(results_dir, paste0(i, "_clk1_ctrl_morpho_kegg_gsva_scores.tsv"))
    }
    else if (geneset_name == "dna_repair_list"){
      out_file <- file.path(results_dir, paste0(i, "_clk1_ctrl_morpho_dna_repair_gsva_scores.tsv"))
    }
    
  #We then calculate the Gaussian-distributed scores
   gsea_scores_param <- gsvaParam(expression_data_each_log2_matrix,
                                geneSets = geneset,
                                kcdf = "Gaussian",
                                assay = NA_character_,
                                annotation = NA_character_,
                                tau = 1,
                                minSize = 1, 
                                maxSize = 1500, ## Arguments from K. Rathi
                                maxDiff = TRUE) ## Setting this argument to TRUE computes Gaussian-distributed scores (bimodal score distribution if FALSE)
  gsea_scores_each <- gsva(gsea_scores_param, verbose = TRUE)
    
  ### Clean scoring into tidy format
  gsea_scores_each_df <- as.data.frame(gsea_scores_each) %>%
    rownames_to_column(var = "geneset") %>%
    tidyr::pivot_longer(cols = -geneset, 
                        names_to = "sample_id", 
                        values_to = "gsva_score") %>%
    dplyr::select(sample_id, geneset, gsva_score) %>%
    unique()
  
  #### Export GSEA scores to TSV --------------------------------------------------------------------
  write_tsv(gsea_scores_each_df, out_file)
    }
}


#### Session info
sessionInfo()

