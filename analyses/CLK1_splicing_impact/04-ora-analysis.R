################################################################################
# 04-ora-analysis.R
# Over-represenative analysis of mis-spliced genes mediated by CLK1
#
# authors: Ammar Naqvi 
################################################################################

## load libraries
suppressPackageStartupMessages({
  library("clusterProfiler")
  library("msigdbr")
  library("org.Hs.eg.db")
  library("ggplot2")
  library("ggridges")
  library("DOSE")
  library("vroom")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1_splicing_impact")

input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")

## create folder if non-existent 
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive=TRUE)
}

## outplut file for plot
ora_dotplot_path <- file.path(plots_dir, "CLK1_targets_ora_dotplot.tiff")

## get gene sets relevant to H. sapiens
hs_msigdb_df <- msigdbr(species = "Homo sapiens")

## filter for hallmark pathways that are included in the curated gene sets
hs_hm_df <- hs_msigdb_df %>%
  dplyr::filter(
    gs_cat == "H"
  )

## extract splicing changes 
rmats_merged_file  <- file.path(analysis_dir,"input","morpholno.merged.rmats.tsv")
splicing_df  <-  vroom(rmats_merged_file, comment = "#", delim="\t") %>% 
  filter(FDR < 0.05 & PValue < 0.05) 

## create a background set of mis-spliced genes
background_set <- splicing_df %>% dplyr::select(geneSymbol) %>% unique()

## create a genes of interest list based on strong splicing 
genes_of_interest <- splicing_df %>% 
                              filter( (IncLevelDifference >= .10) | (IncLevelDifference <= -.10)) %>% 
                              dplyr::select(geneSymbol) %>% 
                             unique()


## run enrichR to compute and identify significant over-repr pathways
ora_results <- enricher(
  gene = genes_of_interest$geneSymbol, # A vector of your genes of interest
  pvalueCutoff = 0.05, 
  pAdjustMethod = "BH", 
  TERM2GENE = dplyr::select(
    hs_hm_df,
    gs_name,
    human_gene_symbol
  )
)

## plot enrichment using dotplot
ora_result_df <- data.frame(ora_results@result)
enrich_plot <- enrichplot::dotplot(ora_results)
enrich_plot

## save ORA dotplot as tiff
ggplot2::ggsave(ora_dotplot_path,
                width=7,
                height=5,
                device="tiff",
                dpi=300)
                
  