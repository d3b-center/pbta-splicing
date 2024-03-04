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
  library("DOSE")
  library("vroom")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing-impact-morpholino")

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

## theme for all plots
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## outplut file for plot
ora_dotplot_inc_path <- file.path(plots_dir, "CLK1_targets_ora_dotplot.incl-in-high-Ex.pdf")
ora_dotplot_skp_path <- file.path(plots_dir, "CLK1_targets_ora_dotplot.skip-in-high-Ex.pdf")

## get gene sets relevant to H. sapiens
hs_msigdb_df <- msigdbr(species = "Homo sapiens")
pathway_df <- hs_msigdb_df %>%
  dplyr::filter(gs_cat == "H" | gs_subcat %in% c("CP:KEGG"))

## extract splicing changes 
rmats_merged_file  <- file.path(analysis_dir,"input","morpholno.merged.rmats.tsv")
events_func_file <- file.path(results_dir,"splicing_events.morpho.intersectUnip.ggplot.txt")

splicing_df  <-  vroom(rmats_merged_file, comment = "#", delim="\t") %>% 
  filter(FDR < 0.05 & PValue < 0.05) 



## create a genes of interest list based on strong splicing 
genes_of_interest_inc <- splicing_df %>% 
                              filter( (IncLevelDifference <= -.10)) %>% ## only get events that are high in Ex 4 (untreated) cells
                              dplyr::select(geneSymbol) %>% 
                              unique()

## create a genes of interest list based on strong splicing 
genes_of_interest_skip <- splicing_df %>% 
  filter( (IncLevelDifference >= .10)) %>% ## only get events that are high in Ex 4 (untreated) cells
  dplyr::select(geneSymbol) %>% 
  unique()

genes_of_interest_inc_uniq <- anti_join(genes_of_interest_inc, genes_of_interest_skip, by='geneSymbol' )
genes_of_interest_skip_uniq <- anti_join(genes_of_interest_skip,genes_of_interest_inc, by='geneSymbol')

bg_set<-  splicing_df %>% 
  dplyr::select(geneSymbol) %>% 
  unique()


## run enrichR to compute and identify significant over-repr pathways
ora_results <- enricher(
  gene = genes_of_interest_inc_uniq$geneSymbol, # A vector of your genes of interest
  pvalueCutoff = 0.05, 
  #universe = bg_set$geneSymbol,
  pAdjustMethod = "BH", 
  TERM2GENE = dplyr::select(
    pathway_df,
    gs_name,
    human_gene_symbol
  )
)

## plot enrichment using dotplot
ora_result_df <- data.frame(ora_results@result)
enrich_plot_incl<- enrichplot::dotplot(ora_results) +   
  theme_Publication() +
  scale_color_gradient(name = "Adjusted p-value", 
                       low = "orange", high = "#0C7BDC") +  # Modify color range
  labs(color = "B-H adj p-value")  # Modify legend title 


## save ORA dotplot as tiff
ggplot2::ggsave(ora_dotplot_inc_path,
                plot=enrich_plot_incl,
                width=9,
                height=5,
                device="pdf",
                dpi=300)

## run enrichR to compute and identify significant over-repr pathways
ora_results <- enricher(
  gene = genes_of_interest_skip$geneSymbol, # A vector of your genes of interest
  pvalueCutoff = 0.05, 
  #universe = bg_set$geneSymbol,
  pAdjustMethod = "BH", 
  TERM2GENE = dplyr::select(
    pathway_df,
    gs_name,
    human_gene_symbol
  )
)

ora_result_df <- data.frame(ora_results@result)
enrich_plot_skip<- enrichplot::dotplot(ora_results) +   
  theme_Publication() +
  scale_color_gradient(name = "Adjusted p-value", 
                       low = "orange", high = "#0C7BDC") +  # Modify color range
  labs(color = "B-H adj p-value")  # Modify legend title 
                
ggplot2::ggsave(ora_dotplot_skp_path,
                plot=enrich_plot_skip,
                width=9,
                height=5,
                device="pdf",
                dpi=300)

## focus on functional impactful events
events_func_df  <-  vroom(events_func_file, comment = "#", delim="\t") %>% 
  mutate(geneSymbol=str_match(SpliceID, "(\\w+[\\.\\d]*)\\:")[, 2]) %>%
  dplyr::select(geneSymbol) %>% 
  unique()

## run enrichR to compute and identify significant over-repr pathways
ora_results <- enricher(
  gene = events_func_df$geneSymbol, # A vector of your genes of interest
  pvalueCutoff = 0.05, 
  #universe = bg_set$geneSymbol,
  pAdjustMethod = "BH", 
  TERM2GENE = dplyr::select(
    pathway_df,
    gs_name,
    human_gene_symbol
  )
)

ora_result_df <- data.frame(ora_results@result)
enrich_plot_func<- enrichplot::dotplot(ora_results) +   
  theme_Publication() +
  scale_color_gradient(name = "Adjusted p-value", 
                       low = "orange", high = "#0C7BDC") +  # Modify color range
  labs(color = "B-H adj p-value")  # Modify legend title 

ggplot2::ggsave(ora_dotplot_func_path,
                plot=enrich_plot_func,
                width=9,
                height=5,
                device="pdf",
                dpi=300)
