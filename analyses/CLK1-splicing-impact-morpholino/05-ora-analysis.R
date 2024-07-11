################################################################################
# 05-ora-analysis.R
# Over-represenative analysis of mis-spliced genes that have splicing variants
# impacting functional sitesmediated by CLK1
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
  library("tidyverse")
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

## theme for all plots
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## outplut file for plot
ora_dotplot_func_path <- file.path(plots_dir, "CLK1_targets_ora_dotplot.func-sites.pdf")

## get gene sets relevant to H. sapiens
hs_msigdb_df <- msigdbr(species = "Homo sapiens")
pathway_df <- hs_msigdb_df %>%
  dplyr::filter(gs_cat == "H" | gs_subcat %in% c("CP:KEGG"))

## extract splicing changes 
events_func_file <- file.path(results_dir,"splicing_events.morpho.intersectUnip.ggplot.txt")

## ORA on functionally relevant splice variants
events_func_df  <-  vroom(events_func_file, comment = "#", delim="\t") %>% 
  mutate(geneSymbol=str_match(SpliceID, "(\\w+[\\.\\d]*)\\:")[, 2]) %>%
  dplyr::select(geneSymbol) %>% 
  unique()

## run enrichR to compute and identify significant over-repr pathways
ora_results <- enricher(
  gene = events_func_df$geneSymbol, # A vector of your genes of interest
  pvalueCutoff = 1, 
  pAdjustMethod = "BH", 
  TERM2GENE = dplyr::select(
    pathway_df,
    gs_name,
    human_gene_symbol
  )
)

ora_result_df <- data.frame(ora_results@result)
<<<<<<< HEAD
enrich_plot_func<- enrichplot::dotplot(ora_results) +   
  theme_Publication() +
  scale_color_gradient(name = "Adjusted p-value", 
                       low = "orange", high = "#0C7BDC") +  # Modify color range
  labs(color = "B-H adj p-value")  # Modify legend title 

ggplot2::ggsave(ora_dotplot_func_path,
                plot=enrich_plot_func,
                width=9,
                height=5,
=======
options(enrichplot.colours = c("darkorange","blue"))
enrich_plot_func <- enrichplot::dotplot(ora_results,
                                        x = "geneRatio",
                                        size = "Count",
                                        color = "p.adjust",
                                        label_format = 30,
                                        showCategory = 20) +   
  labs(y = "Pathway",
       x = "Gene Ratio") +
  theme_Publication() +
  scale_size(name = "Gene Count") +  
  scale_fill_gradient(low = "darkorange", high = "blue", name = "B-H p-value") +
  guides(
    fill = guide_colorbar(title = "B-H p-value", label.position = "right", barwidth = 1, barheight = 4)
  )

ggplot2::ggsave(ora_dotplot_func_path,
                plot=enrich_plot_func,
                width=8,
                height=4,
>>>>>>> main
                device="pdf",
                dpi=300)
