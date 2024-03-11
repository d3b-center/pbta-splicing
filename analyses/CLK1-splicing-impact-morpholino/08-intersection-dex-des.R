################################################################################
# 08-intersection-dex-des.R
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
  library('ggVennDiagram')
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

## theme for all plots
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

##input files
de_file <- file.path(results_dir,"ctrl_vs_treated.de.tsv")
rmats_merged_file  <- file.path(data_dir,"morpholno.merged.rmats.tsv")

## outplut file for plot
ora_dotplot_func_path <- file.path(plots_dir, "CLK1_ds-dex-targets_ora_dotplot.pdf")
venn_output_file <- file.path(plots_dir, "des-dex-venn.pdf")

## extract strong splicing changes
splicing_df  <-  vroom(rmats_merged_file, comment = "#", delim="\t") %>%
  filter(FDR < 0.05 & PValue < 0.05)

## extract strong differential splicing cases (dPSI >= |.10|)
psi_comb <- splicing_df %>% 
  mutate(Preference = case_when(IncLevelDifference  >= .10 ~ "Skipping",
                                IncLevelDifference <= -.10 ~ "Inclusion",
                                TRUE ~ NA_character_),
         abs_IncLevelDifference = abs(IncLevelDifference)) %>%
  filter(!is.na(Preference))

dex_comb  <-  read_tsv(de_file) %>%
  mutate(Preference = case_when(log2FoldChange > 2 & padj < 0.05 ~ "Up",
                                log2FoldChange < 2 & padj < 0.05 ~ "Down",
                                TRUE ~ NA_character_)) %>%
  filter(!is.na(Preference)) %>%
  dplyr::rename(geneSymbol = Gene_Symbol)

intersect <- psi_comb %>%
  inner_join(dex_comb, by='geneSymbol', relationship = "many-to-many",
             suffix = c("_psi", "_de")) %>%
  dplyr::select(geneSymbol, Preference_psi, Preference_de) %>%
  unique

total_events <- psi_comb %>%
  full_join(dex_comb, by='geneSymbol', relationship = "many-to-many",
            suffix = c("_psi", "_de")) %>%
  dplyr::select(geneSymbol, Preference_psi, Preference_de) %>% 
  unique()

## plot venn diagram
venn_diag<- ggVennDiagram(x=list(dex_comb$geneSymbol, psi_comb$geneSymbol), 
                          edge_lty = "dashed", 
                          edge_size = 1,
                          label_size = 6,
                          set_size = 5,
                          category.names = c("DE" , "DS"),
                          label_percent_digit = 1) +  
  scale_fill_distiller(palette = "Blues", direction = 1, name = expression(bold("Gene count"))) + 
  labs(title = expression(bold("Genes dysregulated by splicing and expression")))

ggplot2::ggsave(venn_output_file,
                plot=venn_diag,
                width=6,
                height=5,
                device="pdf",
                dpi=300)


## get gene sets relevant to H. sapiens
hs_msigdb_df <- msigdbr(species = "Homo sapiens")
pathway_df <- hs_msigdb_df %>%
  dplyr::filter(gs_cat == "H" | gs_subcat %in% c("CP:KEGG"))


## run enrichR to compute and identify significant over-repr pathways
ora_results <- enricher(
  gene = unique(intersect$geneSymbol), # A vector of your genes of interest
  pvalueCutoff = 1, 
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

