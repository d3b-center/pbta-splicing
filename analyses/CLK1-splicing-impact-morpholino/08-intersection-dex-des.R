################################################################################
# 08-intersection-dex-des.R
# Over-represenative analysis of mis-spliced genes that have splicing variants
# impacting functional sitesmediated by CLK1
#
# authors: Ammar Naqvi, Jo Lynne Rokita
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
categories_file <- file.path(results_dir, "gene_categories.tsv")
rmats_merged_file  <- file.path(data_dir,"morpholno.merged.rmats.tsv")
file_psi_func <- file.path(results_dir,"differential_splice_by_goi_category.tsv")

## output file for plot
ora_dotplot_path <- file.path(plots_dir, "CLK1_ds-dex-targets_ora_dotplot.pdf")
ora_dotplot_func_path <- file.path(plots_dir, "CLK1_ds-dex-targets_ora_dotplot-func.pdf")

venn_output_file <- file.path(plots_dir, "des-dex-venn.pdf")
venn_output_func_file <- file.path(plots_dir, "des-dex-venn-func.pdf")

# read in categories file
categories_df <- read_tsv(categories_file)

## extract strong splicing changes
splicing_df  <-  vroom(rmats_merged_file, comment = "#", delim="\t") %>%
  filter(FDR < 0.05 & PValue < 0.05)

splice_func_df <- read_tsv(file_psi_func)

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
  labs(title = expression(bold("Differentially expressed and spliced genes"))) +
  coord_flip()

ggplot2::ggsave(venn_output_file,
                plot=venn_diag,
                width=5.5,
                height=4,
                device="pdf",
                dpi=300)

## get gene sets relevant to H. sapiens
hs_msigdb_df <- msigdbr(species = "Homo sapiens")
pathway_df <- hs_msigdb_df %>%
  dplyr::filter(gs_cat == "H" | gs_subcat %in% c("CP:KEGG"))


## run enrichR to compute and identify significant over-repr pathways
ora_results <- enricher(
  gene = unique(intersect$geneSymbol), # A vector of your genes of interest
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  pAdjustMethod = "BH", 
  TERM2GENE = dplyr::select(
    pathway_df,
    gs_name,
    human_gene_symbol
  )
)

ora_result_df <- data.frame(ora_results@result)
options(enrichplot.colours = c("darkorange","blue"))
enrich_plot <- enrichplot::dotplot(ora_results,
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

ggplot2::ggsave(ora_dotplot_path,
                plot=enrich_plot,
                width=7,
                height=3,
                device="pdf",
                dpi=300)

## venn for functional splice sites
## read table of recurrent functional splicing (skipping)
splice_func_df <- splice_func_df %>%
  dplyr::rename(geneSymbol = gene) %>%
  # subset for goi
  filter(geneSymbol %in% categories_df$gene)

# subset DEX
dex_comb_subset <- dex_comb %>%
  filter(geneSymbol %in% categories_df$gene)

intersect_func <- splice_func_df %>%
  inner_join(dex_comb_subset, by='geneSymbol', relationship = "many-to-many") %>%
  pull(geneSymbol) %>%
  unique()

total_events_func <- splice_func_df %>%
  full_join(dex_comb_subset, by='geneSymbol', relationship = "many-to-many") %>%
  dplyr::select(geneSymbol) %>% 
  unique()

sign_targets_crispr_file <- file.path(results_dir,"sign-crispr-cbtn-lines.txt")
sign_targets_crispr <- readLines(sign_targets_crispr_file) 

# Define header names
header_names <- c("geneSymbol")

# Create a data frame with header names
df_sign_targets_crispr <- data.frame(geneSymbol = sign_targets_crispr)

## plot venn diagram
venn_diag<- ggVennDiagram(x=list(dex_comb_subset$geneSymbol, splice_func_df$geneSymbol,df_sign_targets_crispr$geneSymbol), 
                          edge_lty = "dashed", 
                          edge_size = 1,
                          label_size = 5,
                          set_size = 4,
                          category.names = c("DE" , "DS", "DG"),
                          label_percent_digit = 1) +  
  scale_fill_distiller(palette = "Blues", direction = 1, name = expression(bold("Gene count"))) + 
  labs(title = expression(bold("Differentially expressed and spliced genes and dependent genes"))) +
  coord_flip()

ggplot2::ggsave(venn_output_func_file,
                plot=venn_diag,
                width=5.5,
                height=4,
                device="pdf",
                dpi=300)


common <- intersect(unique(dex_comb_subset$geneSymbol), unique(splice_func_df$geneSymbol), unique(df_sign_targets_crispr$geneSymbol)) %>%
  sort() %>%
  write_lines(file.path(results_dir, "common_genes_de_ds_functional.txt"))

## run enrichR to compute and identify significant over-repr pathways
ora_results <- enricher(
  gene = intersect_func, # A vector of your genes of interest
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.1,
  pAdjustMethod = "BH", 
  TERM2GENE = dplyr::select(
    pathway_df,
    gs_name,
    human_gene_symbol
  )
)

ora_results@result$`B-H Adj p-value` <- ora_results@result$p.adjust

ora_result_df <- data.frame(ora_results@result)

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
  )+
  xlim(0,0.20)

ggplot2::ggsave(ora_dotplot_func_path,
                plot=enrich_plot_func,
                width=8,
                height=4,
                device="pdf",
                dpi=300)

