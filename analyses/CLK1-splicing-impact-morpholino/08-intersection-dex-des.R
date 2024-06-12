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

# construct a Venn object
# Define the lists of gene symbols
list_dex <- unique(dex_comb$geneSymbol)
list_spl <- unique(psi_comb$geneSymbol)

# Create the Venn diagram
venn.plot <- venn.diagram(
  x = list(
    DE = list_dex,
    DS = list_spl
  ),
  filename = NULL,  # Save to file as PDF
  height = 4, 
  width = 4,
  resolution = 300,
  
  # Customize the appearance
  col = "black",
  fill = c("lightblue", "blue"),
  alpha = c(0.1, 0.1),  # Increased transparency for each circle
  lty = "dashed",
  lwd = 2,
  cex = 2,
  fontface = "bold",
  cat.cex = 2,
  cat.fontface = "bold",
  cat.col = c("black", "black"),
  
  # Category label positions
  cat.pos = c(-20, 20),  # Adjust positions for horizontal layout
  cat.dist = c(0.05, 0.05),
  
  # Customize title and legend
  #main = expression(bold("Genes dysregulated by splicing and expression")),
  main.cex = 1.5,
  main.fontface = "bold",
  main.col = "black",
  label.cex = 1.5,
  label.percent.digits = 1
)

# Draw the Venn diagram
pdf(venn_output_file, height = 4, width = 4)
grid.draw(venn.plot)
dev.off()

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

ggplot2::ggsave(ora_dotplot_path,
                plot=enrich_plot_func,
                width=8.5,
                height=7,
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


# construct a Venn object
# Define the lists of gene symbols
list_dex_sub <- unique(dex_comb_subset$geneSymbol)
list_spl_sub <- unique(splice_func_df$geneSymbol)

# Create the Venn diagram
venn.plot <- venn.diagram(
  x = list(
    DE = list_dex_sub,
    DS = list_spl_sub
  ),
  filename = NULL,  # Save to file as PDF
  height = 4, 
  width = 4,
  resolution = 300,
  
  # Customize the appearance
  col = "black",
  fill = c("lightblue", "blue"),
  alpha = c(0.1, 0.1),  # Increased transparency for each circle
  lty = "dashed",
  lwd = 2,
  cex = 2,
  fontface = "bold",
  cat.cex = 2,
  cat.fontface = "bold",
  cat.col = c("black", "black"),
  
  # Category label positions
  cat.pos = c(-20, 20),  # Adjust positions for horizontal layout
  cat.dist = c(0.05, 0.05),
  
  # Customize title and legend
  #main = expression(bold("Genes dysregulated by splicing (functional) and expression")),
  main.cex = 1.5,
  main.fontface = "bold",
  main.col = "black",
  label.cex = 1.5,
  label.percent.digits = 1
)

# Draw the Venn diagram
pdf(venn_output_func_file, height = 4, width = 4)
grid.draw(venn.plot)
dev.off()


common <- intersect(unique(dex_comb_subset$geneSymbol), unique(splice_func_df$geneSymbol)) %>%
  sort() %>%
  write_lines(file.path(results_dir, "common_genes_de_ds_functional.txt"))

## run enrichR to compute and identify significant over-repr pathways
ora_results <- enricher(
  gene = intersect_func, # A vector of your genes of interest
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
                width=8.5,
                height=7,
                device="pdf",
                dpi=300)

