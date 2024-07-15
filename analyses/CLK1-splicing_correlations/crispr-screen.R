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

clin_df <- vroom(file.path(data_dir,"histologies.tsv")) %>%
  dplyr::filter(short_histology == "HGAT",
                cohort == "PBTA") 

clk1_targets <- read_lines(file.path(results_dir,"common_genes_de_ds_functional.txt"))

crispr <- read_csv(file.path(input_dir,"CCMA_crispr_genedependency_042024.csv")) %>%
  filter(wald_p_value < 0.05 & wald_fdr < 0.10 & z < 1.5) %>% 
  mutate(sample_id = str_replace(sample, "_[^_]*$", "")) %>%
  mutate(sample_id = str_replace_all(sample_id, "_", "-")) %>%
  inner_join(clin_df, by="sample_id") 




clk1_targets_crispr <- intersect(clk1_targets, crispr$gene)




crispr_filt <-crispr %>%
  filter(gene %in% intersect(clk1_targets, crispr$gene),
         grepl("1763|1769", sample))

unique(crispr_filt$gene) %>%
  write_lines("~/Downloads/clk1_targets_crispr_cbtn_lines.txt")